"""Numerical audit of the corrected analyses. Writes only under code/validated.

Uses the saved participant frames and re-fits the frozen ridge model as a check.
No outcome-guided model tuning. Canonical comparator: Coxnet, alpha=.1,
l1_ratio=.05, five folds (seed 41), scaler fitted inside each training fold.
Original primary bootstrap streams are reproduced; repeated estimates reuse
that result. All new paired comparisons use two-sided bootstrap tail areas,
with the documented 1/1000 reporting floor, not an exact-test interpretation.
"""
from pathlib import Path
import json, pickle, hashlib, platform, subprocess, shutil, importlib.metadata
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import cumulative_dynamic_auc
from lifelines.utils import concordance_index
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
R = OUT / 'results/repro'
R.mkdir(parents=True, exist_ok=True)
for p in (ROOT / 'results/repro').glob('*'):
    if p.is_file() and p.suffix in ['.csv', '.json']:
        shutil.copy2(p, R / p.name)
st = pickle.load(open(ROOT / 'results/intermediate/02_cohort.pkl', 'rb'))
J = json.load(open(ROOT / 'results/repro/section16_stats.json'))
SBR, DEMO = st['SBR_COLS'], st['DEMO']
df, train, prod = st['df'].copy(), st['train'].copy(), st['prod'].copy()
features = SBR + DEMO
tr = train.dropna(subset=features + ['upsit'])
sc = StandardScaler().fit(tr[features].values)
ridge = Ridge(alpha=1).fit(sc.transform(tr[features].values), tr.upsit.values)
valid = df.dropna(subset=features + ['OIS'])
err = np.max(np.abs(ridge.predict(sc.transform(valid[features].values)) - valid.OIS))
assert err < 1e-9, err
assert np.max(np.abs((df.OIS + df.OMI - df.upsit).dropna())) < 1e-10
print('Ridge predictions reproduced, maximum absolute difference', err, flush=True)

# Rebuild the clinical time variables from the saved raw clinical frame.
e4 = st['excel'][st['excel'].COHORT == 4].sort_values(['PATNO', 'visit_date'])
first = e4.groupby('PATNO').visit_date.min()
last = e4.groupby('PATNO').visit_date.max()
diagnosis = e4[e4.PRIMDIAG == 1].groupby('PATNO').visit_date.min()
prevalent = set(e4.groupby('PATNO').first().query('PRIMDIAG == 1').index)
base = df[df.COHORT == 4].copy()
base['converted'] = base.PATNO.isin(diagnosis.index).astype(int)
base['time_years'] = np.where(base.converted == 1,
    np.maximum((base.PATNO.map(diagnosis)-base.PATNO.map(first)).dt.days/365.25, .25),
    (base.PATNO.map(last)-base.PATNO.map(first)).dt.days/365.25)
needs = ['OIS','OMI','upsit','PUTAMEN_REF_CWM','PUTAMEN_L_REF_CWM',
         'PUTAMEN_R_REF_CWM','time_years','converted','age','sex_num']
flow = []
def record(label, g):
    flow.append(dict(stage=label,n=len(g),events=int(g.converted.sum()),
        prevalent=int(g.PATNO.isin(prevalent).sum()),
        incident=int(g.loc[~g.PATNO.isin(prevalent),'converted'].sum())))
for label, b in [('all',base),('hyposmia',base[base.subgroup=='Hyposmia'])]:
    record(label+'_baseline', b)
    b=b[~b.PATNO.isin(prevalent)]; record(label+'_exclude_prevalent',b)
    b=b[(b.converted==1)|(b.time_years>0)]; record(label+'_with_followup',b)
    b=b.dropna(subset=needs); record(label+'_complete',b)
pw=base[(base.converted==1)|(base.time_years>0)].copy()
record('retaining_prevalent_before_complete_case',pw)
pw=pw.dropna(subset=needs);record('retaining_prevalent_complete_case',pw)
pd.DataFrame(flow).to_csv(R/'audit_flow.csv',index=False)
d=base[~base.PATNO.isin(prevalent)].copy()
d=d[(d.converted==1)|(d.time_years>0)].dropna(subset=needs)
assert set(d.PATNO)==set(prod.dropna(subset=needs).PATNO)
assert np.allclose(d.set_index('PATNO').time_years.sort_index(),
                   prod.dropna(subset=needs).set_index('PATNO').time_years.sort_index())
# Preserve the original row order for exact bootstrap replication.
d=prod.dropna(subset=needs).copy()
d['pct_exp']=d[['PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM']].min(axis=1)/(1.8547-.00303*d.age-.2419*d.sex_num)*100
h=d[d.subgroup=='Hyposmia'].copy()
groups={'all prodromal':d,'hyposmia':h,
        'hyposmia + DAT deficit':h[h.pct_exp<65].copy(),
        'hyposmia + no DAT deficit':h[h.pct_exp>=65].copy(),
        'non-hyposmia-enriched':d[d.subgroup!='Hyposmia'].copy()}
assert (len(d),int(d.converted.sum()),len(h),int(h.converted.sum()))==(1759,152,1003,68)
assert not d[features].isna().any().any()

def c(g,v):
    return float(concordance_index(g.time_years,np.asarray(v),g.converted))
def boot(g,a,b=None,rng=None):
    rng=np.random.default_rng(42) if rng is None else rng
    a=np.asarray(a); b=None if b is None else np.asarray(b)
    assert np.isfinite(a).all() and (b is None or np.isfinite(b).all())
    t,e=g.time_years.to_numpy(),g.converted.to_numpy()
    vals=[]
    for _ in range(1000):
        i=rng.choice(len(g),len(g),replace=True)
        if e[i].sum()<5: continue
        v=concordance_index(t[i],a[i],e[i])
        if b is not None: v-=concordance_index(t[i],b[i],e[i])
        vals.append(v)
    vals=np.asarray(vals); lo,hi=np.percentile(vals,[2.5,97.5])
    z=dict(n=len(g),events=int(e.sum()),c=c(g,a),lo=float(lo),hi=float(hi),
           bootstrap_requested=1000,bootstrap_valid=len(vals))
    if b is not None:
        le=int((vals<=0).sum()); ge=int((vals>=0).sum())
        raw=2*min(le,ge)/len(vals)
        z.update(c_a=z.pop('c'),c_b=c(g,b),dC=c(g,a)-c(g,b),
                 p_raw=raw,p=max(raw,.001),n_le_zero=le,n_ge_zero=ge,
                 p_formula='max(2*min(n_le_zero,n_ge_zero)/B_valid,1/1000)',
                 alternative='two-sided')
    else:z['at_chance']=bool(lo<=.5<=hi)
    return z

# Reproduce the original six primary comparisons, including their RNG stream.
plan=[('hyposmia','upsit','UPSIT'),('hyposmia','PUTAMEN_REF_CWM','putamen'),
      ('hyposmia + DAT deficit','upsit','UPSIT'),
      ('hyposmia + no DAT deficit','PUTAMEN_REF_CWM','putamen'),
      ('hyposmia + no DAT deficit','upsit','UPSIT'),
      ('non-hyposmia-enriched','upsit','UPSIT')]
rng=np.random.default_rng(42); paired=[]
for name,col,label in plan:
    g=groups[name];paired.append(dict(stratum=name,comparison='OIS - '+label,**boot(g,g.OIS,g[col],rng)))
old=pd.read_csv(ROOT/'results/repro/section33_prespecified.csv')
new=pd.DataFrame(paired)
for col in ['c_a','c_b','dC','lo','hi','p']: assert np.allclose(old[col],new[col],atol=1e-12)
new.to_csv(R/'section33_prespecified.csv',index=False);J['prespecified_comparisons']=paired
print('All six primary comparisons reproduced exactly',flush=True)

# Single C-index intervals: exact original section 36 stream, reused everywhere.
readings=[('upsit','UPSIT total'),('OIS','OIS'),('OMI','OMI'),
          ('PUTAMEN_REF_CWM','putamen SBR'),('pct_exp','lowest putamen %expected')]
rng=np.random.default_rng(42); singles=[]
for name,g in groups.items():
    for col,label in readings: singles.append(dict(stratum=name,reading=label,**boot(g,g[col],rng=rng)))
old=pd.read_csv(ROOT/'results/repro/section36_cindex_ci.csv');new=pd.DataFrame(singles)
for col in ['c','lo','hi']: assert np.allclose(old[col],new[col],atol=1e-12)
new.to_csv(R/'section36_cindex_ci.csv',index=False);J['cindex_ci']=singles
lookup={(r['stratum'],r['reading']):r for r in singles}
print('All section 36 estimates and intervals reproduced exactly',flush=True)

# Align repeated continuous scores; binary readings keep their original 1000-resample intervals.
c38=pd.read_csv(R/'section38_continuous_vs_binary.csv')
print('section38 labels',c38.reading.tolist(),flush=True)
for i,r in c38.iterrows():
    lab={'UPSIT, continuous score':'UPSIT total','imaging, continuous (lowest putamen %exp)':'lowest putamen %expected',
         'OIS, continuous':'OIS'}.get(r.reading)
    if lab and (r.stratum,lab) in lookup:
        for col in ['c','lo','hi']:c38.loc[i,col]=lookup[r.stratum,lab][col]
c38.to_csv(R/'section38_continuous_vs_binary.csv',index=False)
J['continuous_vs_binary']=c38.to_dict('records')

# One implementation for all conversion-supervised comparisons.
folds=[]
def cvcox(g,cols,label):
    X=g[cols].to_numpy(float);y=Surv.from_arrays(g.converted.astype(bool),g.time_years)
    pred=np.full(len(g),np.nan)
    for k,(tri,tei) in enumerate(KFold(5,shuffle=True,random_state=41).split(X)):
        scaler=StandardScaler().fit(X[tri])
        m=CoxnetSurvivalAnalysis(l1_ratio=.05,alphas=[.1],fit_baseline_model=False)
        m.fit(scaler.transform(X[tri]),y[tri]);pred[tei]=-m.predict(scaler.transform(X[tei]))
        folds.append(dict(model=label,fold=k+1,n_train=len(tri),n_test=len(tei),
                          events_train=int(g.converted.iloc[tri].sum()),events_test=int(g.converted.iloc[tei].sum()),status='ok'))
    assert np.isfinite(pred).all()
    return pred
scs=StandardScaler().fit(tr[SBR].values)
rs=Ridge(alpha=1).fit(scs.transform(tr[SBR].values),tr.upsit.values)
cal=[];predictions={};rng=np.random.default_rng(42)
for name in ['hyposmia','hyposmia + DAT deficit','hyposmia + no DAT deficit']:
    g=groups[name];g['OIS_sbr']=rs.predict(scs.transform(g[SBR].values))
    vv={'single-region putamen':g.PUTAMEN_REF_CWM,'lowest putamen %expected':g.pct_exp,
        'age alone':-g[DEMO[0]],'demographics only, CV Cox':cvcox(g,DEMO,name+'/demographics'),
        '33 SBR, CV Cox on conversion':cvcox(g,SBR,name+'/SBR'),
        '33 SBR + demographics, CV Cox':cvcox(g,SBR+DEMO,name+'/SBR+demographics'),
        'OIS, imaging only, zero-shot':g.OIS_sbr,'OIS, full, zero-shot':g.OIS,'UPSIT total':g.upsit}
    predictions[name]=vv
    maps={'single-region putamen':'putamen SBR','lowest putamen %expected':'lowest putamen %expected',
          'OIS, full, zero-shot':'OIS','UPSIT total':'UPSIT total'}
    for label,v in vv.items():
        z=dict(lookup[name,maps[label]]) if label in maps else boot(g,v,rng=rng)
        z.update(stratum=name,reading=label)
        z['kind']='supervised in-cohort' if 'CV Cox' in label else 'score'
        cal.append(z)
for r in J['calibration_and_age']:
    if r.get('kind')=='multivariable':cal.append(r)
pd.DataFrame(cal).to_csv(R/'section40_calibration_and_age.csv',index=False);J['calibration_and_age']=cal
head=[]
hv=predictions['hyposmia']
hv['Putamen + age/sex/edu']=cvcox(h,['PUTAMEN_REF_CWM']+DEMO,'hyposmia/putamen+demographics')
hv['Binary DAT-deficit flag']=-(h.pct_exp<65).astype(float)
for label,v in [('UPSIT (raw)',h.upsit),('Binary DAT-deficit flag',hv['Binary DAT-deficit flag']),
                ('Putamen + age/sex/edu',hv['Putamen + age/sex/edu']),
                ('33-SBR only',hv['33 SBR, CV Cox on conversion']),
                ('33-SBR + age/sex/edu',hv['33 SBR + demographics, CV Cox']),
                ('OIS, SBR-only',h.OIS_sbr)]:
    z=dict(paired[0]) if label=='UPSIT (raw)' else boot(h,h.OIS,v)
    z.update(model=label,c=c(h,v),dc=z['dC'])
    z['training']='Coxnet 5-fold CV (in-cohort)' if label.startswith(('33-','Putamen')) else 'fixed reading'
    head.append(z)
head.append(dict(model='OIS, full',c=c(h,h.OIS),dc=np.nan,lo=np.nan,hi=np.nan,p=np.nan,training='Ridge on UPSIT, PD+HC'))
J['head_to_head']=head;pd.DataFrame(head).to_csv(R/'audit_head_to_head.csv',index=False)
pd.DataFrame(folds).to_csv(R/'audit_cv_folds.csv',index=False)
assert np.isclose(head[4]['c'],.8124924960979709)
print('Canonical supervised models fitted, all folds valid',flush=True)

# True landmark: alive/event-free and observed beyond L, outcomes after L.
lands=[]
for L in [0,1,2]:
    g=h[h.time_years>L].copy();g['time_years']-=L
    z=dict(paired[0]) if L==0 else boot(g,g.OIS,g.upsit)
    z.update(landmark_yr=L,c_ois=z['c_a'],c_upsit=z['c_b'],dc=z['dC'],
             excluded_before_or_at_landmark=len(h)-len(g),time_origin=f'{L} years after baseline')
    lands.append(z)
J['landmark']=lands;pd.DataFrame(lands).to_csv(R/'audit_landmark.csv',index=False)
print('Landmarks',[(r['landmark_yr'],r['n'],r['events'],r['c_ois'],r['c_upsit'],r['p']) for r in lands],flush=True)

# Recompute KM median splits independently.
median=[]
for name in ['hyposmia','hyposmia + no DAT deficit']:
    g=groups[name]
    for col,label in [('OIS','OIS'),('upsit','UPSIT total'),('pct_exp','lower putamen %expected')]:
        cut=g[col].median();a=g[g[col]<=cut];b=g[g[col]>cut]
        z=dict(stratum=name,reading=label,median=cut)
        for suffix,x in [('low',a),('high',b)]:
            km=KaplanMeierFitter().fit(x.time_years,x.converted)
            ci=km.confidence_interval_survival_function_;ci=ci.loc[ci.index<=2].iloc[-1]
            z.update({f'n_{suffix}':len(x),f'ev_{suffix}':int(x.converted.sum()),
                      f'km2_{suffix}':100*(1-km.predict(2)),f'km2_{suffix}_lo':100*(1-ci.iloc[1]),f'km2_{suffix}_hi':100*(1-ci.iloc[0])})
        z['p_logrank']=logrank_test(a.time_years,b.time_years,a.converted,b.converted).p_value
        median.append(z)
pd.DataFrame(median).to_csv(R/'supp_median_split.csv',index=False)
old=pd.read_csv(ROOT/'results/repro/supp_median_split.csv');new=pd.DataFrame(median)
for col in ['median','n_low','ev_low','km2_low','n_high','ev_high','km2_high','p_logrank']:
    assert np.allclose(old[col],new[col],atol=1e-10)

# Fixed-horizon point estimates and denominators, independent call to IPCW estimator.
y=Surv.from_arrays(h.converted.astype(bool),h.time_years)
aucs={col:cumulative_dynamic_auc(y,y,-h[col].to_numpy(),[1,1.5,2])[0].tolist() for col in ['OIS','upsit','PUTAMEN_REF_CWM']}
assert np.isclose(aucs['OIS'][2],.8610868386818654)
assert np.isclose(aucs['upsit'][2],.7336376270368187)
assert int((h.time_years>=2).sum())==432 and int(((h.time_years<=2)&(h.converted==1)).sum())==41
J['audit_20260924']={'flow':flow,'ridge_max_abs_error':float(err),'ipcw_point_estimates':aucs,
 'scope':'Recomputed survival comparisons, landmark, KM, frozen model; biological results checked against archived outputs separately.',
 'retaining_prevalent_complete_case':dict(n=len(pw),events=int(pw.converted.sum()),
    prevalent=int(pw.PATNO.isin(prevalent).sum()),incident=int(pw.loc[~pw.PATNO.isin(prevalent),'converted'].sum()),
    c_ois=c(pw,pw.OIS),c_upsit=c(pw,pw.upsit))}
with open(R/'section16_stats.json','w') as f:json.dump(J,f,indent=2,default=float)
files=['results/intermediate/02_cohort.pkl','code/01_fit_model.py','code/02_survival_cohort.py',
       'Xing_Core_Lab_-_Quant_SBR_18Mar2026.csv','PPMI_Curated_Data_Cut_Public_20260223.xlsx']
manifest={'audit_completed_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','scikit-survival','lifelines','statsmodels']},
 'release_commit':subprocess.check_output(['git','-C',str(ROOT/'release/pd-olfactory-decomposition'),'rev-parse','HEAD'],text=True).strip(),
 'inputs':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('AUDIT COMPLETE',flush=True)
