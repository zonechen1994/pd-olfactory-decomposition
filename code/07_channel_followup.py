"""07_channel_followup.py  Follow-up by recruitment cohort (Supplementary Table 2).

Rebuilds the survival frame from the source files with the same rules as 02_survival_cohort.py
(it must reproduce 1,759 / 152) and summarises enrolment year, follow-up, person-years and the
two-year Kaplan-Meier estimate for each recruitment cohort.
Output: results/repro/supp_channel_followup.csv
"""
from common import *
ex=pd.read_excel(CLIN_FILE,sheet_name=CLIN_SHEET)
ex['visit_date']=pd.to_datetime(ex['visit_date'],format='%m/%Y',errors='coerce'); assert ex['visit_date'].notna().all()
e4=ex[ex.COHORT==4].sort_values(['PATNO','visit_date'])
first=e4.groupby('PATNO')['visit_date'].min(); last=e4.groupby('PATNO')['visit_date'].max()
diag=e4[e4.PRIMDIAG==1].groupby('PATNO')['visit_date'].min()
prevalent=set(e4.groupby('PATNO').first().query('PRIMDIAG==1').index)
sub=e4.groupby('PATNO')['subgroup'].first()
sc=pd.read_csv('results/repro/all_mismatch_scores.csv'); d=sc[sc.COHORT==4].copy()
d['sub']=d.PATNO.map(sub).astype(str)
d['arm']=np.where(d['sub']=='Hyposmia','Hyposmia',np.where(d['sub'].str.contains('RBD'),'RBD','Genetic'))
d=d[~d.PATNO.isin(prevalent)]
d['conv']=d.PATNO.isin(diag.index).astype(int)
tconv=(d.PATNO.map(diag)-d.PATNO.map(first)).dt.days/365.25
tcens=(d.PATNO.map(last)-d.PATNO.map(first)).dt.days/365.25
d['t']=np.where(d.conv==1,np.maximum(tconv,0.25),tcens)
d=d[(d.conv==1)|(d.t>0)].dropna(subset=['OIS','OMI','t'])
d['enrol_year']=d.PATNO.map(first).dt.year
print('total',len(d),d.conv.sum())
assert len(d)==1759 and int(d.conv.sum())==152, 'survival frame does not reproduce 1,759 / 152'
rows=[]
for arm,g in list(d.groupby('arm'))+[('Non-hyposmia-enriched',d[d.arm!='Hyposmia']),('All prodromal',d)]:
    km=KaplanMeierFitter().fit(g.t,g.conv)
    rows.append(dict(channel=arm,n=len(g),events=int(g.conv.sum()),median_enrol_year=int(g.enrol_year.median()),
        median_followup_y=round(g.t.median(),2),pct_ge4y=round((g.t>=4).mean()*100,1),person_years=round(g.t.sum(),1),
        rate_per_100py=round(g.conv.sum()/g.t.sum()*100,2),km2y_pct=round((1-float(km.predict(2.0)))*100,1)))
out=pd.DataFrame(rows); print(out.to_string()); out.to_csv('results/repro/supp_channel_followup.csv',index=False)
