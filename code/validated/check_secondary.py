"""Read-only checks of source snapshots and residual analyses; aggregate output only."""
from pathlib import Path
import os,pickle,json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[2];DEST=Path(__file__).resolve().parent
BASE=Path(os.environ.get('PPMI_DATA_DIR', ROOT))
st=pickle.load(open(ROOT/'results/intermediate/02_cohort.pkl','rb'))
df,excel=st['df'],st['excel']
raw=pd.read_excel(BASE/'PPMI_Curated_Data_Cut_Public_20260223.xlsx',sheet_name='20260105')
raw['visit_date']=pd.to_datetime(raw.visit_date,format='%m/%Y',errors='coerce')
pd.testing.assert_frame_equal(raw,excel)
sbr=pd.read_csv(BASE/'Xing_Core_Lab_-_Quant_SBR_18Mar2026.csv')
sbr['visit_rank']=sbr.EVENT_ID.map({'SC':0,'BL':1,'V04':2,'V06':3}).fillna(99)
pd.testing.assert_frame_equal(sbr,st['sbr'])
print('Raw clinical and imaging snapshots exactly match saved state',flush=True)
summary={'raw_snapshots_match_state':True,'missing_covariates':{}}
for label,g in [('training',st['train'].dropna(subset=['upsit'])),
                ('prodromal_baseline',df[df.COHORT==4]),
                ('survival',st['prod'].dropna(subset=['OIS','upsit','age','sex_num']))]:
    summary['missing_covariates'][label]={k:int(v) for k,v in g[['age','SEX','EDUCYRS']].isna().sum().items()}
src=(ROOT/'code/02_survival_cohort.py').read_text()
part=src[src.index("ois_map = df.set_index('PATNO')['OIS']"):src.index("pd.DataFrame(_icc12).to_csv")]
exec(compile(part,'ICC_check','exec'))
summary['icc']=ICC_STATS
J=json.load(open(ROOT/'results/repro/section16_stats.json'))
assert np.isclose(ICC_STATS['omi_icc'],J['icc']['omi_icc']) and ICC_STATS['n_pairs']==1809
BIO_FILE=BASE/'molecular/Current_Biospecimen_Analysis_Results_24Mar2026.csv'
src=(ROOT/'code/04_residual.py').read_text()
part=src[src.index("if 'PROD22' not in globals():"):src.index("print('(a) Three-way table")]
exec(compile(part,'CSF_check','exec'))
old=pd.read_csv(ROOT/'results/repro/section22_csf_threeway.csv')
for col in ['n','b_UPSIT','b_OIS','b_OMI','p_UPSIT','p_OIS','p_OMI']:
    assert np.allclose(old[col],CSF22[col],equal_nan=True)
CSF22.to_csv(DEST/'results/repro/audit_csf_recomputed.csv',index=False)
summary['csf_all_rows_reproduced']=True
summary['tau_ratio_prodromal']=CSF22[(CSF22.cohort=='Prodromal')&(CSF22.key=='pTau181_over_ABeta42_CSF')].to_dict('records')
part=src[src.index('gen20 = {}'):src.index('print(\'  Reading: in prodromals')]
sec16={};exec(compile(part,'genotype_check','exec'))
for a,b in zip(sec16['genotype'],J['genotype']):
    for col in ['beta','lo','hi','p','n_gba','n_lrrk2']:assert np.isclose(a[col],b[col])
summary['genotype_all_rows_reproduced']=True
part=src[src.index('from scipy import stats as _st42'):src.index('# --- (b) manifest PD, GBA vs LRRK2')]
exec(compile(part,'partial_r_check','exec'))
summary['partial_correlations']=rows42
(DEST/'secondary_checks.json').write_text(json.dumps(summary,indent=2,default=float))
print('SECONDARY CHECKS COMPLETE',flush=True)
