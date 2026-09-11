"""06_ph_test.py  Proportional-hazards check for the multivariable Cox models reported in the main text.

Rebuilds the survival frame from the source files with the same rules as 02_survival_cohort.py (prevalent baseline PD
excluded, follow-up from first visit to first PD diagnosis or last visit, 0.25-year floor for
converters only), reads OIS/OMI from the ledger export, and runs lifelines' Schoenfeld-residual
proportional_hazard_test on each reported model. Writes results/repro/supp_ph_test.csv and adds
the key `ph_test` to section16_stats.json.
"""
from common import *
from lifelines.statistics import proportional_hazard_test

excel = pd.read_excel(CLIN_FILE, sheet_name=CLIN_SHEET)
excel['visit_date'] = pd.to_datetime(excel['visit_date'], format='%m/%Y', errors='coerce')
assert excel['visit_date'].notna().all()
sbr = pd.read_csv(SBR_FILE)
sbr['visit_rank'] = sbr['EVENT_ID'].map({'SC': 0, 'BL': 1, 'V04': 2, 'V06': 3}).fillna(99)
sbr_bl = sbr.sort_values(['PATNO', 'visit_rank']).drop_duplicates('PATNO', keep='first')
SBR_COLS = [c for c in sbr_bl.columns if 'REF_CWM' in c]
sbr_bl = sbr_bl.dropna(subset=SBR_COLS)
clin = excel[excel['EVENT_ID'] == 'BL'].drop_duplicates('PATNO')
df = sbr_bl.merge(clin[['PATNO', 'COHORT', 'age', 'SEX', 'EDUCYRS', 'upsit', 'subgroup']], on='PATNO', how='inner')
df['sex_num'] = df['SEX'].astype(float); df['educyrs'] = df['EDUCYRS'].astype(float)
scores = pd.read_csv(f'results/repro/all_mismatch_scores.csv')[['PATNO', 'OIS', 'OMI']]
df = df.merge(scores, on='PATNO', how='left')

# survival frame (identical rules to the notebook)
e4 = excel[excel['COHORT'] == 4].sort_values(['PATNO', 'visit_date'])
first_visit = e4.groupby('PATNO')['visit_date'].min(); last_visit = e4.groupby('PATNO')['visit_date'].max()
pd_diag = e4[e4['PRIMDIAG'] == 1].groupby('PATNO')['visit_date'].min()
prevalent = set(e4.groupby('PATNO').first().query('PRIMDIAG == 1').index)
prod = df[(df['COHORT'] == 4) & ~df['PATNO'].isin(prevalent)].copy()
prod['converted'] = prod['PATNO'].isin(pd_diag.index).astype(int)
fv = pd.to_datetime(prod['PATNO'].map(first_visit)); lv = pd.to_datetime(prod['PATNO'].map(last_visit)); pdd = pd.to_datetime(prod['PATNO'].map(pd_diag))
tconv = (pdd - fv).dt.days / 365.25; tcens = (lv - fv).dt.days / 365.25
prod['time_years'] = np.where(prod['converted'] == 1, np.maximum(tconv, 0.25), tcens)
prod = prod[(prod['converted'] == 1) | (prod['time_years'] > 0)].copy()
prod = prod.dropna(subset=['OIS', 'upsit', 'age', 'sex_num', 'time_years'])
hyp = prod[prod['subgroup'] == 'Hyposmia'].copy()
# imaging stratum (PPMI canonical metric, HC-normative equation from the ledger)
J = json.load(open(f'results/repro/section16_stats.json'))
D = J['dat_metric_validation']
exp_ = D['norm_intercept'] + D['norm_beta_age'] * hyp['age'].astype(float) + D['norm_beta_sex'] * hyp['sex_num']
hyp['pct_exp'] = np.minimum(hyp['PUTAMEN_L_REF_CWM'], hyp['PUTAMEN_R_REF_CWM']) / exp_ * 100
print(f'hyposmia group n={len(hyp)} events={int(hyp.converted.sum())} | non-deficit n={int((hyp.pct_exp>=65).sum())} events={int(hyp[hyp.pct_exp>=65].converted.sum())}')
assert len(hyp) == 1003 and int(hyp.converted.sum()) == 68, 'survival frame does not reproduce 1,003 / 68'

rows = []
def run(d, label, score, covs, scale):
    x = d[['time_years', 'converted', score] + covs].dropna().copy()
    if scale == 'z':
        for c in [score] + [c for c in covs if c != 'sex_num']:
            x[c] = (x[c] - x[c].mean()) / x[c].std(ddof=1)
    cph = CoxPHFitter().fit(x, 'time_years', 'converted')
    ph = proportional_hazard_test(cph, x, time_transform='rank'); s = ph.summary
    model = f"{score} + {' + '.join(covs)} ({'per point / year' if scale=='raw' else 'z-scored'})"
    for term in s.index:
        rows.append(dict(population=label, model=model, term=term, hr=float(np.exp(cph.params_[term])),
                         p_cox=float(cph.summary.loc[term, 'p']), model_c=float(cph.concordance_index_),
                         ph_test_stat=float(s.loc[term, 'test_statistic']), ph_p=float(s.loc[term, 'p']),
                         n=int(len(x)), events=int(x['converted'].sum())))
    gstat = float(s['test_statistic'].sum())
    rows.append(dict(population=label, model=model, term='global', hr=np.nan, p_cox=np.nan, model_c=float(cph.concordance_index_),
                     ph_test_stat=gstat, ph_p=float(1 - __import__('scipy').stats.chi2.cdf(gstat, df=len(s))),
                     n=int(len(x)), events=int(x['converted'].sum())))

for score in ['OIS', 'upsit']:
    for pop, d in [('hyposmia group', hyp), ('non-deficit stratum', hyp[hyp['pct_exp'] >= 65])]:
        run(d, pop, score, ['age'], 'raw')                         # the model reported in the main text
        run(d, pop, score, ['age', 'sex_num', 'educyrs'], 'z')    # fuller specification
out = pd.DataFrame(rows)
out.to_csv(f'results/repro/supp_ph_test.csv', index=False)
pd.set_option('display.width', 200)
print(out.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
J['ph_test'] = dict(method="Schoenfeld residuals, lifelines proportional_hazard_test, time_transform='rank'",
                    rows=out.to_dict('records'))
json.dump(J, open(f'results/repro/section16_stats.json', 'w'), indent=2, default=float)
print('-> results/repro/supp_ph_test.csv ; ledger key ph_test')
