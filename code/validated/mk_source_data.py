"""mk_source_data.py  Source data for every figure panel, one workbook.

Writes Source_Data.xlsx with one sheet per panel of Figs 1d, 2, 3 and 4 and of Supplementary
Figs 1 to 5. Every sheet holds exactly the values that the corresponding panel draws, read from
the ledger and the section tables that the figure scripts read, or recomputed on the same
participant-level frame with the same estimators (Kaplan-Meier tables, decile means).
Participant-level points (the scatter of Fig. 2a and the faint points of Fig. 2c) are not
included because the PPMI Data Use Agreement does not permit their redistribution; for those
panels the summary statistics and the decile means are given instead.
"""
import os, sys, json, pickle
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline

ROOT = os.path.dirname(os.path.abspath(__file__)); os.chdir(ROOT)
R = 'results/repro'
J = json.load(open(f'{R}/section16_stats.json'))
STATE = 'results/intermediate/02_cohort.pkl'
if not os.path.exists(STATE):
    sys.exit('participant-level frames not found, run 02_survival_cohort.py first')
st = pickle.load(open(STATE, 'rb')); df = st['df']; prod = st['prod'].copy(); FEAT = st['FEATURES']; SBR = st['SBR_COLS']
OUTFILE = 'Source_Data.xlsx'
sheets = {}      # name -> (DataFrame, note)


def add(name, frame, note=''):
    sheets[name] = (frame.reset_index(drop=True), note)


def r_ci(r, n, z=1.96):
    zf = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return np.tanh(zf - z * se), np.tanh(zf + z * se)


def km_table(g, label, tmax=2.0):
    """Kaplan-Meier survival table to tmax with Greenwood 95% CI, as drawn (step 'post')."""
    k = KaplanMeierFitter().fit(g.time_years, g.converted)
    sf = k.survival_function_; ci = k.confidence_interval_survival_function_
    m = sf.index <= tmax
    at_risk = [int((g.time_years >= t).sum()) for t in sf.index[m]]
    out = pd.DataFrame({'group': label, 'time_years': sf.index[m], 'free_of_PD': sf.iloc[:, 0][m].values,
                        'ci_lower': ci.iloc[:, 0][m].values, 'ci_upper': ci.iloc[:, 1][m].values, 'at_risk': at_risk})
    idx = ci.index[ci.index <= tmax][-1]
    two = dict(group=label, n=int(len(g)), events=int(g.converted.sum()),
               two_year_conversion_pct=(1 - float(k.predict(tmax))) * 100,
               two_year_ci_lower_pct=(1 - ci.loc[idx].iloc[1]) * 100, two_year_ci_upper_pct=(1 - ci.loc[idx].iloc[0]) * 100)
    return out, two


def hyposmia_frame():
    hh = prod[prod.subgroup == 'Hyposmia'].dropna(subset=['OIS', 'upsit', 'PUTAMEN_REF_CWM', 'PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM',
                                                          'time_years', 'converted', 'age', 'sex_num'] + SBR).copy()
    hh = hh[(hh.converted == 1) | (hh.time_years > 0)]
    hh['pct'] = hh[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (1.8547 - .00303 * hh.age - .2419 * hh.sex_num) * 100
    return hh


# ================================================================== Fig 1d
T1 = pd.DataFrame(J['tau_gradient_fig1d']); T1 = T1[T1.split == 'tertile'].copy()
T1['group'] = T1['group'].map({'G1': 'low', 'G2': 'middle', 'G3': 'high'})
add('Fig1d_middle', T1[['score', 'group', 'n', 'tau_mean', 'adj_mean', 'adj_se', 'p_across_groups']],
    'Adjusted means (age, sex, years of education) of OMI and OIS across tertiles of CSF pTau181/Abeta42, prodromal cohort. Bars are adj_mean, error bars adj_se.')
P1 = pd.DataFrame(J['partial_r_two_tests'])
add('Fig1d_right', P1[['test', 'cohort', 'score', 'n', 'r', 'abs_r', 'p']],
    'Partial correlations adjusted for age, sex and years of education. Bars are abs_r.')

# ================================================================== Fig 2
pool = df[df.COHORT.isin([1, 2])].dropna(subset=FEAT + ['upsit']).copy()
X, y = pool[FEAT].values, pool.upsit.values
cv = cross_val_predict(Ridge(alpha=1.0), StandardScaler().fit_transform(X), y, cv=KFold(5, shuffle=True, random_state=41))
r2a, p2a = stats.pearsonr(cv, y); lo2a, hi2a = r_ci(r2a, len(y)); sl, ic = np.polyfit(cv, y, 1)
ispd = (pool.COHORT == 1).values
dec = pd.qcut(cv, 10, labels=False)
f2a = pd.DataFrame({'decile_of_cross_validated_OIS': np.arange(1, 11),
                    'n': pd.Series(y).groupby(dec).size().values,
                    'mean_cross_validated_OIS': pd.Series(cv).groupby(dec).mean().values,
                    'mean_observed_UPSIT': pd.Series(y).groupby(dec).mean().values,
                    'sd_observed_UPSIT': pd.Series(y).groupby(dec).std().values})
f2a_stats = pd.DataFrame([dict(statistic='n', value=len(y)), dict(statistic='n_PD', value=int(ispd.sum())), dict(statistic='n_HC', value=int((~ispd).sum())),
                          dict(statistic='pearson_r', value=r2a), dict(statistic='r_ci_lower', value=lo2a), dict(statistic='r_ci_upper', value=hi2a),
                          dict(statistic='p', value=p2a), dict(statistic='regression_slope', value=sl), dict(statistic='regression_intercept', value=ic)])
add('Fig2a', pd.concat([f2a_stats, pd.DataFrame([{}]), f2a], axis=0),
    'Five-fold cross-validated OIS against observed UPSIT in the training cohort. Individual points are participant-level PPMI data and are not redistributed; decile means are given.')
loso = pd.read_csv(f'{R}/loso_per_site.csv').sort_values('r').reset_index(drop=True)
loso.insert(0, 'bar_from_bottom', np.arange(1, len(loso) + 1))
loso['site'] = 'site_' + loso['site'].astype(int).astype(str)
add('Fig2b', loso, f'One bar per held-out centre, sorted by r. Mean {loso.r.mean():.3f}, SD {loso.r.std(ddof=0):.3f}, median {loso.r.median():.3f}.')
pr = df[df.COHORT == 4].dropna(subset=['OIS', 'upsit']).copy()
pr['arm'] = np.where(pr.subgroup == 'Hyposmia', 'Hyposmia', np.where(pr.subgroup.astype(str).str.contains('RBD'), 'RBD', 'Variant carriers'))
rows, decs = [], []
for a_, g in pr.groupby('arm'):
    r, pc = stats.pearsonr(g.OIS, g.upsit); s2, i2 = np.polyfit(g.OIS, g.upsit, 1)
    rows.append(dict(recruitment_cohort=a_, n=len(g), slope=s2, intercept=i2, r=r, p=pc))
    d_ = pd.qcut(g.OIS, 10, labels=False, duplicates='drop')
    decs.append(pd.DataFrame({'recruitment_cohort': a_, 'decile': np.arange(1, d_.max() + 2), 'n': g.groupby(d_).size().values,
                              'mean_OIS': g.groupby(d_).OIS.mean().values, 'mean_observed_UPSIT': g.groupby(d_).upsit.mean().values}))
rows.append(dict(recruitment_cohort='Pooled', n=len(pr), slope=np.nan, intercept=np.nan, r=stats.pearsonr(pr.OIS, pr.upsit)[0], p=stats.pearsonr(pr.OIS, pr.upsit)[1]))
add('Fig2c', pd.concat([pd.DataFrame(rows), pd.DataFrame([{}]), pd.concat(decs)], axis=0),
    'Calibration by recruitment cohort: regression lines (slope, intercept), correlations and decile means (filled circles). Faint individual points are not redistributed.')

# ================================================================== Fig 3
H = [d for d in J['decomposition'] if d['cohort'] == 'Hyposmia'][0]
C36 = {r_['reading']: r_ for r_ in J['cindex_ci'] if r_['stratum'] == 'hyposmia'}
f3a = pd.DataFrame([dict(score='UPSIT', c_index=C36['UPSIT total']['c'], ci_lower=C36['UPSIT total']['lo'], ci_upper=C36['UPSIT total']['hi'], variance_share_pct=100.0),
                    dict(score='OIS', c_index=C36['OIS']['c'], ci_lower=C36['OIS']['lo'], ci_upper=C36['OIS']['hi'], variance_share_pct=H['sd_ois'] ** 2 / H['sd_upsit'] ** 2 * 100),
                    dict(score='OMI', c_index=C36['OMI']['c'], ci_lower=C36['OMI']['lo'], ci_upper=C36['OMI']['hi'], variance_share_pct=H['sd_omi'] ** 2 / H['sd_upsit'] ** 2 * 100)])
pu = [r_ for r_ in J['prespecified_comparisons'] if r_['stratum'] == 'hyposmia' and r_['comparison'] == 'OIS - UPSIT'][0]
f3a_b = pd.DataFrame([dict(item='n', value=H['n']), dict(item='events', value=68), dict(item='sd_UPSIT', value=H['sd_upsit']), dict(item='sd_OIS', value=H['sd_ois']),
                      dict(item='sd_OMI', value=H['sd_omi']), dict(item='r_OIS_OMI', value=H['r_ois_omi']),
                      dict(item='paired_delta_C_OIS_minus_UPSIT', value=pu['dC']), dict(item='delta_C_ci_lower', value=pu['lo']), dict(item='delta_C_ci_upper', value=pu['hi']), dict(item='delta_C_p', value=pu['p'])])
add('Fig3a', pd.concat([f3a, pd.DataFrame([{}]), f3a_b], axis=0), 'Concordance with 95% bootstrap CI and share of UPSIT variance, hyposmia group.')
td = pd.DataFrame(J['td_auc_ipcw']['rows'])
add('Fig3b', td[td.horizon_yr.isin(J['td_auc_ipcw']['times'])][['horizon_yr', 'reading', 'auc', 'lo', 'hi', 'n_at_risk', 'events_by_t', 'p']],
    'Cumulative/dynamic AUC with IPCW, 95% CI from 1,000 bootstrap resamples. The row OIS minus UPSIT carries the paired P.')
CA = pd.DataFrame([r for r in J['calibration_and_age'] if 'c' in r and r['c'] == r['c']])
METH = ['single-region putamen', 'lowest putamen %expected', '33 SBR + demographics, CV Cox', 'age alone', 'UPSIT total', 'OIS, imaging only, zero-shot', 'OIS, full, zero-shot']
LAB = dict(zip(METH, ['Putamen SBR', 'Lower putamen, % expected', '33 SBR + demographics, supervised Cox', 'Age', 'UPSIT', 'OIS, imaging only', 'OIS']))
ST = ['hyposmia', 'hyposmia + DAT deficit', 'hyposmia + no DAT deficit']
f3c = CA[CA.reading.isin(METH) & CA.stratum.isin(ST)][['stratum', 'reading', 'n', 'events', 'c', 'lo', 'hi']].copy()
f3c['reading_label'] = f3c.reading.map(LAB); f3c['interval_includes_0.5'] = f3c.lo < 0.5
f3c_p = pd.DataFrame(J['prespecified_comparisons'])[['stratum', 'comparison', 'n', 'events', 'dC', 'lo', 'hi', 'p']]
add('Fig3c', pd.concat([f3c, pd.DataFrame([{}]), f3c_p], axis=0), 'Seven readings in three imaging strata (bars with 95% CI), followed by the paired comparisons printed beneath the groups.')
hh = hyposmia_frame(); ndf = hh[hh.pct >= 65]
for name, d, lab_med in [('Fig3d', hh, 'whole hyposmia group'), ('Fig3e', ndf, 'non-deficit stratum')]:
    med = d.OIS.median(); tabs, summ = [], []
    for k_, g in [('high risk (OIS at or below median)', d[d.OIS <= med]), ('low risk (OIS above median)', d[d.OIS > med])]:
        t, two = km_table(g, k_); tabs.append(t); summ.append(two)
    a_, b_ = d[d.OIS <= med], d[d.OIS > med]
    plr = logrank_test(a_.time_years, b_.time_years, a_.converted, b_.converted).p_value
    S = pd.DataFrame(summ); S['OIS_median'] = med; S['log_rank_p'] = plr
    if name == 'Fig3e':
        extra = []
        for col, lab in [('upsit', 'UPSIT'), ('pct', 'lower putamen % expected')]:
            m2 = d[col].median(); x1, x2 = d[d[col] <= m2], d[d[col] > m2]
            k1 = KaplanMeierFitter().fit(x1.time_years, x1.converted); k2 = KaplanMeierFitter().fit(x2.time_years, x2.converted)
            extra.append(dict(split_by=lab, median=m2, two_year_pct_at_or_below=(1 - float(k1.predict(2.))) * 100, two_year_pct_above=(1 - float(k2.predict(2.))) * 100,
                              log_rank_p=logrank_test(x1.time_years, x2.time_years, x1.converted, x2.converted).p_value))
        S = pd.concat([S, pd.DataFrame([{}]), pd.DataFrame(extra)], axis=0)
    add(name, pd.concat([S, pd.DataFrame([{}]), pd.concat(tabs)], axis=0),
        f'Kaplan-Meier curves with Greenwood 95% bands to two years, {lab_med} divided at its OIS median. Summary rows first, then the survival table.')

# ================================================================== Fig 4
m = J['csf_tau']['main']
keys = [('Prodromal|pTau181_over_ABeta42_CSF', 'CSF pTau181/Abeta42'), ('Prodromal|pTau181_CSF', 'CSF pTau181'), ('Prodromal|eMTBR_TAU243_CSF', 'CSF eMTBR-tau243'),
        ('Prodromal|NfL_serum', 'serum NfL'), ('Prodromal|NfL_plasma', 'plasma NfL')]
rows = []
for k, lab in keys:
    if k not in m: continue
    for sc in ['UPSIT', 'OIS', 'OMI']:
        b, p = m[k][f'b_{sc}'], m[k][f'p_{sc}']; se = abs(b) / norm.isf(p / 2) if 0 < p < 1 else 0
        rows.append(dict(analyte=lab, n=m[k]['n'], score=sc, standardised_beta=b, ci_lower=b - 1.96 * se, ci_upper=b + 1.96 * se, p=p))
add('Fig4a', pd.DataFrame(rows), 'Standardised coefficients adjusted for age, sex and years of education, prodromal cohort. CI reconstructed from beta and P as drawn.')
G = pd.DataFrame(J['genotype'])
for name, coh in [('Fig4b', 'Prodromal'), ('Fig4c', 'PD')]:
    g = G[(G.cohort == coh) & G.measure.isin(['OIS', 'OMI', 'UPSIT total'])][['cohort', 'measure', 'n_gba', 'n_lrrk2', 'mean_gba', 'mean_lrrk2', 'beta', 'lo', 'hi', 'p']]
    g = g.rename(columns={'n_gba': 'n_GBA1', 'n_lrrk2': 'n_LRRK2', 'mean_gba': 'mean_GBA1', 'mean_lrrk2': 'mean_LRRK2', 'beta': 'GBA1_minus_LRRK2_adjusted', 'lo': 'ci_lower', 'hi': 'ci_upper'})
    add(name, g, 'Adjusted GBA1 minus LRRK2 difference with 95% CI and P.')

# ================================================================== Supplementary figures
CH = pd.read_csv(f'{R}/supp_channel_followup.csv'); LF = J['layer_flow']
flow = pd.DataFrame([dict(box='Training cohort', n=1398, detail='PD 1,163 + HC 235'),
                     dict(box='Prodromal cohort with imaging and baseline clinical data', n=2453, detail='rule A 1,558 / 534 / 361'),
                     dict(box='Correlation analyses (OIS and UPSIT both present)', n=2441, detail='separate correlation branch from n=2453; excludes 12 without measured UPSIT'),
                     dict(box='After excluding prevalent Parkinson\'s disease at baseline', n=2375, detail='survival branch starts at n=2453; excludes 78 prevalent cases'),
                     dict(box='After excluding no follow-up time', n=1770, detail='excluded 605, 0 converters'),
                     dict(box='After excluding missing analysis variable', n=1759, detail='excluded 11, 152 conversions'),
                     dict(box='DAT deficit stratum (< 65%)', n=265, detail='46 conversions'),
                     dict(box='Non-deficit stratum', n=738, detail='22 conversions')])
add('FigS1', pd.concat([flow, pd.DataFrame([{}]), CH], axis=0), 'Participant flow. Group boxes use the follow-up table (second block).')
add('FigS2a', pd.DataFrame(J['calibration_deciles']), 'Decile calibration in the RBD and variant-carrier group.')
rows = []
d = df[df.COHORT == 4].dropna(subset=['OIS', 'upsit']).copy(); sg = d.subgroup.astype(str)
d['arm'] = np.where(sg == 'Hyposmia', 'Hyposmia', np.where(sg.str.contains('RBD'), 'RBD', 'Genetic'))
for arm, lab in [('ALL', 'All prodromal'), ('Hyposmia', 'Hyposmia cohort'), ('RBD', 'RBD cohort'), ('Genetic', 'Pathogenic-variant cohort')]:
    g = d if arm == 'ALL' else d[d.arm == arm]; s_ = stats.linregress(g.OIS, g.upsit)
    rows.append(dict(recruitment_cohort=lab, n=len(g), slope=s_.slope, slope_ci_lower=s_.slope - 1.96 * s_.stderr, slope_ci_upper=s_.slope + 1.96 * s_.stderr, r=s_.rvalue, p=s_.pvalue))
add('FigS2b', pd.DataFrame(rows), 'Slope of measured UPSIT on OIS by recruitment cohort with 95% CI.')
L = {r['landmark_yr']: r for r in J['landmark']}; TS = J['temporal_split']; CE = J['centre_split']
rows = [dict(analysis='Full hyposmia group', n=L[0]['n'], events=L[0]['events'], delta_C=L[0]['dc'], ci_lower=L[0]['lo'], ci_upper=L[0]['hi'], p=L[0]['p']),
        dict(analysis='Landmark 1 year', n=L[1]['n'], events=L[1]['events'], delta_C=L[1]['dc'], ci_lower=L[1]['lo'], ci_upper=L[1]['hi'], p=L[1]['p']),
        dict(analysis='Landmark 2 years', n=L[2]['n'], events=L[2]['events'], delta_C=L[2]['dc'], ci_lower=L[2]['lo'], ci_upper=L[2]['hi'], p=L[2]['p']),
        dict(analysis='Model retrained on pre-2017 enrolment', n=TS['n'], events=TS['events'], delta_C=TS['dc'], ci_lower=TS['lo'], ci_upper=TS['hi'], p=TS['p'])]
add('FigS3', pd.DataFrame(rows), 'Delta C (OIS minus UPSIT) with 95% CI from 1,000 paired bootstrap resamples and two-sided P values with a reporting floor of 0.001. Landmark risk sets retain participants observed event-free beyond the landmark, with follow-up restarted there.')


def km_three(d, groups, labels, name, note):
    tabs, summ = [], []
    for g_, lab in zip(groups, labels):
        t, two = km_table(d[d.grp == g_], lab); tabs.append(t); summ.append(two)
    S = pd.DataFrame(summ); S['three_group_log_rank_p'] = multivariate_logrank_test(d.time_years, d.grp, d.converted).p_value
    pw = []
    for a, b in [(0, 1), (1, 2), (0, 2)]:
        ga, gb = d[d.grp == groups[a]], d[d.grp == groups[b]]
        pw.append(dict(comparison=f'{labels[a]} vs {labels[b]}', log_rank_p=logrank_test(ga.time_years, gb.time_years, ga.converted, gb.converted).p_value))
    add(name, pd.concat([S, pd.DataFrame([{}]), pd.DataFrame(pw), pd.DataFrame([{}]), pd.concat(tabs)], axis=0), note)


for name, d, ttl in [('FigS4a', hh, 'whole hyposmia group'), ('FigS4b', ndf, 'non-deficit stratum')]:
    d = d.copy(); d['grp'] = pd.qcut(d.OIS, 3, labels=['T1', 'T2', 'T3']).astype(str)
    km_three(d, ['T1', 'T2', 'T3'], ['lowest', 'middle', 'highest'], name, f'OIS tertiles, {ttl}. Summary, pairwise log-rank, then survival tables.')
S44 = pd.read_csv(f'{R}/section44_three_groups.csv'); r44 = S44[(S44.stratum == 'hyposmia') & (S44.reading == 'OIS')].iloc[0]
d = hh.copy(); d['grp'] = np.where(d.OIS <= r44.cut1, 'low', np.where(d.OIS <= r44.cut2, 'mid', 'high'))
km_three(d, ['low', 'mid', 'high'], ['high risk', 'intermediate', 'low risk'], 'FigS5',
         f'OIS bands at {r44.cut1:.4f} and {r44.cut2:.4f} (permutation P {r44.p_permutation}, {int(r44.n_perm)} permutations, bootstrap 95% for the cuts {r44.cut1_boot_lo:.2f} to {r44.cut1_boot_hi:.2f} and {r44.cut2_boot_lo:.2f} to {r44.cut2_boot_hi:.2f}).')
# ================================================================== write
readme = pd.DataFrame([dict(sheet=k, description=v[1]) for k, v in sheets.items()])
readme = pd.concat([pd.DataFrame([dict(sheet='About', description='Source data for the figures of "Decomposition of the olfactory score by dopamine transporter imaging improves Parkinson\'s disease risk stratification in hyposmic individuals". One sheet per panel. Participant-level values are not included (PPMI Data Use Agreement); the underlying data are available on application at www.ppmi-info.org.')]), readme], axis=0)
with pd.ExcelWriter(OUTFILE, engine='openpyxl') as w:
    readme.to_excel(w, sheet_name='README', index=False)
    for k, (fr, note) in sheets.items():
        fr = fr.replace('UPSIT total', 'UPSIT', regex=True)     # ledger label -> manuscript wording (display only)
        fr.to_excel(w, sheet_name=k, index=False)
    for ws in w.book.worksheets:
        for col in ws.columns:
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)
print(f'-> {OUTFILE}: {len(sheets)} panel sheets')
# consistency check against the median-split table used by the text
M = pd.read_csv(f'{R}/supp_median_split.csv')
chk = sheets['Fig3d'][0].iloc[:2]; ref = M[(M.stratum == 'hyposmia') & (M.reading == 'OIS')].iloc[0]
assert abs(chk.two_year_conversion_pct.iloc[0] - ref.km2_low) < 1e-9 and abs(chk.two_year_conversion_pct.iloc[1] - ref.km2_high) < 1e-9, 'Fig3d does not match supp_median_split.csv'
print('   Fig3d two-year rates match supp_median_split.csv')
