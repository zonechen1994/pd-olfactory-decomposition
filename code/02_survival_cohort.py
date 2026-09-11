"""02_survival_cohort.py  Prodromal survival cohort, imaging state, robustness, decomposition

Builds the prodromal survival frame (prevalent cases excluded, 1,759 participants and 152
conversions), the temporal and centre splits, the age/sex-expected putamen metric with its
validation against the official PPMI staging field, test-retest ICC, landmark and time-dependent
analyses in the hyposmia group, the head-to-head against Cox models trained on conversion labels,
calibration in the prodromal cohort, the PARS high-risk subgroup, and the variance decomposition.
Writes the first version of the ledger results/repro/section16_stats.json and
results/intermediate/02_cohort.pkl.
Notebook provenance: sections 6, 14, 12, 16.1, 16.2, 16.3 (validation block only), 16.4, 16.5 and 17. Code is the notebook cell text, unchanged except for
the file paths, which are resolved through common.py.
"""
from common import *
st = load_state('01_model'); globals().update(st)

# ============================================================================
# Section 6. Prodromal conversion survival frame
# ============================================================================
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

_e4 = excel[excel['COHORT'] == 4].sort_values(['PATNO', 'visit_date'])
first_visit = _e4.groupby('PATNO')['visit_date'].min()
last_visit  = _e4.groupby('PATNO')['visit_date'].max()
pd_diag     = _e4[_e4['PRIMDIAG'] == 1].groupby('PATNO')['visit_date'].min()

# PREVALENT CASES. Participants whose baseline visit already carries PRIMDIAG == 1 had manifest
# Parkinson's disease at entry, not an incident conversion. In the hyposmia arm their baseline
# MDS-UPDRS III averages 19.4 against 4.4 in everyone else, so this is manifest motor disease and
# not a coding artefact (PRIMDIAG varies visit to visit, so it is not back-filled). Counting them
# as events would both inflate the event count and, because they carry the most abnormal baseline
# imaging while being assigned the shortest survival time, inflate every C-index by construction.
# They are excluded from the primary survival cohort and retained only for a sensitivity analysis.
prevalent_at_baseline = set(_e4.groupby('PATNO').first().query('PRIMDIAG == 1').index)
print(f'Prevalent PD at the baseline visit, excluded from the survival cohort: '
      f'{len(prevalent_at_baseline)} participants')

prod_with_prevalent = df[df['COHORT'] == 4].copy()   # kept for the sensitivity analysis
prod = df[df['COHORT'] == 4].copy()   # df already contains PUTAMEN_REF_CWM (one of the 33 SBR cols)
prod = prod[~prod['PATNO'].isin(prevalent_at_baseline)].copy()
prod['converted'] = prod['PATNO'].isin(pd_diag.index).astype(int)
prod['first_visit'] = pd.to_datetime(prod['PATNO'].map(first_visit))
prod['last_visit']  = pd.to_datetime(prod['PATNO'].map(last_visit))
prod['pd_date']     = pd.to_datetime(prod['PATNO'].map(pd_diag))
_tconv = (prod['pd_date'] - prod['first_visit']).dt.days / 365.25
_tcens = (prod['last_visit'] - prod['first_visit']).dt.days / 365.25
# 转化者即使诊断在首访同期(t<=0)也保留(=快速转化,下限 0.25 年),保住全部转化事件;删失者需有随访(t>0)
prod['time_years'] = np.where(prod['converted'] == 1, np.maximum(_tconv, 0.25), _tcens)
prod = prod[(prod['converted'] == 1) | (prod['time_years'] > 0)].copy()
print(f'Prodromal with survival data: {len(prod)} | converted: {prod["converted"].sum()} '
      f'({prod["converted"].mean()*100:.1f}%)')

# Sensitivity: the same cohort with prevalent baseline PD retained, i.e. the earlier definition.
_pw = prod_with_prevalent
_pw['converted'] = _pw['PATNO'].isin(pd_diag.index).astype(int)
_pw['first_visit'] = pd.to_datetime(_pw['PATNO'].map(first_visit))
_pw['last_visit']  = pd.to_datetime(_pw['PATNO'].map(last_visit))
_pw['pd_date']     = pd.to_datetime(_pw['PATNO'].map(pd_diag))
_pw['time_years'] = np.where(_pw['converted'] == 1,
                             np.maximum((_pw['pd_date'] - _pw['first_visit']).dt.days / 365.25, 0.25),
                             (_pw['last_visit'] - _pw['first_visit']).dt.days / 365.25)
prod_with_prevalent = _pw[(_pw['converted'] == 1) | (_pw['time_years'] > 0)].copy()
print(f'  sensitivity cohort retaining prevalent PD: {len(prod_with_prevalent)} | '
      f'converted: {prod_with_prevalent["converted"].sum()}')

# ============================================================================
# Section 6 (continued). Concordance by score and by enrolment subgroup
# ============================================================================
# Sign convention: lifelines concordance_index treats the score as survival-time-like
# (higher score -> longer survival). OIS/UPSIT/PUTAMEN are all protective (higher = healthier
# = later conversion), so we pass the score with a POSITIVE sign.
def c_index(d, col):
    v = d[[col, 'time_years', 'converted']].dropna()
    return concordance_index(v['time_years'], v[col], v['converted']), len(v)

print('--- Overall prodromal C-index ---')
for col, lab in [('OIS','DaTscan-OIS'), ('PUTAMEN_REF_CWM','PUTAMEN SBR'), ('upsit','UPSIT (observed)')]:
    c, n = c_index(prod, col)
    print(f'  {lab:20s}: C={c:.3f}  (N={n})')

print('\n--- Subgroup C-index (OIS) ---')
for sg in ['Hyposmia','RBD','LRRK2','GBA']:
    sub = prod[prod['subgroup'] == sg]
    v = sub[['OIS','time_years','converted']].dropna()
    if len(v) >= 20 and v['converted'].sum() >= 3:
        c = concordance_index(v['time_years'], v['OIS'], v['converted'])
        print(f'  {sg:10s}: C={c:.3f}  N={len(v)}  events={v["converted"].sum()}')

# ============================================================================
# Section 6 (continued). Paired bootstrap of C(OIS) minus C(UPSIT) by group
# ============================================================================
# Interaction effect: UPSIT vs OIS C-index in Hyposmia vs Non-Hyposmia + bootstrap dC test
prod['grp'] = np.where(prod['subgroup'] == 'Hyposmia', 'Hyposmia', 'Non-Hyposmia')

def boot_dc(sub, n_boot=1000, seed=42):
    # C(OIS) - C(UPSIT) with bootstrap 95% CI and two-sided p-value
    v = sub[['OIS','upsit','time_years','converted']].dropna()
    v = v[v['converted'].notna()]
    obs = (concordance_index(v['time_years'], v['OIS'], v['converted'])
           - concordance_index(v['time_years'], v['upsit'], v['converted']))
    rng = np.random.RandomState(seed); diffs = []
    idx = np.arange(len(v))
    for _ in range(n_boot):
        b = rng.choice(idx, len(idx), replace=True)
        vb = v.iloc[b]
        if vb['converted'].sum() < 3 or vb['converted'].nunique() < 2:
            continue
        try:
            diffs.append(concordance_index(vb['time_years'], vb['OIS'], vb['converted'])
                         - concordance_index(vb['time_years'], vb['upsit'], vb['converted']))
        except Exception:
            pass
    diffs = np.array(diffs)
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return obs, np.percentile(diffs, 2.5), np.percentile(diffs, 97.5), p, len(v)

print(f"{'Group':14s}{'C(UPSIT)':>10s}{'C(OIS)':>9s}{'dC':>8s}{'95% CI':>20s}{'p':>9s}")
res_int = {}
for g in ['Non-Hyposmia','Hyposmia']:
    sub = prod[prod['grp'] == g]
    cu, _ = c_index(sub, 'upsit'); co, _ = c_index(sub, 'OIS')
    dc, lo, hi, p, n = boot_dc(sub)
    res_int[g] = dc
    print(f'{g:14s}{cu:>10.3f}{co:>9.3f}{dc:>+8.3f}   [{lo:+.3f}, {hi:+.3f}]{p:>9.4f}')
print(f"\nInteraction ratio dC(Hyposmia)/dC(Non-Hyposmia) = {res_int['Hyposmia']/res_int['Non-Hyposmia']:.2f}x")

# ============================================================================
# Section 6 (continued). Kaplan-Meier by OIS tertile and time-dependent AUC
# ============================================================================
# KM by OIS tertiles + time-dependent AUC
from sklearn.metrics import roc_auc_score
fig, (axk, axa) = plt.subplots(1, 2, figsize=(12, 4.5))

prod['OIS_tertile'] = pd.qcut(prod['OIS'], 3, labels=['Low OIS','Mid OIS','High OIS'])
for lab, col in zip(['Low OIS','Mid OIS','High OIS'], ['#d62728','#ff7f0e','#2ca02c']):
    s = prod[prod['OIS_tertile'] == lab]
    kmf = KaplanMeierFitter().fit(s['time_years'], s['converted'],
                                  label=f'{lab} (n={len(s)}, {s["converted"].sum()} conv)')
    kmf.plot_survival_function(ax=axk, color=col)
axk.set_xlabel('Years'); axk.set_ylabel('PD-free survival'); axk.set_title('KM by OIS tertile'); axk.set_ylim(0.4, 1.02)

def td_auc(T, E, score, tp):
    cases = (E == 1) & (T <= tp); controls = T > tp
    if cases.sum() < 5 or controls.sum() < 5: return np.nan
    m = cases | controls
    try: return roc_auc_score(cases[m].astype(int), score[m])
    except Exception: return np.nan

vs = prod.dropna(subset=['OIS','upsit','time_years','converted'])
T, E = vs['time_years'].values, vs['converted'].values
tps = list(range(1, 11))
for lab, sc, col in [('OIS', -vs['OIS'].values, '#16a085'), ('UPSIT', -vs['upsit'].values, '#e74c3c')]:
    axa.plot(tps, [td_auc(T, E, sc, t) for t in tps], 'o-', label=lab, color=col)
axa.axhline(0.5, ls=':', color='grey'); axa.set_ylim(0.45, 1.0)
axa.set_xlabel('Time horizon (years)'); axa.set_ylabel('Time-dependent AUC'); axa.set_title('Time-dependent AUC'); axa.legend()
plt.tight_layout(); plt.savefig(FIG / 'stage3_km_tdauc.png', bbox_inches='tight'); plt.savefig(FIG / 'stage3_km_tdauc.pdf', bbox_inches='tight'); plt.show()

# ============================================================================
# Section 14. Temporal split, expected-putamen metric and its validation, PARS 2x2, centre split
# ============================================================================
# Section 14 (R4 robustness): internal-external validation
# Panel a: temporal split, train OIS on pre-2017 PD+HC, validate on 2017+ subjects
# Panel b: PARS 2017 hyposmia × OIS-tertile 2x2 conversion rate reproduction
from lifelines.utils import concordance_index

# --- A. Temporal split ---
fv_map = excel.groupby('PATNO')['visit_date'].min()
df_t = df.copy()
df_t['first_visit'] = df_t['PATNO'].map(fv_map)
df_t['year_enroll'] = pd.to_datetime(df_t['first_visit']).dt.year

CUTOFF = 2017
train_mask = (df_t['COHORT'].isin([1,2])) & (df_t['year_enroll'] < CUTOFF) & df_t['upsit'].notna()
test_mask  = (df_t['COHORT'].isin([1,2])) & (df_t['year_enroll'] >= CUTOFF) & df_t['upsit'].notna()

X_tr = df_t.loc[train_mask, FEATURES].values
y_tr = df_t.loc[train_mask, 'upsit'].values
X_te = df_t.loc[test_mask,  FEATURES].values
y_te = df_t.loc[test_mask,  'upsit'].values

sc_t = StandardScaler().fit(X_tr)
m_t  = Ridge(alpha=1.0).fit(sc_t.transform(X_tr), y_tr)
yp_te = m_t.predict(sc_t.transform(X_te))
r_temp = stats.pearsonr(y_te, yp_te)[0]
print(f'A. Temporal split (train pre-{CUTOFF} PD+HC, test {CUTOFF}+ PD+HC)')
print(f'   N train = {train_mask.sum()}, N test = {test_mask.sum()}')
print(f'   OIS held-out r on 2017+ PD+HC = {r_temp:.3f}')

# Apply model to ALL 2017+ subjects including prodromal, then run R2 conversion analysis
df_t['OIS_temp'] = m_t.predict(sc_t.transform(df_t[FEATURES].values))
prod_late = df_t[(df_t['COHORT']==4) & (df_t['year_enroll']>=CUTOFF)].copy()
prod_late = prod_late[~prod_late['PATNO'].isin(prevalent_at_baseline)].copy()  # prevalent baseline PD
syn_codes = [1]   # PD-only, matching §6 main R2 analysis (PRIMDIAG=1)
syn_diag = excel[(excel['COHORT']==4) & excel['PRIMDIAG'].isin(syn_codes)].groupby('PATNO')['visit_date'].min()
last_visit = excel.groupby('PATNO')['visit_date'].max()
sg_map = excel[['PATNO','subgroup']].drop_duplicates('PATNO').set_index('PATNO')['subgroup']
prod_late['subgroup'] = prod_late['PATNO'].map(sg_map)
prod_late['converted'] = prod_late['PATNO'].isin(syn_diag.index).astype(int)
prod_late['syn_date']   = pd.to_datetime(prod_late['PATNO'].map(syn_diag))
prod_late['last_visit'] = pd.to_datetime(prod_late['PATNO'].map(last_visit))
prod_late['first_visit'] = pd.to_datetime(prod_late['first_visit'])
_tc = (prod_late['syn_date']  - prod_late['first_visit']).dt.days/365.25
_te2= (prod_late['last_visit']- prod_late['first_visit']).dt.days/365.25
prod_late['time_years'] = np.where(prod_late['converted']==1, np.maximum(_tc, 0.25), _te2)
prod_late = prod_late[((prod_late['converted']==1) | (prod_late['time_years']>0)) & prod_late['OIS_temp'].notna()].copy()
prod_late['grp'] = np.where(prod_late['subgroup']=='Hyposmia','Hyposmia','Non-Hyposmia')
print(f'   2017+ prodromal evaluation cohort: N = {len(prod_late)}, events = {int(prod_late["converted"].sum())}')

c_results = {}
for grp_name, sub in [('Overall', prod_late), ('Non-Hyposmia', prod_late[prod_late['grp']=='Non-Hyposmia']), ('Hyposmia', prod_late[prod_late['grp']=='Hyposmia'])]:
    vu = sub[['upsit','time_years','converted']].dropna()
    vo = sub[['OIS_temp','time_years','converted']].dropna()
    c_u = concordance_index(vu['time_years'], vu['upsit'], vu['converted']) if len(vu)>=20 else np.nan
    c_o = concordance_index(vo['time_years'], vo['OIS_temp'], vo['converted']) if len(vo)>=20 else np.nan
    c_results[grp_name] = (c_u, c_o, len(sub), int(sub['converted'].sum()))
    print(f'   {grp_name:14s}: N={len(sub):4d}, ev={int(sub["converted"].sum()):3d} | C(UPSIT)={c_u:.3f}, C(OIS_temporal)={c_o:.3f}')

# 配对自助法给出 dC 的 95% CI 与 P,口径与 §16 的 landmark 行一致(1000 次,seed 42),
# 否则这一行在 Supplementary Fig 3 的森林图上是个没有区间的裸点估计。
_hyp_late = prod_late[prod_late['grp'] == 'Hyposmia'].dropna(subset=['upsit', 'OIS_temp', 'time_years', 'converted'])
_rngT = np.random.default_rng(42)
_dT = []
for _ in range(1000):
    _b = _hyp_late.iloc[_rngT.choice(len(_hyp_late), len(_hyp_late), replace=True)]
    if _b['converted'].sum() < 3:
        continue
    try:
        _dT.append(concordance_index(_b['time_years'], _b['OIS_temp'], _b['converted'])
                   - concordance_index(_b['time_years'], _b['upsit'], _b['converted']))
    except Exception:
        pass
_dT = np.array(_dT)
_cuT = concordance_index(_hyp_late['time_years'], _hyp_late['upsit'], _hyp_late['converted'])
_coT = concordance_index(_hyp_late['time_years'], _hyp_late['OIS_temp'], _hyp_late['converted'])
_loT, _hiT = np.percentile(_dT, [2.5, 97.5])
_pT = max(2 * min((_dT <= 0).mean(), (_dT >= 0).mean()), 1 / len(_dT))
# 与主分析(1003/68)的人群重叠,写正文时必须交代
_overlapT = len(set(_hyp_late['PATNO']) & set(prod[prod['subgroup'] == 'Hyposmia']['PATNO']))
TEMPORAL_SPLIT = dict(cutoff=int(CUTOFF), n_train=int(train_mask.sum()), n_test_pdhc=int(test_mask.sum()),
                      r_heldout=float(r_temp), n=int(len(_hyp_late)),
                      events=int(_hyp_late['converted'].sum()), c_upsit=float(_cuT), c_ois=float(_coT),
                      dc=float(_coT - _cuT), lo=float(_loT), hi=float(_hiT), p=float(_pT),
                      overlap_with_main=int(_overlapT))
print(f'   Hyposmia dC = {_coT - _cuT:+.3f} [{_loT:+.3f}, {_hiT:+.3f}], P = {_pT:.4f}'
      f'  (training set {int(train_mask.sum())} of 1,398)')
print(f'   test set overlaps the main hyposmia group in {_overlapT} of its {len(_hyp_late)} participants')

# --- B. PARS-style 2x2 reproduction (uses PARS-original putamen DAT-deficit definition) ---
# Canonical PPMI operationalisation of "percent of expected DAT binding".
# Chahine et al., Ann Clin Transl Neurol 2020;8(1):201-212 (PMID 33321002) define it as the
# LOWEST putamen SBR (the lower of left and right) divided by the value expected for that
# participant's AGE AND SEX, the expected value coming from a linear regression fitted in the
# PPMI healthy controls. PPMI publishes the resulting flag as the curated Stage_D field, but
# only for the PD and SWEDD cohorts; Stage_D is empty for every prodromal and healthy-control
# participant, so it has to be recomputed here. The implementation below is therefore first
# validated against the official field in the cohorts where PPMI does release it.
print(f'\nB. PARS-style 2x2 reproduction using PARS-original putamen DAT-deficit definition')
age_col = 'age' if 'age' in df.columns else 'age_at_visit'

def _put_low(frame):
    # Lower of the left/right putamen SBR (the PPMI 'lowest putamen')
    return frame[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1)

_hcn = df[df['COHORT'] == 2].copy()
_hcn['put_low'] = _put_low(_hcn)
_hcn = _hcn.dropna(subset=['put_low', age_col, 'sex_num'])
_Xn = np.column_stack([np.ones(len(_hcn)), _hcn[age_col].astype(float).values,
                       _hcn['sex_num'].astype(float).values])
NORM_BETA, *_ = np.linalg.lstsq(_Xn, _hcn['put_low'].astype(float).values, rcond=None)

def pct_expected_putamen(frame, age_override=None):
    # Percent of age/sex-expected lowest putamen SBR (PPMI canonical metric)
    a = (frame[age_col].astype(float) if age_override is None
         else pd.Series(np.asarray(age_override, dtype=float), index=frame.index))
    exp = NORM_BETA[0] + NORM_BETA[1] * a + NORM_BETA[2] * frame['sex_num'].astype(float)
    return _put_low(frame) / exp

# Validation against the official curated Stage_D field (released for PD and SWEDD only).
# Stage_D is defined in the PPMI data dictionary as age/sex-expected DaT ratio < 0.75.
_sd = excel[excel['EVENT_ID'] == 'BL'][['PATNO', 'Stage_D']].drop_duplicates('PATNO')
_val = df[df['COHORT'].isin([1, 3])].merge(_sd, on='PATNO', how='left').copy()
_val['pct_exp'] = pct_expected_putamen(_val)
_val = _val.dropna(subset=['pct_exp', 'Stage_D'])
_disc = int(((_val['pct_exp'] < 0.75).astype(int) != _val['Stage_D']).sum())
_pos_max = float(_val.loc[_val['Stage_D'] == 1, 'pct_exp'].max())
_neg_min = float(_val.loc[_val['Stage_D'] == 0, 'pct_exp'].min())
DAT_VALIDATION = {
    'n_validated': int(len(_val)),
    'agreement_vs_official_Stage_D': float(1 - _disc / len(_val)),
    'n_discordant': _disc,
    'official_deficit_max_pct_exp': _pos_max,
    'official_normal_min_pct_exp': _neg_min,
    'n_hc_normative': int(len(_hcn)),
    'norm_intercept': float(NORM_BETA[0]),
    'norm_beta_age': float(NORM_BETA[1]),
    'norm_beta_sex': float(NORM_BETA[2]),
}
print(f'   %-expected metric: lowest putamen / (age, sex)-expected value, HC normative n={len(_hcn)}')
print(f'   expected = {NORM_BETA[0]:.4f} {NORM_BETA[1]:+.5f}*age {NORM_BETA[2]:+.4f}*sex')
print(f'   validated against official PPMI Stage_D (cut 0.75): '
      f'{DAT_VALIDATION["agreement_vs_official_Stage_D"]*100:.1f}% agreement in {len(_val)} '
      f'PD/SWEDD participants, {_disc} discordant')
print(f'   official deficit max = {_pos_max:.4f}, official normal min = {_neg_min:.4f} '
      f'(discordance confined to a {_pos_max - _neg_min:.3f}-wide band around the cut)')
print(f'   Stage_D is not released for prodromal or HC participants, so the same validated '
      f'formula is applied there')

# PARS threshold: DAT deficit = below 65% of expected (Jennings 2017, JAMA Neurol 74:933-940).
# PARS additionally treats 65-80% as indeterminate; the binary rule below follows the PARS
# primary contrast (deficit vs not deficit) used for the 2x2 reproduction.
PARS_DAT_CUT = 0.65
# Reconstruct full prodromal survival (mirrors §6)
prod_full = df[df['COHORT']==4][['PATNO','OIS','upsit','PUTAMEN_REF_CWM',
                                 'PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM','sex_num',age_col]].copy()
prod_full = prod_full[~prod_full['PATNO'].isin(prevalent_at_baseline)].copy()  # prevalent baseline PD
prod_full['subgroup'] = prod_full['PATNO'].map(sg_map)
prod_full['converted'] = prod_full['PATNO'].isin(syn_diag.index).astype(int)
prod_full['first_visit'] = pd.to_datetime(prod_full['PATNO'].map(excel.groupby('PATNO')['visit_date'].min()))
prod_full['last_visit']  = pd.to_datetime(prod_full['PATNO'].map(last_visit))
prod_full['syn_date']    = pd.to_datetime(prod_full['PATNO'].map(syn_diag))
_tc = (prod_full['syn_date']  - prod_full['first_visit']).dt.days/365.25
_te2= (prod_full['last_visit']- prod_full['first_visit']).dt.days/365.25
prod_full['time_years'] = np.where(prod_full['converted']==1, np.maximum(_tc, 0.25), _te2)
prod_full = prod_full[(prod_full['converted']==1) | (prod_full['time_years']>0)].copy()
prod_full['grp'] = np.where(prod_full['subgroup']=='Hyposmia','Hyposmia','Non-Hyposmia')

# PARS DAT-deficit = lowest putamen SBR < 65% of the age/sex-expected HC value
prod_full['put_pct'] = pct_expected_putamen(prod_full)
prod_full['dat_deficit'] = (prod_full['put_pct'] < PARS_DAT_CUT).astype(int)

# 4-year conversion rate per cell (truncate follow-up to 4 yr to match PARS reporting window)
T_CUT = 4.0
def rate_4yr(sub):
    if len(sub) < 20: return None, len(sub), 0
    sub2 = sub.copy()
    sub2['event_4y'] = ((sub2['converted']==1) & (sub2['time_years']<=T_CUT)).astype(int)
    follow = sub2[(sub2['time_years']>=T_CUT) | (sub2['converted']==1)]
    rate = follow['event_4y'].mean()*100
    return rate, len(follow), int(follow['event_4y'].sum())

print(f'   Age/sex-expected HC lowest putamen SBR: '
      f'SBR = {NORM_BETA[0]:.3f} {NORM_BETA[1]:+.5f}*age {NORM_BETA[2]:+.4f}*sex   (n_HC={len(_hcn)})')
print(f'   DAT-deficit = lowest putamen SBR < {PARS_DAT_CUT:.0%} of the age/sex-expected value')
print(f'   4-yr conversion rate by hyposmia × putamen-defined DAT-deficit:')
cells_2x2 = {}
for h in ['Hyposmia','Non-Hyposmia']:
    for d, dlab in [(1,'DAT-deficit'), (0,'DAT-normal')]:
        sub = prod_full[(prod_full['grp']==h) & (prod_full['dat_deficit']==d)]
        rate, n_follow, n_ev = rate_4yr(sub)
        cells_2x2[(h,dlab)] = (rate, n_follow, n_ev)
        rate_str = f'{rate:.1f}%' if rate is not None else 'NA'
        print(f'     {h:13s} + {dlab:11s}: 4-yr rate = {rate_str:6s}  (N_followup={n_follow}, events={n_ev})')

# PARS published numbers for direct comparison
print(f'\n   PARS 2017 reported (4-yr conversion):')
print(f'     Hyposmia + DAT-deficit: ~67% (Jennings et al. JAMA Neurol 2017)')
print(f'     Hyposmia + DAT-normal:  ~2.8%')

# --- C. Single-shot prodromal site split: discovery vs external validation ---
# OIS trained on FULL PD+HC (matches §3 / §6), prodromal never seen by OIS.
# Among 11 large prodromal sites (n>=60, same threshold as §3 per-site QC),
# randomly assign 5 to "discovery" cohort and 6 to "external validation" cohort (seed=41).
print(f'\nC. Single-shot prodromal site split: discovery vs external validation')
prod_surv = prod_full[['PATNO','grp','converted','time_years','OIS','upsit']].copy()
prod_site = prod_surv.merge(df[['PATNO','SITE']], on='PATNO', how='left')
site_n_raw = df[df['COHORT']==4]['SITE'].value_counts()
large_sites = sorted(site_n_raw[site_n_raw >= 60].index.tolist())
print(f'   Candidate large prodromal sites (n>=60 raw): {len(large_sites)}')

# Reproducible single-shot random split
rng = np.random.RandomState(14)
shuffled = rng.permutation(large_sites)
discovery_sites = sorted(shuffled[:5].tolist())
external_sites  = sorted(shuffled[5:].tolist())
print(f'   Discovery sites (5): {discovery_sites}')
print(f'   External sites  (6): {external_sites}')

def cohort_R2(prod_subset, label):
    hyp = prod_subset[prod_subset['grp']=='Hyposmia']
    nh  = prod_subset[prod_subset['grp']=='Non-Hyposmia']
    v_h = hyp[['OIS','upsit','time_years','converted']].dropna()
    v_n = nh[['OIS','upsit','time_years','converted']].dropna()
    c_u_h = concordance_index(v_h['time_years'], v_h['upsit'], v_h['converted'])
    c_o_h = concordance_index(v_h['time_years'], v_h['OIS'],   v_h['converted'])
    c_u_n = concordance_index(v_n['time_years'], v_n['upsit'], v_n['converted'])
    c_o_n = concordance_index(v_n['time_years'], v_n['OIS'],   v_n['converted'])
    # bootstrap dC CI for hyposmia
    rng_b = np.random.RandomState(42); diffs = []
    idx = np.arange(len(v_h))
    for _ in range(1000):
        b = rng_b.choice(idx, len(idx), replace=True)
        vb = v_h.iloc[b]
        if vb['converted'].sum() < 3 or vb['converted'].nunique() < 2: continue
        try:
            diffs.append(concordance_index(vb['time_years'], vb['OIS'],   vb['converted'])
                       - concordance_index(vb['time_years'], vb['upsit'], vb['converted']))
        except Exception: pass
    diffs = np.array(diffs)
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    print(f'\n   {label}: {len(prod_subset)} prodromals, {int(prod_subset["converted"].sum())} events')
    print(f'     Hyposmia    : N={len(v_h):3d}, ev={int(v_h["converted"].sum()):3d} | C(UPSIT)={c_u_h:.3f} → C(OIS)={c_o_h:.3f}, ΔC=+{c_o_h-c_u_h:.3f}'
          f' [95% CI {np.percentile(diffs,2.5):+.3f}, {np.percentile(diffs,97.5):+.3f}], p={p:.4f}')
    print(f'     Non-Hyp     : N={len(v_n):3d}, ev={int(v_n["converted"].sum()):3d} | C(UPSIT)={c_u_n:.3f} → C(OIS)={c_o_n:.3f}, ΔC=+{c_o_n-c_u_n:.3f}')
    return dict(label=label, n=len(prod_subset), ev=int(prod_subset['converted'].sum()),
                hyp_n=len(v_h), hyp_ev=int(v_h['converted'].sum()),
                c_ups_hyp=c_u_h, c_ois_hyp=c_o_h, dC_hyp=c_o_h-c_u_h,
                ci_lo=np.percentile(diffs,2.5), ci_hi=np.percentile(diffs,97.5), p_boot=p,
                c_ups_nh=c_u_n, c_ois_nh=c_o_n, dC_nh=c_o_n-c_u_n)

disc = cohort_R2(prod_site[prod_site['SITE'].isin(discovery_sites)], 'Discovery cohort (5 sites)')
extv = cohort_R2(prod_site[prod_site['SITE'].isin(external_sites)],  'External validation (6 sites)')

# Per-site breakdown within each cohort (hyposmics)
def per_site_breakdown(cohort_sites, label):
    print(f'\n   {label}, per-site hyposmia R2 interaction:')
    rows = []
    for s in cohort_sites:
        sub = prod_site[prod_site['SITE']==s]
        hyp = sub[sub['grp']=='Hyposmia']
        v = hyp[['OIS','upsit','time_years','converted']].dropna()
        n, ev = len(v), int(v['converted'].sum())
        if n < 5 or ev < 1:
            print(f'     site={int(s):4d}: n_hyp={n:3d}, ev={ev:2d} | (skipped, n<5 or no events)')
            rows.append({'site': int(s), 'n_hyp': n, 'ev_hyp': ev,
                         'c_ups': np.nan, 'c_ois': np.nan, 'dC': np.nan})
            continue
        c_u = concordance_index(v['time_years'], v['upsit'], v['converted'])
        c_o = concordance_index(v['time_years'], v['OIS'],   v['converted'])
        print(f'     site={int(s):4d}: n_hyp={n:3d}, ev={ev:2d} | '
              f'C(UPSIT)={c_u:.3f}, C(OIS)={c_o:.3f}, ΔC={c_o-c_u:+.3f}')
        rows.append({'site': int(s), 'n_hyp': n, 'ev_hyp': ev,
                     'c_ups': c_u, 'c_ois': c_o, 'dC': c_o-c_u})
    return pd.DataFrame(rows)

disc_per = per_site_breakdown(discovery_sites, 'Discovery (5 sites)')
extv_per = per_site_breakdown(external_sites,  'External validation (6 sites)')

# Save tables for reference
disc_per.to_csv(OUT / 'fig4h_discovery_per_site.csv', index=False)
extv_per.to_csv(OUT / 'fig4h_external_per_site.csv', index=False)

# Summary across External sites
ev_ok = extv_per.dropna(subset=['dC'])
print(f'\n   External site-level summary (n_evaluable={len(ev_ok)} of 6):')
print(f'     mean ΔC across External sites = {ev_ok["dC"].mean():+.3f} (range [{ev_ok["dC"].min():+.3f}, {ev_ok["dC"].max():+.3f}])')
print(f'     fraction External sites with ΔC > 0: {(ev_ok["dC"]>0).mean()*100:.0f}% ({int((ev_ok["dC"]>0).sum())}/{len(ev_ok)})')

# === Fig 4h ===
fig, ax = plt.subplots(1, 3, figsize=(18.5, 5.0))

# Panel a: temporal split C-index bars
groups_a = ['Overall', 'Non-Hyposmia', 'Hyposmia']
xs = np.arange(len(groups_a))
c_u = [c_results[g][0] for g in groups_a]
c_o = [c_results[g][1] for g in groups_a]
w = 0.36
ax[0].bar(xs - w/2, c_u, w, label='UPSIT (raw)', color='#C44E52')
ax[0].bar(xs + w/2, c_o, w, label='OIS (trained on pre-2017)', color='#16A085')
for i, (cu, co) in enumerate(zip(c_u, c_o)):
    if not np.isnan(cu): ax[0].text(i-w/2, cu+0.01, f'{cu:.3f}', ha='center', fontsize=9)
    if not np.isnan(co): ax[0].text(i+w/2, co+0.01, f'{co:.3f}', ha='center', fontsize=9, fontweight='bold')
ax[0].set_xticks(xs); ax[0].set_xticklabels(groups_a)
ax[0].axhline(0.5, ls=':', color='grey')
ax[0].set_ylim(0.45, 0.85)
ax[0].set_ylabel('C-index, prodromal conversion (PD+DLB+MSA)')
ax[0].set_title(f'a  Temporal split: OIS trained pre-{CUTOFF}, tested on {CUTOFF}+ prodromals\n'
                f'    (training N={train_mask.sum()}, test N={len(prod_late)}, events={int(prod_late["converted"].sum())})',
                loc='left', fontsize=10.5, fontweight='bold')
ax[0].legend(loc='upper left', fontsize=9, frameon=False)

# Panel b: OIS continues to stratify WITHIN PARS-defined hyposmia + DAT-deficit high-risk group
# (PPMI prodromal enrichment makes absolute PARS rate reproduction impossible, so we instead
#  show that OIS provides continuous stratification on top of the PARS binary scheme)
from lifelines import KaplanMeierFitter
hyp_dd = prod_full[(prod_full['grp']=='Hyposmia') & (prod_full['dat_deficit']==1)].copy()
hyp_dd['OIS_t'] = pd.qcut(hyp_dd['OIS'], 3, labels=['Q1 low OIS','Q2','Q3 high OIS'], duplicates='drop')
kmf = KaplanMeierFitter()
tert_cols = {'Q1 low OIS':'#7B241C', 'Q2':'#C0392B', 'Q3 high OIS':'#F5B7B1'}
for t, sub in hyp_dd.groupby('OIS_t', observed=True):
    kmf.fit(sub['time_years'], event_observed=sub['converted'],
            label=f'{t} (n={len(sub)}, ev={int(sub["converted"].sum())})')
    kmf.plot_survival_function(ax=ax[1], color=tert_cols[t], ci_show=False, lw=2.2)
# C-index within hyposmia + DAT-deficit
v_dd = hyp_dd[['OIS','upsit','time_years','converted']].dropna()
c_ois_dd = concordance_index(v_dd['time_years'], v_dd['OIS'], v_dd['converted'])
c_ups_dd = concordance_index(v_dd['time_years'], v_dd['upsit'], v_dd['converted'])
ax[1].set_xlim(0, 8); ax[1].set_ylim(0.0, 1.02)
ax[1].set_xlabel('Years from baseline'); ax[1].set_ylabel('Conversion-free survival')
ax[1].legend(loc='lower left', fontsize=9, frameon=False)
ax[1].set_title(f'b  OIS tertile within hyposmia + DAT-deficit (PARS high-risk)\n'
                f'    C(UPSIT)={c_ups_dd:.3f} → C(OIS)={c_ois_dd:.3f}, N={len(hyp_dd)}, ev={int(hyp_dd["converted"].sum())}',
                loc='left', fontsize=10.5, fontweight='bold')

# Panel c: forest plot (per-site ΔC for >10 events sites + pooled cohort estimates)
# Only show per-site rows where events > 10 (EPV threshold for Cox C-index stability).
disc_per_ok = disc_per[disc_per['ev_hyp'] > 10].copy()
extv_per_ok = extv_per[extv_per['ev_hyp'] > 10].copy()
print(f'\n   Per-site rows for forest plot (events > 10 only):')
print(f'     Discovery: {len(disc_per_ok)} sites')
print(f'     External : {len(extv_per_ok)} sites')

forest_rows = []
for _, r in disc_per_ok.iterrows():
    forest_rows.append({'label': f'site {int(r["site"])}  (n={int(r["n_hyp"])}, ev={int(r["ev_hyp"])})',
                        'cohort': 'Discovery', 'c_ups': r['c_ups'], 'c_ois': r['c_ois'],
                        'dC': r['dC'], 'is_pool': False, 'ci_lo': np.nan, 'ci_hi': np.nan, 'p': np.nan})
forest_rows.append({'label': f'Discovery pooled  (5 sites, ev={disc["hyp_ev"]})',
                    'cohort': 'Discovery', 'c_ups': disc['c_ups_hyp'], 'c_ois': disc['c_ois_hyp'],
                    'dC': disc['dC_hyp'], 'is_pool': True, 'ci_lo': disc['ci_lo'], 'ci_hi': disc['ci_hi'],
                    'p': disc['p_boot']})
for _, r in extv_per_ok.iterrows():
    forest_rows.append({'label': f'site {int(r["site"])}  (n={int(r["n_hyp"])}, ev={int(r["ev_hyp"])})',
                        'cohort': 'External', 'c_ups': r['c_ups'], 'c_ois': r['c_ois'],
                        'dC': r['dC'], 'is_pool': False, 'ci_lo': np.nan, 'ci_hi': np.nan, 'p': np.nan})
forest_rows.append({'label': f'External pooled  (6 sites, ev={extv["hyp_ev"]})',
                    'cohort': 'External', 'c_ups': extv['c_ups_hyp'], 'c_ois': extv['c_ois_hyp'],
                    'dC': extv['dC_hyp'], 'is_pool': True, 'ci_lo': extv['ci_lo'], 'ci_hi': extv['ci_hi'],
                    'p': extv['p_boot']})
forest_df = pd.DataFrame(forest_rows)

ys_f = np.arange(len(forest_df))[::-1]
for yi, (_, row) in zip(ys_f, forest_df.iterrows()):
    lw_, sz_, fw_ = (2.5, 130, 'bold') if row['is_pool'] else (1.2, 60, 'normal')
    ax[2].plot([row['c_ups'], row['c_ois']], [yi, yi], '-', color='grey', lw=lw_, zorder=1)
    ax[2].scatter(row['c_ups'], yi, color='#C44E52', s=sz_, zorder=2,
                  label='UPSIT' if yi == ys_f[0] else None)
    ax[2].scatter(row['c_ois'], yi, color='#16A085', s=sz_, zorder=2,
                  label='OIS' if yi == ys_f[0] else None)
    if row['is_pool']:
        ax[2].text(0.985, yi + 0.18,
                   f'ΔC=+{row["dC"]:.3f}  [95% CI {row["ci_lo"]:+.3f}, {row["ci_hi"]:+.3f}]  p={row["p"]:.4f}',
                   ha='right', va='bottom', fontsize=8.8, fontweight='bold')
    else:
        ax[2].text(0.985, yi, f'ΔC={row["dC"]:+.3f}',
                   ha='right', va='center', fontsize=8.5)
# horizontal separator between Discovery block and External block
sep_y = ys_f[(forest_df['cohort'] == 'Discovery').sum() - 1] - 0.5
ax[2].axhline(sep_y, color='grey', lw=0.6, alpha=0.5)
ax[2].axvline(0.5, ls=':', color='grey')
ax[2].set_yticks(ys_f)
ax[2].set_yticklabels(forest_df['label'].values, fontsize=8.5)
ax[2].set_xlim(0.40, 1.00); ax[2].set_ylim(min(ys_f) - 0.6, max(ys_f) + 0.8)
ax[2].set_xlabel('C-index, hyposmic prodromal conversion (PD+DLB+MSA)')
ax[2].legend(loc='lower left', fontsize=9, frameon=False)
ax[2].set_title(f'c  Multi-site external validation: per-site (events ≥10) + pooled cohort\n'
                f'    OIS trained on full PD+HC (n≈1398, no prodromal exposure)',
                loc='left', fontsize=10.5, fontweight='bold')

plt.tight_layout()
plt.savefig(FIG / 'Fig4h_R4_internal_external.png', dpi=300, bbox_inches='tight')
plt.savefig(FIG / 'Fig4h_R4_internal_external.pdf', bbox_inches='tight')
plt.show()
print(f'\nFig 4h saved. Headline numbers:')
print(f'  a) Temporal split: hyposmia C(UPSIT)={c_results["Hyposmia"][0]:.3f} -> C(OIS_temp)={c_results["Hyposmia"][1]:.3f}')
print(f'  b) Within PARS high-risk (hyp+DAT-def, N={len(hyp_dd)}): C(UPSIT)={c_ups_dd:.3f} -> C(OIS)={c_ois_dd:.3f}')
print(f'  c) Discovery (5 sites) hyposmia ΔC = +{disc["dC_hyp"]:.3f} [{disc["ci_lo"]:+.3f}, {disc["ci_hi"]:+.3f}], p={disc["p_boot"]:.4f}')
print(f'     External  (6 sites) hyposmia ΔC = +{extv["dC_hyp"]:.3f} [{extv["ci_lo"]:+.3f}, {extv["ci_hi"]:+.3f}], p={extv["p_boot"]:.4f}')

# ============================================================================
# Centre split summary (values previously copied by hand into the Supplementary generator)
# ============================================================================
CENTRE_SPLIT = {}
for _lab, _r, _sites in [('discovery', disc, discovery_sites), ('external', extv, external_sites)]:
    CENTRE_SPLIT[_lab] = dict(sites=[int(x) for x in _sites], n=int(_r['n']), ev=int(_r['ev']),
                              n_hyp=int(_r['hyp_n']), ev_hyp=int(_r['hyp_ev']),
                              c_upsit=float(_r['c_ups_hyp']), c_ois=float(_r['c_ois_hyp']),
                              dc=float(_r['dC_hyp']), lo=float(_r['ci_lo']), hi=float(_r['ci_hi']),
                              p=float(_r['p_boot']),
                              c_upsit_nh=float(_r['c_ups_nh']), c_ois_nh=float(_r['c_ois_nh']))
print('CENTRE_SPLIT:', json.dumps(CENTRE_SPLIT, indent=1, default=float))

# ============================================================================
# Section 12. Test-retest reliability of the residual
# ============================================================================
# Section 12 (R4): OMI longitudinal test-retest stability (trait-like vs noise)
ois_map = df.set_index('PATNO')['OIS']
# 全文排除 SWEDD(COHORT 3),重测配对必须同口径,否则合并 ICC 里会混进 SWEDD 的配对
_keep12 = set(df.loc[df['COHORT'].isin([1, 2, 4]), 'PATNO'])
lon = excel[['PATNO', 'EVENT_ID', 'upsit']].dropna(subset=['upsit']).copy()
lon = lon[lon['PATNO'].isin(ois_map.index) & lon['PATNO'].isin(_keep12)]
lon['OMI'] = lon['upsit'] - lon['PATNO'].map(ois_map)
nvis = lon.groupby('PATNO').size()
lon = lon[lon['PATNO'].isin(nvis[nvis >= 2].index)]
print(f'有 >=2 次 UPSIT 且具 OIS 的受试者: {lon["PATNO"].nunique()} (记录 {len(lon)})')

def test_retest(long, val):
    w = long.sort_values(['PATNO', 'EVENT_ID']).groupby('PATNO')[val].apply(lambda s: s.iloc[:2].tolist())
    w = w[w.apply(len) == 2]
    a = np.array([x[0] for x in w]); b = np.array([x[1] for x in w])
    M = np.vstack([a, b]).T
    MSB = 2 * np.var(M.mean(1), ddof=1)
    MSW = np.mean([np.var(r, ddof=1) for r in M])
    return stats.pearsonr(a, b)[0], (MSB - MSW) / (MSB + MSW), a, b

r_u, icc_u, _, _ = test_retest(lon, 'upsit')
r_o, icc_o, a_o, b_o = test_retest(lon, 'OMI')

# 分组 ICC 进总账。此前只有一份 2026-09-04 的静态 CSV,不由 notebook 生成,属复现性缺口
_cm12 = df.set_index('PATNO')[['COHORT', 'subgroup']]
def _grp12(pn):
    r = _cm12.loc[pn]
    if r['COHORT'] == 1: return 'PD'
    if r['COHORT'] == 2: return 'HC'
    s = str(r['subgroup'])
    return 'Hyposmia' if s == 'Hyposmia' else ('RBD' if 'RBD' in s else 'Genetic')
lon['grp12'] = [_grp12(x) for x in lon['PATNO']]
_icc12 = []
for g, sub in lon.groupby('grp12'):
    if sub['PATNO'].nunique() < 20:
        continue
    ru, iu, _, _ = test_retest(sub, 'upsit')
    ro, io, aa, _ = test_retest(sub, 'OMI')
    _icc12.append(dict(group=g, n_pairs=int(len(aa)), upsit_r=float(ru), upsit_icc=float(iu),
                       omi_r=float(ro), omi_icc=float(io)))
print('  ICC by group (SWEDD excluded):')
for r in sorted(_icc12, key=lambda x: -x['n_pairs']):
    print(f"    {r['group']:<10}{r['n_pairs']:>6} pairs   UPSIT ICC {r['upsit_icc']:.3f}   OMI ICC {r['omi_icc']:.3f}")
ICC_STATS = dict(n_pairs=int(len(a_o)), upsit_r=float(r_u), upsit_icc=float(icc_u),
                 omi_r=float(r_o), omi_icc=float(icc_o), by_group=_icc12)
pd.DataFrame(_icc12).to_csv(OUT / 'supp_icc_by_group.csv', index=False)
print(f'UPSIT : first-vs-second r={r_u:.3f}, ICC(1,1)={icc_u:.3f}')
print(f'OMI   : first-vs-second r={r_o:.3f}, ICC(1,1)={icc_o:.3f}  (n={len(a_o)} pairs)')
print('=> OMI 重测信度接近 UPSIT 本身,OMI 是稳定的个体特征(trait-like),而非随机噪声。')

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
ax[0].scatter(a_o, b_o, s=8, alpha=0.3, color='#8E44AD')
lim = [min(a_o.min(), b_o.min()), max(a_o.max(), b_o.max())]
ax[0].plot(lim, lim, 'k--', lw=1, alpha=0.6)
ax[0].set_xlabel('OMI(首次随访)'); ax[0].set_ylabel('OMI(第二次随访)')
ax[0].set_title(f'Fig 4d a  OMI 重测一致性 (r={r_o:.2f}, ICC={icc_o:.2f})', loc='left', fontsize=10.5, fontweight='bold')
bars = ax[1].bar(['UPSIT', 'OMI'], [icc_u, icc_o], color=['#C44E52', '#8E44AD'], width=0.6)
ax[1].axhline(0.75, ls=':', color='grey'); ax[1].set_ylim(0, 1); ax[1].set_ylabel('ICC(1,1) 重测信度')
for bb, v in zip(bars, [icc_u, icc_o]):
    ax[1].text(bb.get_x() + bb.get_width() / 2, v + 0.02, f'{v:.2f}', ha='center', fontweight='bold')
ax[1].set_title('Fig 4d b  OMI 重测信度接近 UPSIT 本身', loc='left', fontsize=10.5, fontweight='bold')
plt.tight_layout(); plt.savefig(FIG / 'Fig4d_R4_OMI_stability.png', dpi=150, bbox_inches='tight'); plt.savefig(FIG / 'Fig4d_R4_OMI_stability.pdf', bbox_inches='tight'); plt.show()

# ============================================================================
# Section 16.1. Landmark sensitivity and time-dependent AUC in the hyposmia group; ledger initialised
# ============================================================================
# Section 16.1 (R2.2): landmark sensitivity against reverse causation + time-dependent AUC,
# both inside the hyposmia stratum where the R2 interaction lives.
# These analyses previously existed only as throwaway /tmp scripts; they are the support for
# manuscript R2.2 and are now reproducible from this notebook.
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

sec16 = {}
sec16['temporal_split'] = TEMPORAL_SPLIT      # §14 算出,在此并入总账
sec16['icc'] = ICC_STATS                       # §12 算出,2026-09-10 起进总账(此前是手抄常量 + 静态 CSV)
prod['grp'] = np.where(prod['subgroup'] == 'Hyposmia', 'Hyposmia', 'Non-Hyposmia')
hyp16 = prod[prod['grp'] == 'Hyposmia'].dropna(subset=['OIS', 'upsit', 'time_years', 'converted']).copy()
print(f'§16 hyposmia stratum: N={len(hyp16)}, events={int(hyp16["converted"].sum())}')


def _c16(d, col):
    return concordance_index(d['time_years'], d[col], d['converted'])


def _pfmt(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


def _boot_dc16(d, colA, colB, n_boot=1000, seed=42, one_sided=False):
    # C(colA) - C(colB), bootstrap percentile CI. one_sided=True reports p(dC <= 0).
    v = d[[colA, colB, 'time_years', 'converted']].dropna()
    obs = _c16(v, colA) - _c16(v, colB)
    rng = np.random.RandomState(seed)
    idx = np.arange(len(v)); diffs = []
    for _ in range(n_boot):
        vb = v.iloc[rng.choice(idx, len(idx), replace=True)]
        if vb['converted'].sum() < 3 or vb['converted'].nunique() < 2:
            continue
        try:
            diffs.append(_c16(vb, colA) - _c16(vb, colB))
        except Exception:
            pass
    diffs = np.array(diffs)
    p = (diffs <= 0).mean() if one_sided else 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return obs, np.percentile(diffs, 2.5), np.percentile(diffs, 97.5), p, len(v)


print('\n--- §16.1a Landmark: drop hyposmic converters within L years, recompute dC ---')
print(f"{'Landmark':<10}{'N':>6}{'events':>8}{'C(UPSIT)':>10}{'C(OIS)':>9}{'dC':>9}{'95% CI':>21}{'p':>9}")
land16 = []
for L in [0, 1, 2]:
    d = hyp16[~((hyp16['converted'] == 1) & (hyp16['time_years'] < L))]
    dc, lo, hi, p, n = _boot_dc16(d, 'OIS', 'upsit')
    land16.append(dict(landmark_yr=L, n=n, events=int(d['converted'].sum()),
                       c_upsit=_c16(d, 'upsit'), c_ois=_c16(d, 'OIS'), dc=dc, lo=lo, hi=hi, p=p))
    print(f"{L:<10}{n:>6}{int(d['converted'].sum()):>8}{_c16(d,'upsit'):>10.3f}{_c16(d,'OIS'):>9.3f}"
          f"{dc:>+9.3f}   [{lo:+.3f}, {hi:+.3f}]{p:>9.4f}")
land16 = pd.DataFrame(land16)
sec16['landmark'] = land16.to_dict('records')

print('\n--- §16.1b Time-dependent AUC inside the hyposmia stratum ---')


def _td_auc16(d, score, tp):
    cases = (d['converted'] == 1) & (d['time_years'] <= tp)
    controls = d['time_years'] > tp
    m = cases | controls
    if cases.sum() < 5 or controls.sum() < 5:
        return np.nan, int(cases.sum()), int(m.sum())
    return roc_auc_score(cases[m].astype(int), -d.loc[m, score]), int(cases.sum()), int(m.sum())


td16 = []
print(f"{'Horizon':<10}{'N at risk':>11}{'events<=T':>11}{'AUC(UPSIT)':>12}{'AUC(OIS)':>10}{'dAUC':>9}")
for tp in [1, 2, 3]:
    au, ev, n_at = _td_auc16(hyp16, 'upsit', tp)
    ao, _, _ = _td_auc16(hyp16, 'OIS', tp)
    td16.append(dict(horizon_yr=tp, n_at_risk=n_at, events=ev, auc_upsit=au, auc_ois=ao, d_auc=ao - au))
    print(f"{str(tp)+' yr':<10}{n_at:>11}{ev:>11}{au:>12.3f}{ao:>10.3f}{ao-au:>+9.3f}")
td16 = pd.DataFrame(td16)
sec16['td_auc'] = td16.to_dict('records')

# Fig 2c (v4 numbering): time-dependent AUC in the hyposmia stratum
fig, ax = plt.subplots(figsize=(6.2, 4.4))
xs = np.arange(len(td16)); w = 0.36
ax.bar(xs - w/2, td16['auc_upsit'], w, color=COL_3W['upsit'], edgecolor='black', label='UPSIT (raw)')
ax.bar(xs + w/2, td16['auc_ois'], w, color=COL_3W['OIS'], edgecolor='black', label='OIS')
for x, r in zip(xs, td16.itertuples()):
    ax.text(x, max(r.auc_upsit, r.auc_ois) + 0.02, f'{r.d_auc:+.3f}', ha='center', fontsize=9.5,
            family='monospace', fontweight='bold')
ax.axhline(0.5, color='grey', lw=0.8, ls='--')
ax.set_xticks(xs)
ax.set_xticklabels([f'{int(r.horizon_yr)} yr\n(n={r.n_at_risk}, {r.events} events)' for r in td16.itertuples()],
                   fontsize=9.5)
ax.set_ylim(0.4, 1.0); ax.set_ylabel('Time-dependent AUC', fontsize=11)
ax.set_title('Fig 2c. OIS leads UPSIT at 1 to 2 year horizons (hyposmia stratum)',
             fontsize=10.5, loc='left', fontweight='bold')
ax.legend(frameon=False, fontsize=9.5, loc='upper left')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Fig2c_v4_tdAUC_hyposmia.png', dpi=200)
plt.savefig(f'{FIG_DIR}/Fig2c_v4_tdAUC_hyposmia.pdf')
plt.show()

# Fig S8: landmark sensitivity forest
fig, ax = plt.subplots(figsize=(6.6, 3.4))
ys = np.arange(len(land16))[::-1]
ax.errorbar(land16['dc'], ys, xerr=[land16['dc'] - land16['lo'], land16['hi'] - land16['dc']],
            fmt='o', color=COL_3W['OIS'], capsize=5, ms=7, lw=1.6)
ax.axvline(0, color='black', lw=0.9, ls='--')
ax.set_yticks(ys)
ax.set_yticklabels([f'Landmark {int(r.landmark_yr)} yr\nn={r.n}, {r.events} events' for r in land16.itertuples()],
                   fontsize=9.5)
for y, r in zip(ys, land16.itertuples()):
    ax.text(r.hi + 0.012, y, f'dC = {r.dc:+.3f}, p = {_pfmt(r.p)}', va='center', fontsize=9, family='monospace')
ax.set_xlim(-0.05, 0.45)
ax.set_xlabel('dC-index (OIS minus UPSIT), hyposmia stratum', fontsize=10.5)
ax.set_title('Fig S8. Landmark sensitivity: the OIS advantage is not detection of imminent diagnosis',
             fontsize=10, loc='left', fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/FigS8_v4_landmark.png', dpi=200)
plt.savefig(f'{FIG_DIR}/FigS8_v4_landmark.pdf')
plt.show()
print('\n-> Fig2c_v4_tdAUC_hyposmia + FigS8_v4_landmark saved')
sec16['centre_split'] = CENTRE_SPLIT

# ============================================================================
# Section 16.2. OIS against Cox models trained on conversion labels
# ============================================================================
# Section 16.2 (R2.3): OIS versus fully supervised DaTscan Cox models.
# Every opponent is fit by 5-fold CV INSIDE the prodromal hyposmia cohort with access to
# conversion labels (home advantage); OIS is zero-shot from PD+HC with UPSIT as target.
# Covariates are aligned to the OIS input set, and the decisive control strips age/sex/edu
# from OIS while leaving them with the opponent.

# OIS variant trained on SBR only (no covariates)
_m_tr = train['upsit'].notna()
_sc_sbr = StandardScaler().fit(train.loc[_m_tr, SBR_COLS].values)
_ridge_sbr = Ridge(alpha=1.0).fit(_sc_sbr.transform(train.loc[_m_tr, SBR_COLS].values),
                                  train.loc[_m_tr, 'upsit'].values)
df['OIS_sbr'] = _ridge_sbr.predict(_sc_sbr.transform(df[SBR_COLS].values))
hyp16['OIS_sbr'] = hyp16['PATNO'].map(df.set_index('PATNO')['OIS_sbr'])

# Binary DAT-deficit flag (PARS rule applied to the validated PPMI %-expected metric)
hyp16['dat_deficit'] = (pct_expected_putamen(hyp16) < PARS_DAT_CUT).astype(float)
hyp16['neg_flag'] = -hyp16['dat_deficit']

COV16 = ['age', 'sex_num', 'educyrs']
h2h = hyp16.dropna(subset=['OIS', 'OIS_sbr', 'upsit', 'dat_deficit'] + SBR_COLS + COV16).copy()
print(f'§16.2 head-to-head cohort: N={len(h2h)}, events={int(h2h["converted"].sum())}')


def _coxnet_oof(d, feats, seed=41, l1_ratio=0.05, alpha=0.1):
    # 5-fold CV out-of-fold risk predictions from an elastic-net Cox model (handles collinear SBR).
    X = d[feats].values.astype(float)
    y = Surv.from_arrays(event=d['converted'].astype(bool).values, time=d['time_years'].values)
    oof = np.full(len(d), np.nan); n_fail = 0
    for tr_i, te_i in KFold(n_splits=5, shuffle=True, random_state=seed).split(X):
        sc = StandardScaler().fit(X[tr_i])
        try:
            m = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[alpha], fit_baseline_model=False)
            m.fit(sc.transform(X[tr_i]), y[tr_i])
            oof[te_i] = m.predict(sc.transform(X[te_i]))
        except Exception as e:
            n_fail += 1
            print('   fold failed:', e)
    return oof, n_fail


_fails = 0
for _name, _feats in [('put_cov', ['PUTAMEN_REF_CWM'] + COV16),
                      ('sbr33', SBR_COLS),
                      ('sbr33_cov', SBR_COLS + COV16)]:
    _risk, _f = _coxnet_oof(h2h, _feats)
    h2h['neg_risk_' + _name] = -_risk        # Coxnet returns risk (higher = worse), flip to protective sign
    _fails += _f
print(f'  Coxnet fold failures: {_fails} / 15')

MODELS16 = [
    ('upsit',            'UPSIT (raw)',            'behavioural baseline'),
    ('neg_flag',         'Binary DAT-deficit flag', 'threshold rule'),
    ('neg_risk_put_cov', 'Putamen + age/sex/edu',  'Coxnet 5-fold CV (in-cohort)'),
    ('neg_risk_sbr33',   '33-SBR only',            'Coxnet 5-fold CV (in-cohort)'),
    ('neg_risk_sbr33_cov', '33-SBR + age/sex/edu', 'Coxnet 5-fold CV (in-cohort)'),
    ('OIS_sbr',          'OIS, SBR-only',          'Ridge on UPSIT, PD+HC zero-shot'),
    ('OIS',              'OIS, full',              'Ridge on UPSIT, PD+HC zero-shot'),
]
print(f"\n{'Model':<24}{'C':>7}{'dC vs OIS':>11}{'95% CI':>21}{'p(<=0)':>9}  training")
h2h_rows = []
for col, lab, how in MODELS16:
    c = _c16(h2h, col)
    if col == 'OIS':
        h2h_rows.append(dict(model=lab, c=c, dc=np.nan, lo=np.nan, hi=np.nan, p=np.nan, training=how))
        print(f'{lab:<24}{c:>7.3f}{"-":>11}{"-":>21}{"-":>9}  {how}')
        continue
    dc, lo, hi, p, _n = _boot_dc16(h2h, 'OIS', col, one_sided=True)
    h2h_rows.append(dict(model=lab, c=c, dc=dc, lo=lo, hi=hi, p=p, training=how))
    print(f'{lab:<24}{c:>7.3f}{dc:>+11.3f}   [{lo:+.3f}, {hi:+.3f}]{p:>9.3f}  {how}')
h2h_tab = pd.DataFrame(h2h_rows)
sec16['head_to_head'] = h2h_tab.to_dict('records')

_dcs, _los, _his, _ps, _ = _boot_dc16(h2h, 'OIS_sbr', 'neg_risk_sbr33_cov', one_sided=True)
sec16['decisive_control'] = dict(dc=_dcs, lo=_los, hi=_his, p=_ps)
print(f'\n  Decisive control, covariates given to the opponent only:')
print(f'    OIS(SBR-only) vs 33-SBR+cov Coxnet: dC = {_dcs:+.3f} [{_los:+.3f}, {_his:+.3f}], p = {_ps:.3f}'
      f'   -> {"statistical tie" if _ps > 0.05 else "OIS ahead"}')

# Fig 2d: head-to-head bars
fig, ax = plt.subplots(figsize=(8.4, 4.8))
_ord = h2h_tab.iloc[::-1].reset_index(drop=True)
ys = np.arange(len(_ord))
_cols = ['#d4624a' if 'OIS' in m else ('#6b6b6b' if 'UPSIT' in m else '#2e3a5c') for m in _ord['model']]
ax.barh(ys, _ord['c'], color=_cols, edgecolor='black', height=0.62)
ax.axvline(0.5, color='grey', lw=0.8, ls='--')
for y, r in zip(ys, _ord.itertuples()):
    txt = f'  C = {r.c:.3f}'
    if not np.isnan(r.dc):
        txt += f'   (dC {r.dc:+.3f}, p = {_pfmt(r.p)})'
    ax.text(r.c + 0.005, y, txt, va='center', fontsize=9, family='monospace')
ax.set_yticks(ys)
ax.set_yticklabels([f'{r.model}\n{r.training}' for r in _ord.itertuples()], fontsize=8.8)
ax.set_xlim(0.45, 1.02)
ax.set_xlabel('C-index for prodromal PD conversion (hyposmia stratum)', fontsize=10.5)
ax.set_title(f'Fig 2d. Zero-shot OIS matches fully supervised in-cohort DaTscan Cox models\n'
             f'N = {len(h2h)}, {int(h2h["converted"].sum())} events; dC vs OIS(full), 1000-iteration bootstrap',
             fontsize=10.5, loc='left', fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Fig2d_v4_headtohead.png', dpi=200)
plt.savefig(f'{FIG_DIR}/Fig2d_v4_headtohead.pdf')
plt.show()
print('-> Fig2d_v4_headtohead saved')

# ============================================================================
# Section 16.3 (validation block). Expected-putamen metric validation into the ledger
# ============================================================================
sec16['dat_metric_validation'] = DAT_VALIDATION
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(sec16, f, indent=2, default=float)
print('-> results/repro/section16_stats.json written (incl. dat_metric_validation)')

# ============================================================================
# Section 16.4. Calibration of OIS in the prodromal cohort
# ============================================================================
# Section 16.4 (R1): does an OIS trained on PD+HC extrapolate to the prodromal cohort?
# Calibration = regress observed UPSIT on predicted OIS. Perfect calibration is slope 1,
# intercept 0, mean residual (OMI) 0. Reported for all prodromal and by hyposmia stratum.
prod_cal = df[df['COHORT'] == 4].dropna(subset=['OIS', 'upsit']).copy()
prod_cal['grp'] = np.where(prod_cal['subgroup'] == 'Hyposmia', 'Hyposmia', 'Non-Hyposmia')

print('§16.4 OIS calibration in the prodromal cohort (never seen during training):')
print(f"{'Cohort':<16}{'N':>6}{'slope':>9}{'intercept':>11}{'r':>8}{'mean OMI':>10}")
cal16 = []
for lab, d in [('All prodromal', prod_cal),
               ('Non-Hyposmia', prod_cal[prod_cal['grp'] == 'Non-Hyposmia']),
               ('Hyposmia', prod_cal[prod_cal['grp'] == 'Hyposmia'])]:
    sl, ic, r, p, _se = stats.linregress(d['OIS'], d['upsit'])
    cal16.append(dict(cohort=lab, n=len(d), slope=sl, intercept=ic, r=r, p=p, mean_omi=d['OMI'].mean()))
    print(f'{lab:<16}{len(d):>6}{sl:>9.3f}{ic:>11.2f}{r:>8.3f}{d["OMI"].mean():>10.2f}')
cal16 = pd.DataFrame(cal16)
sec16['calibration'] = cal16.to_dict('records')
print('  Non-Hyposmia slope close to 1 means no systematic scale drift when extrapolating')
print('  from PD+HC to prodromal; the Hyposmia slope below 1 is the UPSIT floor, not model failure.')

# Decile calibration in Non-Hyposmia (the layer where UPSIT variance is preserved)
nh_cal = prod_cal[prod_cal['grp'] == 'Non-Hyposmia'].copy()
nh_cal['dec'] = pd.qcut(nh_cal['OIS'], 10, labels=False) + 1
dec16 = nh_cal.groupby('dec').agg(n=('OIS', 'size'), pred=('OIS', 'mean'), obs=('upsit', 'mean')).reset_index()
dec16['residual'] = dec16['obs'] - dec16['pred']
print('\n  Decile calibration (Non-Hyposmia):')
print(dec16.to_string(index=False, float_format=lambda v: f'{v:.2f}'))
sec16['calibration_deciles'] = dec16.to_dict('records')

# Sup Fig S1: decile calibration plot
fig, ax = plt.subplots(figsize=(5.6, 5.2))
lo = min(dec16['pred'].min(), dec16['obs'].min()) - 2
hi = max(dec16['pred'].max(), dec16['obs'].max()) + 2
ax.plot([lo, hi], [lo, hi], color='lightgrey', lw=1.2, ls='--', label='perfect calibration')
ax.plot(dec16['pred'], dec16['obs'], 'o-', color=COL_3W['OIS'], ms=7, lw=1.5, label='observed decile means')
_sl_nh = cal16.loc[cal16['cohort'] == 'Non-Hyposmia', 'slope'].iloc[0]
_r_nh = cal16.loc[cal16['cohort'] == 'Non-Hyposmia', 'r'].iloc[0]
ax.set_xlabel('Predicted OIS (decile mean)', fontsize=11)
ax.set_ylabel('Observed UPSIT (decile mean)', fontsize=11)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_title(f'Fig S1. OIS calibration in non-hyposmic prodromals\nslope = {_sl_nh:.3f}, r = {_r_nh:.3f}, n = {len(nh_cal)}',
             fontsize=10.5, loc='left', fontweight='bold')
ax.legend(frameon=False, fontsize=9.5, loc='upper left')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/FigS1_v4_calibration_decile.png', dpi=200)
plt.savefig(f'{FIG_DIR}/FigS1_v4_calibration_decile.pdf')
plt.show()
print('\n-> FigS1_v4_calibration_decile saved')

# ============================================================================
# Section 16.5. PARS high-risk subgroup
# ============================================================================
# Section 16.5 (R2.4): the PARS extension, tested inside the PARS-defined high-risk layer
# (hyposmia AND putamen DAT deficit). This is the subgroup PARS 2017 reported at ~67%
# four-year conversion, i.e. exactly where a continuous score has to earn its place.
pars16 = hyp16[hyp16['dat_deficit'] == 1].copy()
_dc, _lo, _hi, _p, _n = _boot_dc16(pars16, 'OIS', 'upsit')
print(f'§16.5 PARS high-risk layer (hyposmia + putamen DAT deficit): N={_n}, events={int(pars16["converted"].sum())}')
print(f'  C(UPSIT) = {_c16(pars16, "upsit"):.3f}  ->  C(OIS) = {_c16(pars16, "OIS"):.3f}')
print(f'  dC = {_dc:+.3f}  95% CI [{_lo:+.3f}, {_hi:+.3f}]  p = {_p:.4f}')
sec16['pars_extension'] = dict(n=int(_n), events=int(pars16['converted'].sum()),
                               c_upsit=_c16(pars16, 'upsit'), c_ois=_c16(pars16, 'OIS'),
                               dc=_dc, lo=_lo, hi=_hi, p=_p)

with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(sec16, f, indent=2, default=float)
print('\n-> results/repro/section16_stats.json rewritten with calibration + PARS extension')

# ============================================================================
# Section 17. Variance decomposition of the UPSIT total
# ============================================================================
# Section 17 (R1, motivation): the raw UPSIT score is a mixture. Splitting it against the
# dopaminergic signal separates a small component that carries the conversion signal from a
# large residual that carries none. This is the quantitative basis for the whole decomposition.
prod['grp'] = np.where(prod['subgroup'] == 'Hyposmia', 'Hyposmia', 'Non-Hyposmia')


def _c17(d, col):
    v = d[[col, 'time_years', 'converted']].dropna()
    return concordance_index(v['time_years'], v[col], v['converted']), len(v), int(v['converted'].sum())


def _boot_c17(d, col, n_boot=1000, seed=42):
    v = d[[col, 'time_years', 'converted']].dropna()
    rng = np.random.RandomState(seed); idx = np.arange(len(v)); out = []
    for _ in range(n_boot):
        vb = v.iloc[rng.choice(idx, len(idx), replace=True)]
        if vb['converted'].sum() < 3:
            continue
        try:
            out.append(concordance_index(vb['time_years'], vb[col], vb['converted']))
        except Exception:
            pass
    out = np.array(out)
    return np.percentile(out, 2.5), np.percentile(out, 97.5)


print('§17 Decomposition of the smell score against the dopaminergic signal')
print('   UPSIT = OIS + OMI by construction. C-index uses the protective sign convention.\n')
dec17 = []
_DEC_COLS = ['upsit', 'OIS', 'OMI', 'PUTAMEN_REF_CWM', 'time_years', 'converted']
for lab, d in [('All prodromal', prod), ('Hyposmia', prod[prod['grp'] == 'Hyposmia']),
               ('Non-Hyposmia', prod[prod['grp'] != 'Hyposmia'])]:
    # One common complete-case set for all four scores, so every row of the table shares an n
    # and matches the paired comparisons in section 33. Scoring each column on its own
    # complete-case set would silently change the denominator between rows.
    d = d.dropna(subset=_DEC_COLS)
    v = d[['upsit', 'OIS', 'OMI']].dropna()
    sd_u, sd_o, sd_m = v['upsit'].std(), v['OIS'].std(), v['OMI'].std()
    r_om = stats.pearsonr(v['OIS'], v['OMI'])[0]
    print(f'  {lab}  (n={len(v)})')
    print(f'    SD: UPSIT {sd_u:.2f} | OIS {sd_o:.2f} | OMI {sd_m:.2f}   '
          f'(the residual carries {sd_m**2/sd_u**2*100:.0f}% of the score variance, r(OIS,OMI)={r_om:+.3f})')
    row = dict(cohort=lab, n=len(v), sd_upsit=sd_u, sd_ois=sd_o, sd_omi=sd_m, r_ois_omi=r_om)
    for col, nm in [('upsit', 'UPSIT (raw)'), ('OIS', 'OIS'), ('OMI', 'OMI'), ('PUTAMEN_REF_CWM', 'Putamen SBR')]:
        c, n, e = _c17(d, col)
        lo, hi = _boot_c17(d, col)
        row[f'c_{col}'] = c; row[f'c_{col}_lo'] = lo; row[f'c_{col}_hi'] = hi
        flag = '   <- at chance' if lo < 0.5 < hi else ''
        print(f'    C({nm:12s}) = {c:.3f}  [{lo:.3f}, {hi:.3f}]  (n={n}, {e} events){flag}')
    dec17.append(row)
    print()
dec17 = pd.DataFrame(dec17)
sec16['decomposition'] = dec17.to_dict('records')
print('  Reading: the component that carries almost all the conversion signal (OIS) is the smaller')
print('  part of the score, and the larger part (OMI) is at chance for conversion. Raw UPSIT, being')
print('  their sum, lands in between. This is why the smell test underperforms where it is applied.')

# Figure: decomposition of the score in the hyposmia stratum
hyp17 = prod[prod['grp'] == 'Hyposmia']
fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.4, 4.3))
row_h = dec17[dec17['cohort'] == 'Hyposmia'].iloc[0]
labs = ['UPSIT (raw)', 'OIS', 'OMI']
cs = [row_h['c_upsit'], row_h['c_OIS'], row_h['c_OMI']]
los = [row_h['c_upsit_lo'], row_h['c_OIS_lo'], row_h['c_OMI_lo']]
his = [row_h['c_upsit_hi'], row_h['c_OIS_hi'], row_h['c_OMI_hi']]
cols = [COL_3W['upsit'], COL_3W['OIS'], COL_3W['OMI']]
axA.bar(np.arange(3), cs, yerr=[np.array(cs) - np.array(los), np.array(his) - np.array(cs)],
        capsize=6, color=cols, edgecolor='black', width=0.6)
axA.axhline(0.5, color='grey', lw=0.9, ls='--')
axA.text(2.42, 0.505, 'chance', fontsize=8.5, color='grey', va='bottom', ha='right')
for x, c in zip(np.arange(3), cs):
    axA.text(x, c + 0.02, f'{c:.3f}', ha='center', fontsize=10, family='monospace', fontweight='bold')
axA.set_xticks(np.arange(3)); axA.set_xticklabels(labs, fontsize=10)
axA.set_ylim(0.40, 0.88); axA.set_ylabel('C-index for PD conversion', fontsize=10.5)
axA.set_title('a  Conversion signal by component', fontsize=10.5, loc='left', fontweight='bold')

sds = [row_h['sd_upsit'], row_h['sd_ois'], row_h['sd_omi']]
axB.bar(np.arange(3), sds, color=cols, edgecolor='black', width=0.6)
for x, v_ in zip(np.arange(3), sds):
    axB.text(x, v_ + 0.12, f'{v_:.2f}', ha='center', fontsize=10, family='monospace', fontweight='bold')
axB.set_xticks(np.arange(3)); axB.set_xticklabels(labs, fontsize=10)
axB.set_ylabel('SD (UPSIT points)', fontsize=10.5)
axB.set_title('b  Where the variance sits', fontsize=10.5, loc='left', fontweight='bold')
fig.suptitle(f'The smell score is a mixture, hyposmic prodromals (n={len(hyp17)}, '
             f'{int(hyp17["converted"].sum())} events)', fontsize=11, fontweight='bold', x=0.01, ha='left')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Fig2_v4_decomposition.png', dpi=200)
plt.savefig(f'{FIG_DIR}/Fig2_v4_decomposition.pdf')
plt.show()
print('\n-> Fig2_v4_decomposition saved')
merge_ledger(sec16)

# ============================================================================
# State for the next scripts
# ============================================================================
save_state('02_cohort', df=df, train=train, excel=excel, sbr=sbr, SBR_COLS=SBR_COLS, FEATURES=FEATURES, DEMO=DEMO,
           prod=prod, prod_with_prevalent=prod_with_prevalent, prevalent_at_baseline=prevalent_at_baseline,
           first_visit=first_visit, last_visit=last_visit, pd_diag=pd_diag, hyp16=hyp16,
           NORM_BETA=NORM_BETA, age_col=age_col, PARS_DAT_CUT=PARS_DAT_CUT, DAT_VALIDATION=DAT_VALIDATION,
           TEMPORAL_SPLIT=TEMPORAL_SPLIT, ICC_STATS=ICC_STATS, CENTRE_SPLIT=CENTRE_SPLIT)
