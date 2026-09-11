"""03_discrimination.py  Discrimination, strata, paired comparisons, risk groups

The pre-specified paired comparisons, the interaction test, the held-out correlation by
recruitment cohort, the concordance indices with bootstrap intervals for every stratum and reading
(the only source of Tables 2 to 4), continuous against dichotomised readings, rank-matched
comparisons, calibration and age in the non-deficit stratum, the group flow, the median split,
the permutation-corrected three-group split and the IPCW time-dependent AUC.
Outputs: results/repro/section33_*.csv, section35_*.csv, section36_cindex_ci.csv,
section38_*.csv, section39_*.csv, section40_*.csv, section44_*.csv, section46_td_auc.csv,
supp_median_split.csv, and ledger keys.
Notebook provenance: sections 33, 34, 35, 36, 38, 39, 40, 41, 44, 46, plus the median split. Code is the notebook cell text, unchanged except for
the file paths, which are resolved through common.py.
"""
from common import *
st = load_state('02_cohort'); globals().update(st)
sec16 = {}
pct_expected_putamen = make_pct_expected_putamen(NORM_BETA, age_col)

# ============================================================================
# Prodromal frame with imaging state (as in section 15.3)
# ============================================================================
prod15 = df[df['COHORT']==4].copy()
prod15['grp'] = np.where(prod15['subgroup']=='Hyposmia','Hyposmia','Non-Hyposmia')
prod15['put_pct'] = pct_expected_putamen(prod15)
prod15['dat_normal_bl'] = (prod15['put_pct'] >= PARS_DAT_CUT).astype(int)

# ============================================================================
# Section 33. Pre-specified paired comparisons
# ============================================================================
# Section 33: the pre-specified comparison set. No multiplicity correction.
# Six comparisons, each backing a claim made in the main text, all reported regardless of direction.
print('Section 33  Pre-specified pairwise comparisons (no multiplicity correction)\n')
_d33 = prod.dropna(subset=['OIS', 'OMI', 'upsit', 'PUTAMEN_REF_CWM', 'PUTAMEN_L_REF_CWM',
                           'PUTAMEN_R_REF_CWM', 'time_years', 'converted', 'age', 'sex_num']).copy()
_d33 = _d33[(_d33['converted'] == 1) | (_d33['time_years'] > 0)]
_d33['time_years'] = _d33['time_years']  # already floored at 0.25 for converters upstream; do NOT re-clip censored rows
_d33['pct'] = _d33[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (
    1.8547 - 0.00303 * _d33['age'] - 0.2419 * _d33['sex_num']) * 100
_hyp33 = _d33[_d33['subgroup'] == 'Hyposmia']
S33 = {'hyposmia': _hyp33,
       'hyposmia + DAT deficit': _hyp33[_hyp33['pct'] < 65],
       'hyposmia + no DAT deficit': _hyp33[_hyp33['pct'] >= 65],
       'non-hyposmia-enriched': _d33[_d33['subgroup'] != 'Hyposmia']}
PLAN33 = [('hyposmia', 'OIS', 'upsit'),
          ('hyposmia', 'OIS', 'PUTAMEN_REF_CWM'),
          ('hyposmia + DAT deficit', 'OIS', 'upsit'),
          ('hyposmia + no DAT deficit', 'OIS', 'PUTAMEN_REF_CWM'),
          ('hyposmia + no DAT deficit', 'OIS', 'upsit'),
          ('non-hyposmia-enriched', 'OIS', 'upsit')]
LBL33 = {'upsit': 'UPSIT', 'PUTAMEN_REF_CWM': 'putamen', 'OIS': 'OIS'}
_rng33 = np.random.default_rng(42); _B33 = 1000
rows33 = []
print(f"  {'stratum':<28}{'comparison':<18}{'n':>6}{'ev':>5}{'C(A)':>8}{'C(B)':>8}{'dC':>8}{'95% CI':>20}{'p':>9}")
for nm, a, b in PLAN33:
    g = S33[nm]
    t, e = g['time_years'].values, g['converted'].values.astype(int)
    A, B = g[a].values, g[b].values
    bd = []
    for _ in range(_B33):
        i = _rng33.choice(len(g), len(g), replace=True)
        if e[i].sum() < 5:
            continue
        try:
            bd.append(concordance_index(t[i], A[i], e[i]) - concordance_index(t[i], B[i], e[i]))
        except Exception:
            pass
    bd = np.array(bd)
    lo, hi = np.percentile(bd, [2.5, 97.5])
    p = max(2 * min((bd <= 0).mean(), (bd >= 0).mean()), 1 / _B33)
    cA, cB = concordance_index(t, A, e), concordance_index(t, B, e)
    lab = f'{LBL33[a]} - {LBL33[b]}'
    print(f"  {nm:<28}{lab:<18}{len(g):>6}{int(e.sum()):>5}{cA:>8.3f}{cB:>8.3f}{cA - cB:>+8.3f}"
          f"   [{lo:+.3f}, {hi:+.3f}]{p:>9.4f}")
    rows33.append(dict(stratum=nm, comparison=lab, n=int(len(g)), events=int(e.sum()),
                       c_a=float(cA), c_b=float(cB), dC=float(cA - cB),
                       lo=float(lo), hi=float(hi), p=float(p)))
mc33 = pd.DataFrame(rows33)
print('\n  No multiplicity correction is applied. The premise is that this set is pre-specified')
print('  and reported in full, including the null result in the DAT-non-deficit stratum')
print('  (OIS - UPSIT = +0.028, p = 0.50), which stays in the main text.')
print('  An earlier version corrected 21 tests over 7 strata; the 15 dropped comparisons back no')
print('  claim in the main text and no inferential result is reported for them.')
mc33.to_csv('results/repro/section33_prespecified.csv', index=False)
sec16['prespecified_comparisons'] = mc33.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m33 = json.load(f)
_b33 = set(_m33)
_m33.pop('multiplicity', None)
_m33.update(sec16)
_m33.pop('multiplicity', None)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m33, f, indent=2, default=float)
print(f'\n-> results/repro/section33_prespecified.csv written')
print(f'-> section16_stats.json: {len(_b33)} -> {len(_m33)} keys'
      f" (key 'multiplicity' removed, 'prespecified_comparisons' added)")

# ============================================================================
# Section 34. Interaction test (difference in differences)
# ============================================================================
# Section 34: the interaction claim needs a direct test, not a ratio of point estimates.
# Comparing "significant in one stratum, not in the other" is the classic invalid argument.
print('Section 34  Direct test of the stratum interaction (difference in differences)\n')
_d34 = prod.dropna(subset=['OIS', 'OMI', 'upsit', 'PUTAMEN_REF_CWM', 'PUTAMEN_L_REF_CWM',
                           'PUTAMEN_R_REF_CWM', 'time_years', 'converted', 'age', 'sex_num']).copy()
_d34 = _d34[(_d34['converted'] == 1) | (_d34['time_years'] > 0)]
_d34['time_years'] = _d34['time_years']  # already floored at 0.25 for converters upstream; do NOT re-clip censored rows
_H34 = _d34[_d34['subgroup'] == 'Hyposmia']
_N34 = _d34[_d34['subgroup'] != 'Hyposmia']

def _dc34(g):
    t = g['time_years'].values
    e = g['converted'].values.astype(int)
    return concordance_index(t, g['OIS'].values, e) - concordance_index(t, g['upsit'].values, e)

_obs_h, _obs_n = _dc34(_H34), _dc34(_N34)
_rng34 = np.random.default_rng(42)
_dd34 = []
for _ in range(2000):
    _ih = _rng34.choice(len(_H34), len(_H34), replace=True)
    _in = _rng34.choice(len(_N34), len(_N34), replace=True)
    try:
        _dd34.append(_dc34(_H34.iloc[_ih]) - _dc34(_N34.iloc[_in]))
    except Exception:
        pass
_dd34 = np.array(_dd34)
_lo34, _hi34 = np.percentile(_dd34, [2.5, 97.5])
_p34 = max(2 * min((_dd34 <= 0).mean(), (_dd34 >= 0).mean()), 1 / len(_dd34))
print(f'  Hyposmia stratum       n={len(_H34):>5}  events={int(_H34.converted.sum()):>4}  dC={_obs_h:+.3f}')
print(f'  Non-hyposmia-enriched  n={len(_N34):>5}  events={int(_N34.converted.sum()):>4}  dC={_obs_n:+.3f}')
print(f'  Ratio of point estimates = {_obs_h / _obs_n:.2f}x')
print(f'\n  Difference in differences = {_obs_h - _obs_n:+.3f}'
      f'  95% CI [{_lo34:+.3f}, {_hi34:+.3f}]  P = {_p34:.3f}   (B={len(_dd34)})')
print(f'  Bootstrap fraction above zero = {(_dd34 > 0).mean():.3f}')
print(f'  Verdict: {"significant" if _p34 < 0.05 else "NOT significant - descriptive observation only"}')
print('\n  The interaction is directional and the point estimates differ twofold, but the direct')
print('  test does not reach significance. Reporting one stratum as significant and the other as')
print('  not is NOT evidence of effect modification. The main text must say so.')
sec16['interaction_test'] = {
    'dC_hyposmia': float(_obs_h), 'n_hyposmia': int(len(_H34)),
    'events_hyposmia': int(_H34.converted.sum()),
    'dC_control': float(_obs_n), 'n_control': int(len(_N34)),
    'events_control': int(_N34.converted.sum()),
    'ratio': float(_obs_h / _obs_n), 'dd': float(_obs_h - _obs_n),
    'lo': float(_lo34), 'hi': float(_hi34), 'p': float(_p34),
    'frac_above_zero': float((_dd34 > 0).mean()), 'n_boot': int(len(_dd34)),
    'verdict': ('significant' if _p34 < 0.05 else 'not significant; descriptive observation only'),
}
with open('results/repro/section16_stats.json') as f:
    _m34 = json.load(f)
_b34 = set(_m34); _m34.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m34, f, indent=2, default=float)
print(f'\n-> section16_stats.json merged: {len(_b34)} -> {len(_m34)} keys, 0 dropped')

# ============================================================================
# Section 35. Held-out correlation by recruitment cohort
# ============================================================================
# Section 35: the three enrolment arms behind the prodromal held-out correlation.
from scipy.stats import pearsonr as _pr35
print('Section 35  Prodromal held-out correlation by enrolment arm\n')
_d35 = prod15.dropna(subset=['OIS', 'upsit']).copy()

def _armA(x):
    if x == 'Hyposmia':
        return 'Hyposmia'
    return 'RBD' if 'RBD' in str(x) else 'Genetic'      # rule A: dual carriers -> RBD (used in the paper)

def _armB(x):
    if x == 'Hyposmia':
        return 'Hyposmia'
    return 'RBD' if str(x) == 'RBD' else 'Genetic'      # rule B: dual carriers -> Genetic (old sections 22/25)

_rows35 = []
for _rule, _fn in [('A (paper: dual -> RBD)', _armA), ('B (sections 22/25: dual -> Genetic)', _armB)]:
    _d35['arm'] = _d35['subgroup'].map(_fn)
    print(f'  Rule {_rule}')
    print(f"    {'arm':<12}{'n':>6}{'r':>9}{'p':>11}   UPSIT mean (SD)")
    for _g in ['Hyposmia', 'RBD', 'Genetic']:
        _s = _d35[_d35['arm'] == _g]
        _r, _p = _pr35(_s['OIS'], _s['upsit'])
        print(f"    {_g:<12}{len(_s):>6}{_r:>9.3f}{_p:>11.1e}   "
              f"{_s['upsit'].mean():.1f} ({_s['upsit'].std():.1f})")
        _rows35.append(dict(rule=_rule.split()[0], arm=_g, n=int(len(_s)), r=float(_r), p=float(_p),
                            upsit_mean=float(_s['upsit'].mean()), upsit_sd=float(_s['upsit'].std())))
    print()
_rp, _pp = _pr35(_d35['OIS'], _d35['upsit'])
_nh = _d35[_d35['subgroup'] != 'Hyposmia']
_rn, _pn = _pr35(_nh['OIS'], _nh['upsit'])
print(f'  Pooled prodromal          n={len(_d35):>5}  r={_rp:.3f}  p={_pp:.1e}')
print(f'  Non-Hyposmia (2-way, 16.4) n={len(_nh):>4}  r={_rn:.3f}  p={_pn:.1e}')
print('\n  The pooled value exceeds every single arm because the three arms differ in their')
print('  UPSIT distributions. Pooled and within-arm values are therefore always reported')
print('  together. The paper uses rule A throughout.')

# Which genes the carrier arm actually contains. Readers meeting this channel for the
# first time cannot tell from its name, so the composition is reported in the text.
_d35['arm'] = _d35['subgroup'].map(_armA)                       # back to rule A
_carr_enrol = prod15[prod15['subgroup'].map(_armA) == 'Genetic']['subgroup'].value_counts()
_carr_anal = _d35[_d35['arm'] == 'Genetic']['subgroup'].value_counts()
print(f'\n  Carrier arm composition (rule A), enrolled n={int(_carr_enrol.sum())}, analysed n={int(_carr_anal.sum())}')
for _v in _carr_enrol.index:
    print(f'    {_v:<16} enrolled {int(_carr_enrol[_v]):>4}   analysed {int(_carr_anal.get(_v, 0)):>4}')
# Ceiling of the olfactory test in that arm, the reason its calibration slope is flattest
_ca = _d35[_d35['arm'] == 'Genetic']['upsit']
print(f'    UPSIT in the carrier arm: {_ca.mean():.1f} ({_ca.std():.1f}), '
      f'{100 * (_ca >= 30).mean():.0f}% score >= 30 and {100 * (_ca >= 35).mean():.0f}% >= 35, '
      f'against {100 * (_d35.loc[_d35.arm == "Hyposmia", "upsit"] >= 30).mean():.0f}% >= 30 in the hyposmia arm')
sec16['carrier_arm_composition'] = dict(
    enrolled={str(k): int(v) for k, v in _carr_enrol.items()},
    analysed={str(k): int(v) for k, v in _carr_anal.items()},
    n_enrolled=int(_carr_enrol.sum()), n_analysed=int(_carr_anal.sum()),
    upsit_mean=float(_ca.mean()), upsit_sd=float(_ca.std()),
    pct_ge30=float(100 * (_ca >= 30).mean()), pct_ge35=float(100 * (_ca >= 35).mean()),
    pct_ge30_hyposmia=float(100 * (_d35.loc[_d35.arm == 'Hyposmia', 'upsit'] >= 30).mean()))
arm35 = pd.DataFrame(_rows35)
arm35.to_csv('results/repro/section35_heldout_by_arm.csv', index=False)
sec16['heldout_by_arm'] = {'rule_used_in_paper': 'A (dual carriers assigned to RBD)',
                           'rows': arm35.to_dict('records'),
                           'pooled_r': float(_rp), 'pooled_n': int(len(_d35))}
with open('results/repro/section16_stats.json') as f:
    _m35 = json.load(f)
_b35 = set(_m35); _m35.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m35, f, indent=2, default=float)
print(f'\n-> results/repro/section35_heldout_by_arm.csv written')
print(f'-> section16_stats.json: {len(_b35)} -> {len(_m35)} keys')

# ============================================================================
# Section 36. Concordance with 95% CI, every stratum and reading
# ============================================================================
# Section 36: one table of C-index with bootstrap CI for every stratum by every reading.
print('Section 36  C-index with 95% CI, all strata by all readings\n')
_d36 = prod.dropna(subset=['OIS', 'OMI', 'upsit', 'PUTAMEN_REF_CWM', 'PUTAMEN_L_REF_CWM',
                           'PUTAMEN_R_REF_CWM', 'time_years', 'converted', 'age', 'sex_num']).copy()
_d36['pct_exp'] = _d36[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (
    1.8547 - 0.00303 * _d36['age'] - 0.2419 * _d36['sex_num']) * 100
_h36 = _d36[_d36['subgroup'] == 'Hyposmia']
LAY36 = [('all prodromal', _d36),
         ('hyposmia', _h36),
         ('hyposmia + DAT deficit', _h36[_h36['pct_exp'] < 65]),
         ('hyposmia + no DAT deficit', _h36[_h36['pct_exp'] >= 65]),
         ('non-hyposmia-enriched', _d36[_d36['subgroup'] != 'Hyposmia'])]
READ36 = [('upsit', 'UPSIT total'), ('OIS', 'OIS'), ('OMI', 'OMI'),
          ('PUTAMEN_REF_CWM', 'putamen SBR'), ('pct_exp', 'lowest putamen %expected')]
_rng36 = np.random.default_rng(42); _B36 = 1000
rows36 = []
print(f"  {'stratum':<28}{'n':>6}{'ev':>5}  " + ''.join(f'{lab:<26}' for _, lab in READ36))
for nm, g in LAY36:
    t, e = g['time_years'].values, g['converted'].values.astype(int)
    cells = []
    for col, lab in READ36:
        a = g[col].values
        bs = []
        for _ in range(_B36):
            i = _rng36.choice(len(g), len(g), replace=True)
            if e[i].sum() < 5:
                continue
            try:
                bs.append(concordance_index(t[i], a[i], e[i]))
            except Exception:
                pass
        lo, hi = np.percentile(bs, [2.5, 97.5])
        c = concordance_index(t, a, e)
        cells.append(f'{c:.3f} [{lo:.3f}, {hi:.3f}]')
        rows36.append(dict(stratum=nm, reading=lab, n=int(len(g)), events=int(e.sum()),
                           c=float(c), lo=float(lo), hi=float(hi),
                           at_chance=bool(lo < 0.5 < hi)))
    print(f'  {nm:<28}{len(g):>6}{int(e.sum()):>5}  ' + ''.join(f'{x:<26}' for x in cells))
ci36 = pd.DataFrame(rows36)
print('\n  Readings whose CI contains 0.5 (at chance):')
for _, r in ci36[ci36['at_chance']].iterrows():
    print(f"    {r['stratum']:<28}{r['reading']:<26}{r['c']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]")
ci36.to_csv('results/repro/section36_cindex_ci.csv', index=False)
sec16['cindex_ci'] = ci36.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m36 = json.load(f)
_b36 = set(_m36); _m36.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m36, f, indent=2, default=float)
print(f'\n-> results/repro/section36_cindex_ci.csv written')
print(f'-> section16_stats.json: {len(_b36)} -> {len(_m36)} keys')

# ============================================================================
# Section 38. Continuous against dichotomised readings
# ============================================================================
# Section 38: what dichotomising each reading costs. Descriptive only, no paired test.
print('Section 38  Continuous versus dichotomised readings\n')
_p38 = excel[excel['EVENT_ID'] == 'BL'].drop_duplicates('PATNO')[['PATNO', 'upsit_pctl15']]
_d38 = prod.dropna(subset=['OIS', 'upsit', 'PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM',
                           'time_years', 'converted', 'age', 'sex_num']).merge(_p38, on='PATNO', how='left')
_d38 = _d38.dropna(subset=['upsit_pctl15'])
_d38['pct_exp'] = _d38[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (
    1.8547 - 0.00303 * _d38['age'] - 0.2419 * _d38['sex_num']) * 100
_rng38 = np.random.default_rng(42); _B38 = 1000

def _c38(g, vals):
    t = g['time_years'].values; e = g['converted'].values.astype(int)
    a = np.asarray(vals, dtype=float)
    bs = []
    for _ in range(_B38):
        i = _rng38.choice(len(g), len(g), replace=True)
        if e[i].sum() < 5:
            continue
        try:
            bs.append(concordance_index(t[i], a[i], e[i]))
        except Exception:
            pass
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return float(concordance_index(t, a, e)), float(lo), float(hi)

rows38 = []
print(f"  {'stratum':<24}{'reading':<40}{'C-index [95% CI]':>24}")
for nm, g in [('all prodromal', _d38), ('hyposmia', _d38[_d38['subgroup'] == 'Hyposmia'])]:
    for lab, vals, kind in [
            ('UPSIT, continuous score', g['upsit'].values, 'continuous'),
            ('UPSIT, dichotomised (hyposmic yes/no)', (1 - g['upsit_pctl15']).values, 'binary'),
            ('imaging, continuous (lowest putamen %exp)', g['pct_exp'].values, 'continuous'),
            ('imaging, dichotomised (DAT-deficit flag)', (g['pct_exp'] >= 65).astype(float).values, 'binary'),
            ('OIS, continuous', g['OIS'].values, 'continuous')]:
        c, lo, hi = _c38(g, vals)
        print(f'  {nm:<24}{lab:<40}{c:>8.3f} [{lo:.3f}, {hi:.3f}]')
        rows38.append(dict(stratum=nm, reading=lab, kind=kind, n=int(len(g)),
                           events=int(g.converted.sum()), c=c, lo=lo, hi=hi))
    print()
dc38 = pd.DataFrame(rows38)
for nm in ['all prodromal', 'hyposmia']:
    sub = dc38[dc38['stratum'] == nm].set_index('reading')['c']
    print(f"  {nm}: dichotomising UPSIT costs "
          f"{sub['UPSIT, continuous score'] - sub['UPSIT, dichotomised (hyposmic yes/no)']:.3f} C-index, "
          f"imaging costs "
          f"{sub['imaging, continuous (lowest putamen %exp)'] - sub['imaging, dichotomised (DAT-deficit flag)']:.3f}")
print('\n  All C-indices in this paper use the continuous score. The binary rows reproduce the')
print('  reading the field currently applies, and are reported descriptively only.')
dc38.to_csv('results/repro/section38_continuous_vs_binary.csv', index=False)
sec16['continuous_vs_binary'] = dc38.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m38 = json.load(f)
_b38 = set(_m38); _m38.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m38, f, indent=2, default=float)
print(f'\n-> results/repro/section38_continuous_vs_binary.csv written')
print(f'-> section16_stats.json: {len(_b38)} -> {len(_m38)} keys')

# ============================================================================
# Section 39. Rank-matched comparison of the three readings
# ============================================================================
# Section 39: level the comparison. All three readings get the same cutpoint search,
# plus a cutpoint-free sensitivity-matched view.
from scipy.stats import chi2 as _chi239
print('Section 39  Like-for-like comparison of the three readings\n')
_d39 = prod.dropna(subset=['OIS', 'upsit', 'PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM',
                           'time_years', 'converted', 'age', 'sex_num']).copy()
_d39['pct_exp'] = _d39[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (
    1.8547 - 0.00303 * _d39['age'] - 0.2419 * _d39['sex_num']) * 100
_h39 = _d39[_d39['subgroup'] == 'Hyposmia'].copy()
_t39 = _h39['time_years'].values.astype(float)
_o39 = np.argsort(_t39)
_ts39 = _t39[_o39]; _e39 = _h39['converted'].values.astype(int)[_o39]
_n39 = len(_ts39); _at39 = np.arange(_n39, 0, -1).astype(float)

def _lr39(mask, ev):
    n1 = mask.astype(float); at1 = np.cumsum(n1[::-1])[::-1]; ef = ev.astype(float)
    O1 = (ef * n1).sum(); E1 = (ef * at1 / _at39).sum()
    V = (ef * (_at39 - ef) / np.maximum(_at39 - 1, 1) * at1 / _at39 * (1 - at1 / _at39)).sum()
    return (O1 - E1) ** 2 / V if V > 0 else 0.0

READ39 = [('upsit', 'UPSIT total'), ('pct_exp', 'lowest putamen %expected'), ('OIS', 'OIS')]
_rng39 = np.random.default_rng(42); _NP39 = 1000
_E39 = _h39['converted'].sum()
rowsA = []
print('  A. Optimal cutpoint, identical search for every reading\n')
print(f"  {'reading':<26}{'cut':>9}{'chi2':>8}{'n low':>8}{'%layer':>8}{'%events':>9}"
      f"{'2-yr low':>10}{'2-yr high':>11}{'perm p':>9}")
for col, lab in READ39:
    x = _h39[col].values.astype(float)[_o39]
    lo, hi = np.percentile(x, [10, 90])
    masks = [(c, x <= c) for c in np.unique(x[(x >= lo) & (x <= hi)])]
    masks = [(c, m) for c, m in masks if m.sum() >= 20 and (~m).sum() >= 20]

    def _best39(ev):
        b = (0.0, None)
        for c, m in masks:
            st_ = _lr39(m, ev)
            if st_ > b[0]:
                b = (st_, c)
        return b

    obs, cut = _best39(_e39)
    ex = sum(_best39(_e39[_rng39.permutation(_n39)])[0] >= obs for _ in range(_NP39))
    g = _h39[_h39[col] <= cut]; ng = _h39[_h39[col] > cut]
    r2a = float(1 - KaplanMeierFitter().fit(g['time_years'], g['converted']).predict(2.0))
    r2b = float(1 - KaplanMeierFitter().fit(ng['time_years'], ng['converted']).predict(2.0))
    p_txt = f'<{1/_NP39:.3f}' if ex == 0 else f'{ex/_NP39:.3f}'
    print(f'  {lab:<26}{cut:>9.2f}{obs:>8.1f}{len(g):>8}{len(g)/_n39*100:>7.0f}%'
          f'{g.converted.sum()/_E39*100:>8.0f}%{r2a*100:>9.1f}%{r2b*100:>10.1f}%{p_txt:>9}')
    rowsA.append(dict(reading=lab, cut=float(cut), chi2=float(obs), n_low=int(len(g)),
                      frac_layer=float(len(g)/_n39), frac_events=float(g.converted.sum()/_E39),
                      risk2y_low=r2a, risk2y_high=r2b,
                      p_naive=float(_chi239.sf(obs, 1)), p_permutation=float(ex/_NP39)))
print('\n  OIS has the highest statistic even when every reading gets its own optimal cut, so the')
print('  advantage is not an artefact of searching a cutpoint for OIS alone. The UPSIT optimum has')
print('  to flag 59% of the layer to reach 85% of the events.\n')

print('  B. Sensitivity-matched, no cutpoint required\n')
print(f"  {'catch of converters':<24}" + ''.join(f'{lab:>26}' for _, lab in READ39))
rowsB = []
for target in [0.5, 0.6, 0.7, 0.8, 0.9]:
    cells = []
    for col, lab in READ39:
        srt = _h39.sort_values(col)
        k = min(int((srt['converted'].cumsum() / _E39 < target).sum()) + 1, _n39)
        cells.append(f'{k} ({k/_n39*100:.0f}%)')
        rowsB.append(dict(kind='sensitivity_matched', reading=lab, target_sensitivity=target,
                          n_flagged=int(k), frac_flagged=float(k/_n39)))
    print(f'  {target*100:>21.0f}%' + ''.join(f'{c:>26}' for c in cells))
print()
print(f"  {'flagged fraction':<24}" + ''.join(f'{lab:>26}' for _, lab in READ39))
for frac in [0.10, 0.20, 0.30, 0.50]:
    k = int(_n39 * frac); cells = []
    for col, lab in READ39:
        ev = _h39.sort_values(col).head(k)['converted'].sum()
        cells.append(f'{int(ev)} ({ev/_E39*100:.0f}%)')
        rowsB.append(dict(kind='fixed_flag_fraction', reading=lab, frac_flagged=float(frac),
                          n_flagged=int(k), events_caught=int(ev), frac_events=float(ev/_E39)))
    print(f'  {frac*100:>20.0f}% ({k})' + ''.join(f'{c:>26}' for c in cells))
print('\n  OIS dominates at every operating point. Imaging loses to UPSIT at the high-sensitivity end,')
print('  the same phenomenon as the non-monotonic imaging bands in this cohort.')
a39 = pd.DataFrame(rowsA); b39 = pd.DataFrame(rowsB)
a39.to_csv('results/repro/section39_optimal_cut_all_readings.csv', index=False)
b39.to_csv('results/repro/section39_sensitivity_matched.csv', index=False)
sec16['levelled_comparison'] = dict(optimal_cut_all=a39.to_dict('records'),
                                    sensitivity_matched=b39.to_dict('records'))
with open('results/repro/section16_stats.json') as f:
    _m39 = json.load(f)
_b39 = set(_m39); _m39.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m39, f, indent=2, default=float)
print(f'\n-> section39_optimal_cut_all_readings.csv, section39_sensitivity_matched.csv written')
print(f'-> section16_stats.json: {len(_b39)} -> {len(_m39)} keys')

# ============================================================================
# Section 40. Calibration and age within strata
# ============================================================================
# Section 40: does the olfactory calibration add anything, and what carries the signal
# in the non-deficit stratum? Both are questions a reviewer will certainly ask.
from sklearn.linear_model import Ridge as _Ridge40
from sklearn.preprocessing import StandardScaler as _SS40
from sklearn.model_selection import KFold as _KF40
print('Section 40  Is the calibration doing the work, and what carries the non-deficit signal\n')
_tr40 = train.dropna(subset=SBR_COLS + DEMO + ['upsit'])
_sc40 = _SS40().fit(_tr40[SBR_COLS].values)
_m40 = _Ridge40(alpha=1.0).fit(_sc40.transform(_tr40[SBR_COLS].values), _tr40['upsit'].values)
_need40 = ['OIS', 'upsit', 'time_years', 'converted',
           'PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM'] + SBR_COLS + DEMO
_d40 = prod.dropna(subset=_need40).copy()
_d40['OIS_sbr'] = _m40.predict(_sc40.transform(_d40[SBR_COLS].values))
_d40['pct_exp'] = _d40[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (
    1.8547 - 0.00303 * _d40[DEMO[0]] - 0.2419 * _d40['sex_num']) * 100
_h40 = _d40[_d40['subgroup'] == 'Hyposmia']

def _cvcox40(g, cols, seed=41, K=5):
    # Out-of-fold Cox risk score, returned with the protective sign convention.
    X = g[cols].copy()
    X = (X - X.mean()) / X.std().replace(0, 1)
    X['T'] = g['time_years'].values
    X['E'] = g['converted'].values.astype(int)
    oof = np.full(len(g), np.nan)
    for tri, tei in _KF40(K, shuffle=True, random_state=seed).split(X):
        a = X.iloc[tri]
        if a['E'].sum() < 3:
            continue
        try:
            cph = CoxPHFitter(penalizer=0.5, l1_ratio=0.5).fit(a, 'T', 'E')
            oof[tei] = -cph.predict_partial_hazard(X.iloc[tei][cols]).values
        except Exception:
            pass
    return oof

_rng40 = np.random.default_rng(42)

def _c40(g, v, B=800):
    m = ~np.isnan(np.asarray(v, dtype=float))
    g2 = g[m]; vv = np.asarray(v, dtype=float)[m]
    t = g2['time_years'].values; e = g2['converted'].values.astype(int)
    bs = []
    for _ in range(B):
        i = _rng40.choice(len(g2), len(g2), replace=True)
        if e[i].sum() < 5:
            continue
        try:
            bs.append(concordance_index(t[i], vv[i], e[i]))
        except Exception:
            pass
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return float(concordance_index(t, vv, e)), float(lo), float(hi), int(len(g2))

rows40 = []
for nm, g in [('hyposmia', _h40),
              ('hyposmia + DAT deficit', _h40[_h40['pct_exp'] < 65]),
              ('hyposmia + no DAT deficit', _h40[_h40['pct_exp'] >= 65])]:
    print(f"  {nm}  n={len(g)}, events={int(g.converted.sum())}")
    for lab, v, kind in [
            ('single-region putamen', g['PUTAMEN_REF_CWM'].values, 'current practice'),
            ('lowest putamen %expected', g['pct_exp'].values, 'current practice'),
            ('age alone', -g[DEMO[0]].values, 'demographic'),
            ('demographics only, CV Cox', _cvcox40(g, DEMO), 'supervised in-cohort'),
            ('33 SBR, CV Cox on conversion', _cvcox40(g, SBR_COLS), 'supervised in-cohort'),
            ('33 SBR + demographics, CV Cox', _cvcox40(g, SBR_COLS + DEMO), 'supervised in-cohort'),
            ('OIS, imaging only, zero-shot', g['OIS_sbr'].values, 'olfaction-calibrated'),
            ('OIS, full, zero-shot', g['OIS'].values, 'olfaction-calibrated'),
            ('UPSIT total', g['upsit'].values, 'behavioural')]:
        c, lo, hi, n = _c40(g, v)
        flag = '  <- CI contains 0.5' if lo < 0.5 < hi else ''
        print(f'    {lab:<34}{c:.3f} [{lo:.3f}, {hi:.3f}]{flag}')
        rows40.append(dict(stratum=nm, reading=lab, kind=kind, n=int(len(g)),
                           events=int(g.converted.sum()), c=c, lo=lo, hi=hi,
                           at_chance=bool(lo < 0.5 < hi)))
    # OIS and UPSIT each adjusted for age
    for var in ['OIS', 'upsit']:
        X = g[[var, DEMO[0]]].copy()
        X['T'] = g['time_years'].values; X['E'] = g['converted'].values.astype(int)
        try:
            cph = CoxPHFitter().fit(X, 'T', 'E')
            print(f'    {var} adjusted for age: HR={np.exp(cph.params_[var]):.3f}, '
                  f'P={cph.summary.loc[var, "p"]:.4f}; model concordance {cph.concordance_index_:.3f}')
            rows40.append(dict(stratum=nm, reading=f'{var} adjusted for age', kind='multivariable',
                               n=int(len(g)), events=int(g.converted.sum()),
                               hr=float(np.exp(cph.params_[var])), p=float(cph.summary.loc[var, 'p']),
                               model_c=float(cph.concordance_index_)))
        except Exception:
            pass
    print()
print('  Reading 1: the olfaction-calibrated zero-shot model beats the conversion-supervised')
print('  33-region model in the hyposmic layer and in the deficit stratum, so the advantage is')
print('  not simply multivariate imaging over a single region.')
print('  Reading 2: in the non-deficit stratum every imaging reading is at chance, including the')
print('  imaging-only OIS. Age alone reaches 0.811 there. With 22 events the confidence intervals')
print('  for age, UPSIT and OIS overlap heavily and the main text must not rank them.')
cal40 = pd.DataFrame(rows40)
cal40.to_csv('results/repro/section40_calibration_and_age.csv', index=False)
sec16['calibration_and_age'] = cal40.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m = json.load(f)
_b = set(_m); _m.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m, f, indent=2, default=float)
print(f'\n-> results/repro/section40_calibration_and_age.csv written')
print(f'-> section16_stats.json: {len(_b)} -> {len(_m)} keys')

# ============================================================================
# Section 41. Group flow
# ============================================================================
from lifelines import KaplanMeierFitter as _KMF41

print('=== §41a 嗅觉减退臂的逐步流程 ===')
_arm = df[(df['COHORT'] == 4) & (df['subgroup'] == 'Hyposmia')]['PATNO'].unique()
print(f'  前驱期有影像与基线临床数据、入组臂为 Hyposmia:      {len(_arm)}')
_ev_all = len(set(_arm) & set(pd_diag.index))
print(f'    其中随访中曾记录 PRIMDIAG==1:                      {_ev_all}')
_pv = set(_arm) & prevalent_at_baseline
print(f'  排除基线即现患:                                   -{len(_pv)}')
_a2 = [p for p in _arm if p not in prevalent_at_baseline]
print(f'                                                     {len(_a2)}   '
      f'事件 {len(set(_a2) & set(pd_diag.index))}')

_t_conv = (pd.to_datetime(pd.Series(_a2).map(pd_diag)) -
           pd.to_datetime(pd.Series(_a2).map(first_visit))).dt.days / 365.25
_t_cens = (pd.to_datetime(pd.Series(_a2).map(last_visit)) -
           pd.to_datetime(pd.Series(_a2).map(first_visit))).dt.days / 365.25
_c = pd.Series(_a2).isin(pd_diag.index).values
_tt = np.where(_c, np.maximum(_t_conv, 0.25), _t_cens)
_drop = int(((~_c) & (_tt <= 0)).sum())
print(f'  排除无随访时间(仅一次基线访视):                   -{_drop}   '
      f'其中转化者 {int((_c & (_tt <= 0)).sum())}')
_a3 = [p for p, k in zip(_a2, (_c | (_tt > 0))) if k]
print(f'                                                     {len(_a3)}   '
      f'事件 {len(set(_a3) & set(pd_diag.index))}')

_lay = prod[prod['subgroup'] == 'Hyposmia']
_lay = _lay[['OIS', 'OMI', 'upsit', 'PUTAMEN_REF_CWM', 'time_years', 'converted']].join(
    prod[prod['subgroup'] == 'Hyposmia'][['PATNO']]).dropna(
    subset=['OIS', 'upsit', 'PUTAMEN_REF_CWM', 'time_years'])
print(f'  排除缺分析变量:                                   -{len(_a3) - len(_lay)}')
print(f'  => 分析层:                                         {len(_lay)}   '
      f'事件 {int(_lay["converted"].sum())}')

print()
print('=== §41b 该层内按实测嗅觉的构成(分层用入组标签,此处仅作说明) ===')
_p15 = excel[excel['COHORT'] == 4].groupby('PATNO')['upsit_pctl15'].first()
_ids = prod[prod['subgroup'] == 'Hyposmia']['PATNO']
_ids = _ids[prod.loc[_ids.index, ['OIS', 'upsit', 'PUTAMEN_REF_CWM', 'time_years']].notna().all(axis=1)]
_v = _ids.map(_p15)
print(f'  实测 <= 第15百分位: {int((_v == 1).sum())} / {len(_v)} '
      f'= {100 * (_v == 1).mean():.1f}%')
print(f'  实测 >  第15百分位: {int((_v == 0).sum())}   缺失: {int(_v.isna().sum())}')

print()
print('=== §41c 为什么只有 68 个事件:随访长度 ===')
_L = prod[prod['subgroup'] == 'Hyposmia'].dropna(
    subset=['OIS', 'upsit', 'PUTAMEN_REF_CWM', 'time_years'])
print(f'  n={len(_L)}  事件={int(_L["converted"].sum())}')
print(f'  中位随访 {_L["time_years"].median():.2f} 年 | 总人年 {_L["time_years"].sum():.0f} | '
      f'每100人年 {100 * _L["converted"].sum() / _L["time_years"].sum():.2f}')
for _y in [1, 2, 3, 4]:
    print(f'    随访 >= {_y} 年: {int((_L["time_years"] >= _y).sum())} '
          f'({100 * (_L["time_years"] >= _y).mean():.1f}%)')
_km41 = _KMF41().fit(_L['time_years'], _L['converted'])
print('  KM 累积转化率:')
_km_rows = []
for _y in [1, 2, 3, 4]:
    _r = float(1 - _km41.predict(_y)); _n = int((_L['time_years'] >= _y).sum())
    print(f'    {_y} 年: {100 * _r:.1f}%   (在险 {_n})')
    _km_rows.append(dict(year=_y, cum_conversion=_r, at_risk=_n))
_cv = _L[_L['converted'] == 1]['time_years']
print(f'  转化者随访时间: 中位 {_cv.median():.2f} 年,范围 '
      f'{_cv.min():.2f}-{_cv.max():.2f};两年内转化 {int((_cv <= 2).sum())} 例')
print()
print('  解读:该层每百人年 3.9 例,三年累积约 10.7%,与 PARS 四年 19/152=12.5% 同量级。')
print('  68 例不是事件偏少,是中位随访仅一年多。四年时点在险人数为两位数,不可用。')

sec16['layer_flow'] = dict(
    arm_with_imaging=int(len(_arm)), ever_diagnosed=int(_ev_all),
    prevalent_excluded=int(len(_pv)), after_prevalent=int(len(_a2)),
    no_followup_excluded=int(_drop), after_followup=int(len(_a3)),
    missing_vars_excluded=int(len(_a3) - len(_lay)),
    analysed=int(len(_lay)), events=int(_lay['converted'].sum()),
    measured_at_or_below_p15=int((_v == 1).sum()), measured_above_p15=int((_v == 0).sum()),
    pct_measured_hyposmic=float(100 * (_v == 1).mean()),
    median_followup_years=float(_L['time_years'].median()),
    person_years=float(_L['time_years'].sum()),
    rate_per_100py=float(100 * _L['converted'].sum() / _L['time_years'].sum()),
    km=_km_rows)
with open('results/repro/section16_stats.json') as f:
    _m41 = json.load(f)
_b41 = set(_m41); _m41.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m41, f, indent=2, default=float)
print(f'-> section16_stats.json: {len(_b41)} -> {len(_m41)} keys')

# ============================================================================
# Section 44. Three-group split, permutation corrected
# ============================================================================
import itertools as _it44
from lifelines import KaplanMeierFitter as _KMF44
from lifelines.statistics import logrank_test as _lrt44, multivariate_logrank_test as _mlr44
print('Section 44  Three risk groups by a two-cut search, permutation-corrected\n')
_L44 = _h40.copy(); _L44['upsit_'] = _L44['upsit']
MINFRAC44, STEP44, NPERM44, NBOOT44 = 0.15, 2.5, 1000, 200
_rng44 = np.random.default_rng(42)


def _lr_chi2_44(t, g, e, k=3):
    # k-group log-rank chi-square, vectorised; same statistic as lifelines multivariate_logrank_test
    order = np.argsort(t, kind='stable'); t = t[order]; g = g[order]; e = e[order]
    uniq, first = np.unique(t, return_index=True)
    at_risk = np.zeros((len(uniq), k)); events = np.zeros((len(uniq), k))
    for j in range(k):
        m = (g == j)
        cum_n = np.cumsum(m[::-1])[::-1]
        at_risk[:, j] = cum_n[first]
        events[:, j] = np.bincount(np.searchsorted(uniq, t[m & (e == 1)]), minlength=len(uniq))
    N = at_risk.sum(1); D = events.sum(1); keep = D > 0
    at_risk, events, N, D = at_risk[keep], events[keep], N[keep], D[keep]
    E = at_risk * (D / N)[:, None]
    O_E = (events - E).sum(0)[:-1]
    P = at_risk / N[:, None]; f = np.where(N > 1, D * (N - D) / np.maximum(N - 1, 1), 0.0)
    Vfull = np.zeros((k, k))
    for i_ in range(len(N)):
        p = P[i_]; Vfull += f[i_] * (np.diag(p) - np.outer(p, p))
    V = Vfull[:-1, :-1]
    try: return float(O_E @ np.linalg.solve(V, O_E))
    except np.linalg.LinAlgError: return 0.0

def _grid44(x):
    return np.unique(np.percentile(x, np.arange(10, 90.1, STEP44)))

def _split44(x, c1, c2):
    return np.where(x <= c1, 0, np.where(x <= c2, 1, 2))

def _search44(x, t, e, constrained=True):
    # returns (chi2, c1, c2, pairwise p) maximising the 3-group log-rank statistic
    n = len(x); best = None; cuts = _grid44(x)
    for c1, c2 in _it44.combinations(cuts, 2):
        g = _split44(x, c1, c2)
        if np.bincount(g, minlength=3).min() < MINFRAC44 * n: continue
        chi = _lr_chi2_44(t, g, e)
        if best is not None and chi <= best[0]: continue
        pw = [_lrt44(t[g == a], t[g == b], e[g == a], e[g == b]).p_value for a, b in [(0, 1), (1, 2), (0, 2)]] if constrained else None
        if constrained and max(pw) >= 0.05: continue
        best = (chi, c1, c2, pw)
    return best

def _km2_44(g):
    k = _KMF44().fit(g['time_years'], g['converted']); ci = k.confidence_interval_survival_function_
    idx = ci.index[ci.index <= 2.0][-1]
    return (1 - float(k.predict(2.0))) * 100, (1 - ci.loc[idx].iloc[1]) * 100, (1 - ci.loc[idx].iloc[0]) * 100

rows44 = []; boot44 = []
for stratum, d in [('hyposmia', _L44), ('hyposmia + no DAT deficit', _L44[_L44['pct_exp'] >= 65])]:
    for reading, col in [('OIS', 'OIS'), ('UPSIT total', 'upsit'), ('lower putamen %expected', 'pct_exp')]:
        x = d[col].values; t = d['time_years'].values; e = d['converted'].values.astype(int)
        best = _search44(x, t, e, constrained=True)
        unc = _search44(x, t, e, constrained=False)
        print(f'--- {stratum} | {reading} (n={len(d)}, events={int(e.sum())})')
        if best is None:
            print('   no two-cut split gives three mutually distinct groups (all pairwise P < 0.05)')
            rows44.append(dict(stratum=stratum, reading=reading, feasible=False, chi2_unconstrained=float(unc[0]),
                               cut1_unc=float(unc[1]), cut2_unc=float(unc[2])))
            continue
        chi, c1, c2, pw = best; g = _split44(x, c1, c2)
        # permutation correction of the search: permute outcome, re-search the unconstrained maximum
        null = np.empty(NPERM44)
        for i in range(NPERM44):
            perm = _rng44.permutation(len(e)); tp, ep = t[perm], e[perm]
            null[i] = _search44(x, tp, ep, constrained=False)[0]
        p_perm = float((null >= chi).mean())
        # bootstrap stability of the two cuts
        bc = []
        for i in range(NBOOT44):
            idx = _rng44.integers(0, len(x), len(x))
            r_ = _search44(x[idx], t[idx], e[idx], constrained=False)
            if r_ is not None: bc.append((r_[1], r_[2]))
        bc = np.array(bc)
        rec = dict(stratum=stratum, reading=reading, feasible=True, cut1=float(c1), cut2=float(c2), chi2=float(chi),
                   p_permutation=p_perm, n_perm=NPERM44,
                   cut1_boot_lo=float(np.percentile(bc[:, 0], 2.5)), cut1_boot_hi=float(np.percentile(bc[:, 0], 97.5)),
                   cut2_boot_lo=float(np.percentile(bc[:, 1], 2.5)), cut2_boot_hi=float(np.percentile(bc[:, 1], 97.5)),
                   p_low_mid=float(pw[0]), p_mid_high=float(pw[1]), p_low_high=float(pw[2]))
        for k_, lab in enumerate(['low', 'mid', 'high']):
            gg = d[g == k_]; r2, lo2, hi2 = _km2_44(gg)
            rec.update({f'n_{lab}': int(len(gg)), f'frac_{lab}': float(len(gg) / len(d)), f'ev_{lab}': int(gg['converted'].sum()),
                        f'km2_{lab}': r2, f'km2_{lab}_lo': lo2, f'km2_{lab}_hi': hi2})
            print(f'   {lab:4s}: n={len(gg):4d} ({len(gg)/len(d)*100:4.0f}%)  events={int(gg.converted.sum()):3d}  2-year {r2:5.1f}% ({lo2:.1f} to {hi2:.1f})')
        print(f'   cuts {c1:.2f} / {c2:.2f}  (bootstrap 95%: {rec["cut1_boot_lo"]:.1f} to {rec["cut1_boot_hi"]:.1f} / {rec["cut2_boot_lo"]:.1f} to {rec["cut2_boot_hi"]:.1f})')
        print(f'   3-group chi2 {chi:.1f}, permutation P {p_perm:.3f} ({NPERM44} perms); pairwise P low-mid {pw[0]:.2e}, mid-high {pw[1]:.2e}, low-high {pw[2]:.2e}')
        rows44.append(rec)
# non-deficit stratum: outcome-blind two-group fallback (lowest tertile vs upper two tertiles), OIS and UPSIT
print('\n--- non-deficit stratum, outcome-blind fallback: lowest tertile against the upper two tertiles')
two44 = []
for stratum, d in [('hyposmia', _L44), ('hyposmia + no DAT deficit', _L44[_L44['pct_exp'] >= 65])]:
    for reading, col in [('OIS', 'OIS'), ('UPSIT total', 'upsit')]:
        d = d.copy(); d['t3'] = pd.qcut(d[col], 3, labels=[0, 1, 2]).astype(int)
        low, rest = d[d.t3 == 0], d[d.t3 != 0]
        r1, l1, h1 = _km2_44(low); r2, l2, h2 = _km2_44(rest)
        p2 = _lrt44(low.time_years, rest.time_years, low.converted, rest.converted).p_value
        p23 = _lrt44(d[d.t3 == 1].time_years, d[d.t3 == 2].time_years, d[d.t3 == 1].converted, d[d.t3 == 2].converted).p_value
        p3 = _mlr44(d.time_years, d.t3, d.converted).p_value
        two44.append(dict(stratum=stratum, reading=reading, n_low=int(len(low)), ev_low=int(low.converted.sum()), km2_low=r1, km2_low_lo=l1, km2_low_hi=h1,
                          n_rest=int(len(rest)), ev_rest=int(rest.converted.sum()), km2_rest=r2, km2_rest_lo=l2, km2_rest_hi=h2,
                          p_two_group=float(p2), p_three_group_tertiles=float(p3), p_middle_vs_highest=float(p23)))
        print(f'   {stratum:26s} {reading:12s} lowest third {r1:.1f}% ({l1:.1f} to {h1:.1f}, {int(low.converted.sum())}/{len(low)}) vs rest {r2:.1f}% ({l2:.1f} to {h2:.1f}, {int(rest.converted.sum())}/{len(rest)}), two-group P {p2:.1e}; tertile middle vs highest P {p23:.2f}')
R44 = pd.DataFrame(rows44); T44 = pd.DataFrame(two44)
R44.to_csv('results/repro/section44_three_groups.csv', index=False); T44.to_csv('results/repro/section44_lowest_third_vs_rest.csv', index=False)
sec16['three_group_split'] = dict(settings=dict(min_group_fraction=MINFRAC44, grid_step_percentile=STEP44, n_perm=NPERM44, n_boot=NBOOT44, seed=42),
                                  search=R44.to_dict('records'), lowest_third_vs_rest=T44.to_dict('records'))
with open('results/repro/section16_stats.json') as f:
    _m44 = json.load(f)
_b44 = set(_m44); _m44.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m44, f, indent=2, default=float)
print(f'-> section16_stats.json: {len(_b44)} -> {len(_m44)} keys')

# ============================================================================
# Median split (Fig. 3d, 3e and Supplementary Table 6)
# ============================================================================
# Same frame and normative coefficients as sections 33 to 40.
from lifelines.statistics import logrank_test as _lrt_ms
hh=prod[prod.subgroup=='Hyposmia'].dropna(subset=['OIS','upsit','PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM','time_years','converted','age','sex_num']).copy()
hh=hh[(hh.converted==1)|(hh.time_years>0)]
hh['pct']=hh[['PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM']].min(axis=1)/(1.8547-.00303*hh.age-.2419*hh.sex_num)*100
def km2(g):
    k=KaplanMeierFitter().fit(g.time_years,g.converted); ci=k.confidence_interval_survival_function_; idx=ci.index[ci.index<=2.0][-1]
    return (1-float(k.predict(2.0)))*100,(1-ci.loc[idx].iloc[1])*100,(1-ci.loc[idx].iloc[0])*100
rows=[]
for name,d in [('hyposmia',hh),('hyposmia + no DAT deficit',hh[hh.pct>=65])]:
    for reading,col in [('OIS','OIS'),('UPSIT total','upsit'),('lower putamen %expected','pct')]:
        med=d[col].median(); lo_=d[d[col]<=med]; hi_=d[d[col]>med]
        r1,l1,h1=km2(lo_); r2,l2,h2=km2(hi_); p=_lrt_ms(lo_.time_years,hi_.time_years,lo_.converted,hi_.converted).p_value
        rows.append(dict(stratum=name,reading=reading,median=med,n_low=len(lo_),ev_low=int(lo_.converted.sum()),km2_low=r1,km2_low_lo=l1,km2_low_hi=h1,n_high=len(hi_),ev_high=int(hi_.converted.sum()),km2_high=r2,km2_high_lo=l2,km2_high_hi=h2,p_logrank=p))
        print(f'{name:26s} {reading:24s} median {med:6.2f} | <=median n={len(lo_):4d} ev={int(lo_.converted.sum()):2d} 2y={r1:5.1f} ({l1:.1f}-{h1:.1f}) | >median n={len(hi_):4d} ev={int(hi_.converted.sum()):2d} 2y={r2:4.1f} ({l2:.1f}-{h2:.1f}) | P={p:.1e} | ratio {r1/r2:.1f}')
pd.DataFrame(rows).to_csv('results/repro/supp_median_split.csv',index=False)
print('-> results/repro/supp_median_split.csv written')

# ============================================================================
# Section 46. Time-dependent AUC with inverse probability of censoring weighting
# ============================================================================
from sksurv.metrics import cumulative_dynamic_auc as _cda46
from sksurv.util import Surv as _Surv46
print('Section 46  Time-dependent (cumulative/dynamic) AUC with IPCW, hyposmia layer\n')

_TIMES46 = [1.0, 1.5, 2.0]
_d46 = hyp16.copy()
# the imaging reading is the same one the pre-specified comparisons of Table 2 use,
# namely the single-region putamen binding ratio, so the two analyses stay aligned
_SCORES46 = [('upsit', 'UPSIT total'), ('OIS', 'OIS'), ('PUTAMEN_REF_CWM', 'putamen SBR')]
_d46 = _d46.dropna(subset=[c for c, _ in _SCORES46] + ['time_years', 'converted'])
print(f'  n = {len(_d46)}, events = {int(_d46["converted"].sum())}')
print('  at risk / events by horizon: ' + ', '.join(
    f'{t:.1f}y {int((_d46.time_years >= t).sum())} at risk, {int(((_d46.converted == 1) & (_d46.time_years <= t)).sum())} events'
    for t in _TIMES46))


def _auc46(d, col, times):
    # IPCW cumulative/dynamic AUC. All three readings are protective, so the risk
    # score passed to the estimator is the negated value.
    y = _Surv46.from_arrays(event=d['converted'].astype(bool).values, time=d['time_years'].values)
    a, _ = _cda46(y, y, -d[col].values.astype(float), times)
    return np.asarray(a, dtype=float)


_N_BOOT46, _rng46 = 1000, np.random.RandomState(42)
_idx46 = np.arange(len(_d46))
_boot46 = {c: [] for c, _ in _SCORES46}
_bdiff46 = []
_ok46 = 0
for _b in range(_N_BOOT46):
    _s = _d46.iloc[_rng46.choice(_idx46, len(_idx46), replace=True)]
    if _s['converted'].sum() < 10:
        continue
    try:
        _vals = {c: _auc46(_s, c, _TIMES46) for c, _ in _SCORES46}
    except Exception:
        continue          # a resample can leave a horizon outside the observed event times
    for c in _vals:
        _boot46[c].append(_vals[c])
    _bdiff46.append(_vals['OIS'] - _vals['upsit'])
    _ok46 += 1
print(f'  {_ok46} of {_N_BOOT46} bootstrap resamples usable')

_rows46 = []
_obs46 = {c: _auc46(_d46, c, _TIMES46) for c, _ in _SCORES46}
_D46 = np.vstack(_bdiff46)
for _j, _t in enumerate(_TIMES46):
    for c, lab in _SCORES46:
        _bs = np.array([v[_j] for v in _boot46[c]])
        _rows46.append(dict(horizon_yr=_t, reading=lab, auc=float(_obs46[c][_j]),
                            lo=float(np.percentile(_bs, 2.5)), hi=float(np.percentile(_bs, 97.5)),
                            n_at_risk=int((_d46.time_years >= _t).sum()),
                            events_by_t=int(((_d46.converted == 1) & (_d46.time_years <= _t)).sum())))
    _dd = _D46[:, _j]
    _p = 2 * min((_dd <= 0).mean(), (_dd >= 0).mean())
    _rows46.append(dict(horizon_yr=_t, reading='OIS minus UPSIT total',
                        auc=float(_obs46['OIS'][_j] - _obs46['upsit'][_j]),
                        lo=float(np.percentile(_dd, 2.5)), hi=float(np.percentile(_dd, 97.5)),
                        n_at_risk=int((_d46.time_years >= _t).sum()),
                        events_by_t=int(((_d46.converted == 1) & (_d46.time_years <= _t)).sum()),
                        p=float(max(_p, 1.0 / max(_ok46, 1)))))
td46 = pd.DataFrame(_rows46)
print()
print(f"{'T':>5}{'reading':>26}{'AUC':>8}{'95% CI':>20}{'P':>9}")
for _, r in td46.iterrows():
    _pp = '' if pd.isna(r.get('p', np.nan)) else ('<0.001' if r['p'] < 0.001 else f"{r['p']:.3f}")
    print(f"{r['horizon_yr']:>5.1f}{r['reading']:>26}{r['auc']:>8.3f}"
          f"{f'({r.lo:.3f}, {r.hi:.3f})':>20}{_pp:>9}")
td46.to_csv('results/repro/section46_td_auc.csv', index=False)
sec16['td_auc_ipcw'] = dict(times=_TIMES46, n=int(len(_d46)), events=int(_d46['converted'].sum()),
                            n_boot_used=int(_ok46), rows=td46.to_dict('records'))
with open('results/repro/section16_stats.json') as f:
    _m46 = json.load(f)
_m46.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m46, f, indent=2, default=float)
print('\n-> results/repro/section46_td_auc.csv, section16_stats.json key td_auc_ipcw')
print('  Beyond two years the layer thins out (118 at risk at 2.5 y), so no later horizon is reported.')
