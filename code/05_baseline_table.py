"""05_baseline_table.py  Baseline characteristics (Table 1)

Outputs: results/repro/section47_baseline.csv and the ledger key baseline_table.
Notebook provenance: section 47. Code is the notebook cell text, unchanged except for
the file paths, which are resolved through common.py.
"""
from common import *
st = load_state('02_cohort'); globals().update(st)
sec16 = {}
pct_expected_putamen = make_pct_expected_putamen(NORM_BETA, age_col)

# ============================================================================
# Section 47. Baseline table
# ============================================================================
print('Section 47  Baseline characteristics\n')

# df 只带了建模需要的几列,量表与病程另从 excel 的基线访视并进来
_EXTRA47 = ['race', 'upsit_pctl15', 'updrs1_score', 'updrs2_score', 'hy', 'ess', 'gds', 'stai',
            'scopa', 'scopa_gi', 'ageonset', 'duration_yrs', 'LEDD']
_bl47 = excel[excel['EVENT_ID'] == 'BL'].drop_duplicates('PATNO')[['PATNO'] + _EXTRA47]

_tr47 = df[df['COHORT'].isin([1, 2])].dropna(subset=FEATURES + ['upsit']).copy()
# 与生存分析同一帧,prod 尚未剔除缺分析变量的 11 人
_pr47 = prod.dropna(subset=['OIS', 'upsit', 'time_years', 'converted']).copy()
_pr47['grp'] = np.where(_pr47['subgroup'] == 'Hyposmia', 'Hyposmia layer', 'Non-hyposmia-enriched layer')
for _n in ('_tr47', '_pr47'):
    _f = eval(_n).merge(_bl47, on='PATNO', how='left')
    # pct_expected_putamen 返回分数,这里换算成百分比以便与 65% 的判定阈值一致
    _f['pct_exp'] = 100 * pct_expected_putamen(_f)
    _f['put_low'] = _f[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1)
    globals()[_n] = _f

_G47 = [("Parkinson's disease (training)", _tr47[_tr47['COHORT'] == 1], 'train'),
        ('Healthy controls (training)', _tr47[_tr47['COHORT'] == 2], 'train'),
        ('Hyposmia layer', _pr47[_pr47['grp'] == 'Hyposmia layer'], 'prodromal'),
        ('Non-hyposmia-enriched layer', _pr47[_pr47['grp'] != 'Hyposmia layer'], 'prodromal')]


def _cont47(g, col):
    x = pd.to_numeric(g[col], errors='coerce').dropna()
    return (f'{x.mean():.1f} ({x.std():.1f})', len(x)) if len(x) else ('NA', 0)


def _bin47(g, mask):
    m = mask.dropna()
    return (f'{int(m.sum())} ({100 * m.mean():.0f}%)', len(m)) if len(m) else ('NA', 0)


# (section, label, kind, accessor)
_SPEC47 = [
    ('Demographics', 'Age, years', 'c', 'age'),
    ('Demographics', 'Female', 'b', lambda g: (1 - pd.to_numeric(g['sex_num'], errors='coerce')).astype('boolean')),
    ('Demographics', 'White race', 'b', lambda g: (pd.to_numeric(g['race'], errors='coerce') == 1).astype('boolean')),
    ('Demographics', 'Education, years', 'c', 'educyrs'),
    ('Olfaction', 'UPSIT total', 'c', 'upsit'),
    ('Olfaction', 'UPSIT at or below 15th percentile', 'b',
     lambda g: pd.to_numeric(g['upsit_pctl15'], errors='coerce').astype('boolean')),
    ('Imaging', 'Putamen binding ratio, bilateral', 'c', 'PUTAMEN_REF_CWM'),
    ('Imaging', 'Lower putamen binding ratio', 'c', 'put_low'),
    ('Imaging', 'Lower putamen, % of expected', 'c', 'pct_exp'),
    ('Imaging', 'DAT deficit, below 65% of expected', 'b',
     lambda g: (pd.to_numeric(g['pct_exp'], errors='coerce') < 65).astype('boolean')),
    ('Components', 'OIS', 'c', 'OIS'),
    ('Components', 'OMI', 'c', 'OMI'),
    ('Clinical', 'MoCA', 'c', 'moca'),
    ('Clinical', 'MDS-UPDRS part I', 'c', 'updrs1_score'),
    ('Clinical', 'MDS-UPDRS part II', 'c', 'updrs2_score'),
    ('Clinical', 'MDS-UPDRS part III', 'c', 'updrs3_score'),
    ('Clinical', 'Hoehn and Yahr', 'c', 'hy'),
    ('Clinical', 'SCOPA-AUT total', 'c', 'scopa'),
    ('Clinical', 'SCOPA-AUT gastrointestinal', 'c', 'scopa_gi'),
    ('Clinical', 'Epworth sleepiness scale', 'c', 'ess'),
    ('Clinical', 'Geriatric depression scale', 'c', 'gds'),
    ('Clinical', 'State-trait anxiety inventory', 'c', 'stai'),
    ('Parkinson-specific', 'Age at onset, years', 'c', 'ageonset'),
    ('Parkinson-specific', 'Disease duration, years', 'c', 'duration_yrs'),
    ('Parkinson-specific', 'Levodopa equivalent daily dose, mg', 'c', 'LEDD'),
]

_rows47 = []
for _sec, _lab, _kind, _acc in _SPEC47:
    _row = dict(section=_sec, variable=_lab)
    for _gl, _g, _type in _G47:
        if _sec == 'Parkinson-specific' and _type != 'train':
            _row[_gl], _row[_gl + ' n'] = '', 0
            continue
        _v, _n = (_cont47(_g, _acc) if _kind == 'c' else _bin47(_g, _acc(_g)))
        _row[_gl], _row[_gl + ' n'] = _v, _n
    _rows47.append(_row)

# 前驱期特有:入组通道构成与随访
for _lab, _fn in [('Recruited through olfactory screening', lambda g: (g['subgroup'] == 'Hyposmia')),
                  ('Recruited through REM sleep behaviour disorder',
                   lambda g: g['subgroup'].astype(str).str.contains('RBD') & (g['subgroup'] != 'Hyposmia')),
                  ('Recruited through pathogenic variant carriage',
                   lambda g: ~g['subgroup'].astype(str).str.contains('RBD') & (g['subgroup'] != 'Hyposmia'))]:
    _row = dict(section='Recruitment channel', variable=_lab)
    for _gl, _g, _type in _G47:
        _row[_gl], _row[_gl + ' n'] = (_bin47(_g, _fn(_g).astype('boolean')) if _type == 'prodromal' else ('', 0))
    _rows47.append(_row)
for _lab, _fn, _fmt in [('Follow-up, median years', lambda g: g['time_years'].median(), '{:.2f}'),
                        ('Person-years', lambda g: g['time_years'].sum(), '{:,.0f}'),
                        ('Conversions to Parkinson disease', lambda g: g['converted'].sum(), '{:.0f}'),
                        ('Conversions per 100 person-years',
                         lambda g: 100 * g['converted'].sum() / g['time_years'].sum(), '{:.2f}')]:
    _row = dict(section='Follow-up', variable=_lab)
    for _gl, _g, _type in _G47:
        _row[_gl], _row[_gl + ' n'] = ((_fmt.format(_fn(_g)), len(_g)) if _type == 'prodromal' else ('', 0))
    _rows47.append(_row)

base47 = pd.DataFrame(_rows47)
_cols47 = [g[0] for g in _G47]
print('  group sizes: ' + ', '.join(f'{g[0]} {len(g[1]):,}' for g in _G47))
print()
_w = 42
print(f"{'variable':{_w}}" + ''.join(f'{c[:26]:>28}' for c in _cols47))
_prev = None
for _, r in base47.iterrows():
    if r['section'] != _prev:
        print(f"  [{r['section']}]"); _prev = r['section']
    print(f"{r['variable']:{_w}}" + ''.join(f"{r[c]:>28}" for c in _cols47))
base47.to_csv('results/repro/section47_baseline.csv', index=False)
sec16['baseline_table'] = dict(groups=_cols47, group_n={g[0]: int(len(g[1])) for g in _G47},
                               rows=base47.to_dict('records'))
with open('results/repro/section16_stats.json') as f:
    _m47 = json.load(f)
_m47.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m47, f, indent=2, default=float)
print('\n-> results/repro/section47_baseline.csv, ledger key baseline_table')
