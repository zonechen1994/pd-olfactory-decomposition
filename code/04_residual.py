"""04_residual.py  What the residual carries

Genotype contrast (GBA1 against LRRK2), the biological staging axes, CSF tau and
neurofilament, the residual by biological stage, the training-cohort alternatives, prodromal
seed amplification, the genotype premise, cognitive endpoints, age in the prodromal residual,
the partial correlations placing the two tests on one metric, and the tau tertile means of Fig. 1d.
Outputs: results/repro/section22_csf_threeway.csv, section42_partial_r.csv, section45_tau_groups.csv
and ledger keys.
Notebook provenance: sections 20, 21, 22.1, 13 (stage table only), 26, 27, 28, 30, 32, 42 and 45. Code is the notebook cell text, unchanged except for
the file paths, which are resolved through common.py.
"""
from common import *
st = load_state('02_cohort'); globals().update(st)
sec16 = {}

# ============================================================================
# Section 20. Genotype contrast
# ============================================================================
# Section 20 (biological interpretability): the decomposition is tested against genotype.
# GBA encodes a lysosomal enzyme and LRRK2 a kinase, and the two genotypes are known to differ
# in olfactory involvement relative to their dopaminergic state. Both groups enter PPMI by
# genotype rather than by smell, so the comparison is free of olfactory selection.
print('§20 Genotype test of the decomposition (GBA versus LRRK2)\n')
gen20 = {}
rows20 = []
for coh, name in [(4, 'Prodromal'), (1, 'PD')]:
    d = df[(df['COHORT'] == coh) & df['subgroup'].isin(['GBA', 'LRRK2'])].dropna(
        subset=['upsit', 'OIS', 'OMI', 'PUTAMEN_REF_CWM', 'age', 'sex_num', 'educyrs']).copy()
    d['is_GBA'] = (d['subgroup'] == 'GBA').astype(int)
    n_g = int(d['is_GBA'].sum()); n_l = len(d) - n_g
    print(f'  {name}: GBA n={n_g}, LRRK2 n={n_l}   (adjusted for age, sex, education)')
    print(f"    {'measure':<16}{'GBA':>9}{'LRRK2':>9}{'beta':>10}{'95% CI':>20}{'p':>10}")
    for col, lab in [('upsit', 'UPSIT total'), ('PUTAMEN_REF_CWM', 'Putamen SBR'),
                     ('OIS', 'OIS'), ('OMI', 'OMI')]:
        m = sm.OLS(d[col], sm.add_constant(d[['is_GBA', 'age', 'sex_num', 'educyrs']])).fit()
        b, p = m.params['is_GBA'], m.pvalues['is_GBA']
        lo, hi = m.conf_int().loc['is_GBA']
        rows20.append(dict(cohort=name, measure=lab, n_gba=n_g, n_lrrk2=n_l,
                           mean_gba=d[d.is_GBA == 1][col].mean(), mean_lrrk2=d[d.is_GBA == 0][col].mean(),
                           beta=b, lo=lo, hi=hi, p=p))
        star = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        print(f'    {lab:<16}{d[d.is_GBA==1][col].mean():>9.2f}{d[d.is_GBA==0][col].mean():>9.2f}'
              f'{b:>+10.3f}   [{lo:+.3f}, {hi:+.3f}]{p:>10.4f}  {star}')
    print()
gen20 = pd.DataFrame(rows20)
sec16['genotype'] = gen20.to_dict('records')
print('  Reading: in prodromals the two genotypes have the same total score but opposite components.')
print('  In PD they have the same dopaminergic state, and the entire olfactory difference falls on OMI')
print('  while OIS correctly shows none. The decomposition assigns the GBA deficit to the right half.\n')

# (b) direct molecular test: peripheral glucocerebrosidase activity
bio20 = pd.read_csv(BIO_FILE, low_memory=False)
gc = bio20[bio20['TESTNAME'] == 'GCase activity'].copy()
gc['val'] = pd.to_numeric(gc['TESTVALUE'], errors='coerce')
gc = gc.dropna(subset=['val']).groupby('PATNO')['val'].mean().rename('gcase').reset_index()
dg = df.merge(gc, on='PATNO', how='inner')
print(f'§20b GCase activity available for {len(dg)} subjects with a decomposition')
print(f"    {'cohort':<12}{'n':>6}{'score':>8}{'beta per +1SD':>15}{'p':>10}")
gc20 = []
for coh, name in [(4, 'Prodromal'), (1, 'PD')]:
    x = dg[dg['COHORT'] == coh].dropna(subset=['gcase', 'OMI', 'OIS', 'upsit', 'age', 'sex_num', 'educyrs']).copy()
    if len(x) < 50:
        continue
    X = x[['gcase', 'age', 'sex_num', 'educyrs']].astype(float).copy()
    X['gcase'] = (X['gcase'] - X['gcase'].mean()) / X['gcase'].std()
    for col in ['upsit', 'OIS', 'OMI']:
        m = sm.OLS(x[col].astype(float), sm.add_constant(X)).fit()
        gc20.append(dict(cohort=name, score=col, n=len(x), beta=m.params['gcase'], p=m.pvalues['gcase']))
        print(f'    {name:<12}{len(x):>6}{col:>8}{m.params["gcase"]:>+15.3f}{m.pvalues["gcase"]:>10.4f}')
sec16['gcase'] = gc20
print('    Direction is consistent with the lysosomal hypothesis in both samples, and borderline in both.')

# Figure: all three olfactory quantities are in UPSIT points, so the split is directly visible
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3), sharey=True)
for ax, name in zip(axes, ['Prodromal', 'PD']):
    g = gen20[(gen20['cohort'] == name) & (gen20['measure'] != 'Putamen SBR')]
    order = ['UPSIT total', 'OIS', 'OMI']
    g = g.set_index('measure').loc[order].reset_index()
    ys = np.arange(len(g))[::-1]
    cols = [COL_3W['upsit'], COL_3W['OIS'], COL_3W['OMI']]
    ax.errorbar(g['beta'], ys, xerr=[g['beta'] - g['lo'], g['hi'] - g['beta']],
                fmt='o', ms=8, lw=1.8, capsize=5, color='black', zorder=1, ls='none')
    ax.scatter(g['beta'], ys, s=90, c=cols, zorder=2, edgecolors='black')
    ax.axvline(0, color='grey', lw=1, ls='--')
    for y, r in zip(ys, g.itertuples()):
        star = '***' if r.p < 0.001 else ('**' if r.p < 0.01 else ('*' if r.p < 0.05 else 'ns'))
        ax.text(r.hi + 0.35, y, f'p = {r.p:.3g} {star}', va='center', fontsize=9, family='monospace')
    ax.set_yticks(ys); ax.set_yticklabels(order, fontsize=10.5)
    pt = gen20[(gen20['cohort'] == name) & (gen20['measure'] == 'Putamen SBR')].iloc[0]
    ax.set_title(f'{name}  (GBA n={g["n_gba"].iloc[0]}, LRRK2 n={g["n_lrrk2"].iloc[0]})\n'
                 f'putamen SBR difference p = {pt.p:.3g}', fontsize=10.5, loc='left', fontweight='bold')
    ax.set_xlabel('GBA minus LRRK2, adjusted (UPSIT points)', fontsize=10)
axes[0].set_xlim(-3.5, 5.5); axes[1].set_xlim(-9, 4)
fig.suptitle('Genotype assigns the olfactory difference to the expected component',
             fontsize=11.5, fontweight='bold', x=0.01, ha='left')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Fig3_v5_genotype_dissociation.png', dpi=200)
plt.savefig(f'{FIG_DIR}/Fig3_v5_genotype_dissociation.pdf')
plt.show()

merge_ledger(sec16)
print('\n-> Fig3_v5_genotype_dissociation saved; section16_stats.json updated with §20')

# ============================================================================
# Section 21. Biological staging axes
# ============================================================================
# Section 21.1: the two components are tested against the two axes of the MDS biological
# staging system. Stage_S is assigned from CSF alpha-synuclein seed amplification and Stage_D
# from DaTscan. A valid decomposition should map each component onto one axis and neither onto
# the other. Note that OIS versus the D axis is partly circular by construction and serves as a
# positive control; the other cells are not.
from sklearn.metrics import roc_curve
from matplotlib.gridspec import GridSpec

bl21 = excel[excel['EVENT_ID'] == 'BL'][['PATNO', 'COHORT', 'Stage_S', 'Stage_D', 'NHY',
                                         'NSD_STAGE', 'updrs3_score']].drop_duplicates('PATNO')
d21 = bl21.merge(df[['PATNO', 'OIS', 'OMI', 'upsit']], on='PATNO')
_map21 = {'Not NSD': 0, '1a': 1, '1b': 1.5, '2a': 2, '2b': 2.5}
d21['nsd'] = d21['NSD_STAGE'].astype(str).map(lambda v: _map21.get(v, pd.to_numeric(v, errors='coerce')))
SCORES21 = ['upsit', 'OIS', 'OMI']
LAB21 = {'upsit': 'UPSIT (raw)', 'OIS': 'OIS', 'OMI': 'OMI'}
cells21 = []


def _boot_p21(fn, n, seed=42, B=300):
    rng = np.random.RandomState(seed); out = []
    for _ in range(B):
        v = fn(rng.choice(np.arange(n), n, replace=True))
        if v is not None:
            out.append(v)
    return max((np.array(out) <= 0.5).mean(), 1 / B) if out else np.nan


def _axis21(col, name):
    v = d21.dropna(subset=[col] + SCORES21); v = v[v[col].isin([0, 1])]
    y = v[col].astype(int)
    print(f'  {name}: n={len(v)}, positive={int(y.sum())}')
    for s in SCORES21:
        a = roc_auc_score(y, -v[s])
        p = _boot_p21(lambda b: roc_auc_score(y.iloc[b], -v[s].iloc[b]) if y.iloc[b].nunique() > 1 else None, len(v))
        cells21.append(dict(score=s, endpoint=name, stat=f'{a:.3f}', p=p, eff=abs(a - 0.5) * 2))
        print(f'    {LAB21[s]:12s} AUC = {a:.3f}')
    return v, y


print('§21.1 The two axes of the MDS biological staging system\n')
vS21, yS21 = _axis21('Stage_S', 'S axis')
vD21, yD21 = _axis21('Stage_D', 'D axis')

print('\n  Ordinal and clinical endpoints (PD cohort, Spearman or standardised beta):')
for col, name in [('nsd', 'NSD stage'), ('NHY', 'Hoehn-Yahr')]:
    v = d21[d21['COHORT'] == 1].dropna(subset=[col] + SCORES21)
    for s in SCORES21:
        r, p = stats.spearmanr(v[s], v[col])
        cells21.append(dict(score=s, endpoint=name, stat=f'{r:+.3f}', p=p, eff=abs(r)))
    print(f'    {name:12s} ' + ' | '.join(f'{LAB21[s]} {stats.spearmanr(v[s], v[col])[0]:+.3f}' for s in SCORES21))

v = d21[d21['COHORT'] == 1].dropna(subset=['updrs3_score'] + SCORES21)
for s in SCORES21:
    X = v[[s]].astype(float).copy(); X[s] = (X[s] - X[s].mean()) / X[s].std()
    m = sm.OLS(v['updrs3_score'].astype(float), sm.add_constant(X)).fit()
    cells21.append(dict(score=s, endpoint='UPDRS-III', stat=f'{m.params[s]:+.2f}', p=m.pvalues[s],
                        eff=min(abs(m.params[s]) / 2, 1)))
print(f'    UPDRS-III    ' + ' | '.join(
    f'{LAB21[s]} beta={sm.OLS(v["updrs3_score"].astype(float), sm.add_constant((v[[s]].astype(float)-v[s].mean())/v[s].std())).fit().params[s]:+.2f}'
    for s in SCORES21))

prod['grp'] = np.where(prod['subgroup'] == 'Hyposmia', 'Hyposmia', 'Non-Hyposmia')
hyp21 = prod[prod['grp'] == 'Hyposmia']
for s in SCORES21:
    vv = hyp21[[s, 'time_years', 'converted']].dropna()
    c = concordance_index(vv['time_years'], vv[s], vv['converted'])
    p = _boot_p21(lambda b: (concordance_index(vv['time_years'].iloc[b], vv[s].iloc[b], vv['converted'].iloc[b])
                             if vv['converted'].iloc[b].sum() >= 3 else None), len(vv))
    cells21.append(dict(score=s, endpoint='PD conversion', stat=f'{c:.3f}', p=p, eff=abs(c - 0.5) * 2))
C21 = pd.DataFrame(cells21)
C21.to_csv('results/repro/fig4_dissociation_stats.csv', index=False)
sec16['staging_dissociation'] = C21.to_dict('records')
print('\n  Summary (native statistic per cell):')
print(C21.pivot(index='score', columns='endpoint', values='stat').to_string())
print('\n  Each component maps onto one axis and is near chance on the other.')
print('  OIS versus the D axis is partly circular and serves as a positive control.')

# Figure 4
fig = plt.figure(figsize=(14.5, 4.8))
gs21 = GridSpec(1, 3, width_ratios=[1, 1, 1.85], wspace=0.30)
for i, (v, y, nm, note) in enumerate([(vS21, yS21, 'S axis (synuclein)', 'CSF alpha-synuclein positive'),
                                      (vD21, yD21, 'D axis (dopaminergic)', 'DaTscan deficit')]):
    ax = fig.add_subplot(gs21[0, i])
    for s in SCORES21:
        fpr, tpr, _ = roc_curve(y, -v[s])
        ax.plot(fpr, tpr, color=COL_3W[s], lw=2.2 if s != 'upsit' else 1.6,
                ls='-' if s != 'upsit' else '--', label=f'{LAB21[s]}   {roc_auc_score(y, -v[s]):.3f}')
    ax.plot([0, 1], [0, 1], color='lightgrey', lw=1, ls=':')
    ax.set_xlabel('1 - specificity'); ax.set_ylabel('Sensitivity' if i == 0 else '')
    ax.legend(frameon=False, fontsize=9, loc='lower right', title='AUC', title_fontsize=8.5)
    ax.set_title(f'{"ab"[i]}  {nm}\n   {note}, n={len(v)}', fontsize=10.5, loc='left', fontweight='bold')
ax = fig.add_subplot(gs21[0, 2])
eps21 = ['S axis', 'D axis', 'NSD stage', 'Hoehn-Yahr', 'UPDRS-III', 'PD conversion']
short21 = ['S axis', 'D axis', 'NSD\nstage', 'Hoehn\nYahr', 'UPDRS-III', 'conversion']
M21 = np.array([[float(C21[(C21.score == s) & (C21.endpoint == e)]['eff'].iloc[0]) for e in eps21] for s in SCORES21])
im = ax.imshow(M21, cmap='RdYlBu_r', vmin=0, vmax=0.75, aspect='auto')
for i, s in enumerate(SCORES21):
    for j, e in enumerate(eps21):
        r = C21[(C21.score == s) & (C21.endpoint == e)].iloc[0]
        star = '***' if r['p'] < 0.001 else ('**' if r['p'] < 0.01 else ('*' if r['p'] < 0.05 else 'ns'))
        ax.text(j, i - 0.10, r['stat'], ha='center', va='center', fontsize=9.5, fontweight='bold',
                color='white' if M21[i, j] > 0.45 else 'black')
        ax.text(j, i + 0.22, star, ha='center', va='center', fontsize=8,
                color='white' if M21[i, j] > 0.45 else '#555')
ax.set_xticks(range(len(eps21))); ax.set_xticklabels(short21, fontsize=9)
ax.set_yticks(range(3)); ax.set_yticklabels([LAB21[s] for s in SCORES21], fontsize=10.5)
for i in range(3): ax.axhline(i + 0.5, color='white', lw=2)
for j in range(len(eps21)): ax.axvline(j + 0.5, color='white', lw=2)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label('discrimination strength', fontsize=8.5)
ax.set_title('c  Each component maps onto one axis\n   cell shows the native statistic',
             fontsize=10.5, loc='left', fontweight='bold')
plt.savefig(f'{FIG_DIR}/Fig4_v5_double_dissociation.png', dpi=250, bbox_inches='tight')
plt.savefig(f'{FIG_DIR}/Fig4_v5_double_dissociation.pdf', bbox_inches='tight')
plt.show()
print('\n-> Fig4_v5_double_dissociation saved')

# ============================================================================
# Section 22.1. CSF tau and neurofilament, three-way head to head
# ============================================================================
# Section 22.1 (biological interpretability): CSF core biomarkers, three-way head to head.
# Cohort. PPMI baseline visit. Prodromal (COHORT==4) is the primary sample, n=827 with a
#   complete tau panel, and the PD cohort (COHORT==1, n=656) is an independent replication.
# Model. analyte(z) ~ score(z) + age + sex + education, ordinary least squares, so every beta is
#   in analyte SD per score SD. All analytes are log transformed before z scoring.
# Multiplicity. Benjamini-Hochberg across the analytes tested inside each cohort, reported as q.
# Scores. OIS and OMI are the frozen values of section 4, UPSIT = OIS + OMI by construction.
print('22.1 CSF core biomarkers, UPSIT versus OIS versus OMI\n')

if 'PROD22' not in globals():
    _sc = df[['PATNO', 'COHORT', 'OIS', 'OMI']].copy() if 'df' in globals() \
        else pd.read_csv(OUT / 'all_mismatch_scores.csv')
    _sc['PATNO'] = _sc['PATNO'].astype(int)
    _sc['UPSIT'] = _sc['OIS'] + _sc['OMI']
    _cols22 = ['PATNO', 'EVENT_ID', 'YEAR', 'age', 'SEX', 'EDUCYRS', 'subgroup', 'APOE', 'APOE_e4', 'moca']
    _cur22 = (excel[_cols22].copy() if 'excel' in globals()
              else pd.read_excel(CLIN_FILE,
                                 sheet_name=CLIN_SHEET, usecols=_cols22))
    _cur22['PATNO'] = _cur22['PATNO'].astype(int)
    _bl22 = _cur22[_cur22['EVENT_ID'] == 'BL'].drop_duplicates('PATNO').set_index('PATNO')
    _fst22 = _cur22.sort_values(['PATNO', 'YEAR']).drop_duplicates('PATNO').set_index('PATNO')
    _bl22 = _bl22.combine_first(_fst22.loc[[p for p in _fst22.index if p not in _bl22.index]])
    D22 = _sc.merge(_bl22.reset_index()[['PATNO', 'subgroup', 'age', 'SEX', 'EDUCYRS',
                                         'APOE', 'APOE_e4', 'moca']], on='PATNO', how='left')
    D22['sex_m'] = pd.to_numeric(D22['SEX'], errors='coerce')
    D22['educ'] = pd.to_numeric(D22['EDUCYRS'], errors='coerce')
    D22['age'] = pd.to_numeric(D22['age'], errors='coerce')
    D22['APOE_e4_n'] = pd.to_numeric(D22['APOE_e4'], errors='coerce')
    _sg22 = D22['subgroup'].astype(str)
    D22['enr_hyposmia'] = _sg22.str.contains('Hyposmia').astype(int)
    D22['enr_rbd'] = _sg22.str.contains('RBD').astype(int)
    D22['enr_genetic'] = _sg22.str.contains('GBA|LRRK2|SNCA|PRKN|PINK1|VPS35|Other GV').astype(int)
    D22['enr_class'] = np.where(D22.enr_genetic == 1, 'Genetic',
                        np.where(D22.enr_rbd == 1, 'RBD',
                          np.where(D22.enr_hyposmia == 1, 'Hyposmia', 'Other')))
    PROD22 = D22[D22['COHORT'] == 4].copy()
    PDC22 = D22[D22['COHORT'] == 1].copy()
    ENR_SG = ['enr_hyposmia', 'enr_rbd', 'enr_genetic']

    def z22(s):
        s = pd.to_numeric(s, errors='coerce')
        return (s - s.mean()) / s.std(ddof=1)

    def three_way22(frame, ycol, extra=None, min_n=25):
        # y(z) ~ score(z) + age + sex + educ [+ extra], one model per score, plus a joint OIS+OMI model
        extra = extra or []
        base = ['age', 'sex_m', 'educ'] + extra
        res = {}
        for s in ['UPSIT', 'OIS', 'OMI']:
            d = frame[[ycol, s] + base].apply(pd.to_numeric, errors='coerce').dropna()
            if len(d) < min_n:
                res[s] = dict(n=len(d), beta=np.nan, se=np.nan, p=np.nan)
                continue
            X = pd.concat([z22(d[s]), d[base]], axis=1)
            X = X.loc[:, X.std(ddof=0) > 0]
            X = sm.add_constant(X)
            m = sm.OLS(z22(d[ycol]), X).fit()
            res[s] = dict(n=len(d), beta=m.params[s], se=m.bse[s], p=m.pvalues[s])
        d = frame[[ycol, 'OIS', 'OMI'] + base].apply(pd.to_numeric, errors='coerce').dropna()
        if len(d) >= min_n:
            X = pd.concat([z22(d['OIS']).rename('OIS'), z22(d['OMI']).rename('OMI'), d[base]], axis=1)
            X = X.loc[:, X.std(ddof=0) > 0]
            m = sm.OLS(z22(d[ycol]), sm.add_constant(X)).fit()
            res['JOINT'] = dict(n=len(d), b_ois=m.params['OIS'], p_ois=m.pvalues['OIS'],
                                b_omi=m.params['OMI'], p_omi=m.pvalues['OMI'])
        else:
            res['JOINT'] = dict(n=len(d), b_ois=np.nan, p_ois=np.nan, b_omi=np.nan, p_omi=np.nan)
        return res

    print(f'  scores merged with curated baseline: prodromal n={len(PROD22)}, PD n={len(PDC22)}')
    print('  prodromal enrichment class:',
          PROD22['enr_class'].value_counts().to_dict(), '\n')

# --- analytes from the biospecimen file, baseline draw, log scale --------------------------
BIO22 = BIO_FILE
WANT22 = {'ABeta42_CSF': ('ABeta42', 'CSF', 283),
          'pTau181_CSF': ('pTau181', 'CSF', 283),
          'eMTBR_TAU243_CSF': ('eMTBR-TAU243', 'CSF', 307),
          'NfL_CSF': ('NFL', 'CSF', 152),
          'NfL_serum': ('NfL', 'SERUM', 144),
          'NfL_plasma': ('NFL', 'PLASMA', 283)}
_names22 = set(v[0] for v in WANT22.values())
_ch = [c[c['TESTNAME'].isin(_names22)]
       for c in pd.read_csv(BIO22, low_memory=False, chunksize=400000)]
raw22 = pd.concat(_ch, ignore_index=True)
TYPEMAP22 = {'CSF': 'CSF', 'CEREBROSPINAL FLUID': 'CSF', 'PLASMA': 'PLASMA',
             'SERUM': 'SERUM', 'WHOLE BLOOD': 'WHOLE BLOOD'}
raw22['TYPE2'] = raw22['TYPE'].astype(str).str.upper().map(TYPEMAP22)
raw22['PATNO'] = pd.to_numeric(raw22['PATNO'], errors='coerce')
raw22 = raw22.dropna(subset=['PATNO'])
raw22['PATNO'] = raw22['PATNO'].astype(int)
# project 283 plasma reports several unit variants of the same assay, keep the Average rows only
raw22 = raw22[~((raw22['PROJECTID'] == 283) & (raw22['TYPE2'] == 'PLASMA') &
                (~raw22['UNITS'].astype(str).str.contains('Average')))]

def baseline_pick22(sub):
    # BL preferred, then SC, then the earliest remaining draw
    sub = sub.copy()
    sub['rk'] = sub['CLINICAL_EVENT'].map({'BL': 0, 'SC': 1}).fillna(9)
    return sub.sort_values('rk').drop_duplicates('PATNO')

analyte22 = {}
for key, (tn, ty, pid) in WANT22.items():
    s = raw22[(raw22['TESTNAME'] == tn) & (raw22['TYPE2'] == ty) & (raw22['PROJECTID'] == pid)]
    if not len(s):
        continue
    s = baseline_pick22(s)
    v = pd.to_numeric(s['TESTVALUE'].astype(str).str.replace('<', '').str.replace('>', '').str.strip(),
                      errors='coerce')
    ser = pd.Series(v.values, index=s['PATNO'].values).dropna()
    ser = ser[~ser.index.duplicated()]
    if (ser > 0).all():
        ser = np.log(ser)
    analyte22[key] = ser
_cm = analyte22['pTau181_CSF'].index.intersection(analyte22['ABeta42_CSF'].index)
analyte22['pTau181_over_ABeta42_CSF'] = analyte22['pTau181_CSF'][_cm] - analyte22['ABeta42_CSF'][_cm]
print('  analytes loaded:', {k: len(v) for k, v in analyte22.items()}, '\n')

def attach22(frame, keys):
    d = frame.copy()
    for k in keys:
        if k in analyte22:
            d[k] = d['PATNO'].map(analyte22[k])
    return d

TAU_KEYS = ['pTau181_over_ABeta42_CSF', 'pTau181_CSF', 'eMTBR_TAU243_CSF', 'ABeta42_CSF']
NFL_KEYS = ['NfL_CSF', 'NfL_serum', 'NfL_plasma']
LAB22 = {'pTau181_over_ABeta42_CSF': 'CSF pTau181/ABeta42', 'pTau181_CSF': 'CSF pTau181',
         'eMTBR_TAU243_CSF': 'CSF eMTBR-TAU243', 'ABeta42_CSF': 'CSF ABeta42',
         'NfL_CSF': 'CSF NfL', 'NfL_serum': 'Serum NfL', 'NfL_plasma': 'Plasma NfL'}

rows22 = []
for frame, coh in [(PROD22, 'Prodromal'), (PDC22, 'PD')]:
    d = attach22(frame, list(analyte22.keys()))
    for k in TAU_KEYS + NFL_KEYS:
        if k not in d.columns or d[k].notna().sum() < 30:
            continue
        r = three_way22(d, k)
        rows22.append(dict(cohort=coh, analyte=LAB22[k], key=k, n=r['UPSIT']['n'],
                           b_UPSIT=r['UPSIT']['beta'], se_UPSIT=r['UPSIT']['se'], p_UPSIT=r['UPSIT']['p'],
                           b_OIS=r['OIS']['beta'], se_OIS=r['OIS']['se'], p_OIS=r['OIS']['p'],
                           b_OMI=r['OMI']['beta'], se_OMI=r['OMI']['se'], p_OMI=r['OMI']['p'],
                           jb_OIS=r['JOINT']['b_ois'], jp_OIS=r['JOINT']['p_ois'],
                           jb_OMI=r['JOINT']['b_omi'], jp_OMI=r['JOINT']['p_omi']))
CSF22 = pd.DataFrame(rows22)
for coh in CSF22['cohort'].unique():
    m = CSF22['cohort'] == coh
    for s in ['UPSIT', 'OIS', 'OMI']:
        CSF22.loc[m, 'q_' + s] = multipletests(CSF22.loc[m, 'p_' + s], method='fdr_bh')[1]

print('(a) Three-way table, adjusted for age, sex and education')
for coh in ['Prodromal', 'PD']:
    t = CSF22[CSF22['cohort'] == coh]
    print(f'\n  --- {coh} ---')
    print(f"    {'analyte':<22}{'n':>6}{'UPSIT beta':>13}{'p':>10}{'OIS beta':>11}{'p':>9}"
          f"{'OMI beta':>11}{'p':>10}{'q(OMI)':>9}")
    for r in t.itertuples():
        print(f'    {r.analyte:<22}{r.n:>6}{r.b_UPSIT:>+13.4f}{r.p_UPSIT:>10.4g}'
              f'{r.b_OIS:>+11.4f}{r.p_OIS:>9.3g}{r.b_OMI:>+11.4f}{r.p_OMI:>10.4g}{getattr(r, "q_OMI"):>9.3g}')
print('\n  Reading: the tau signal sits entirely on OMI, OIS is flat, and UPSIT carries it only')
print('  because UPSIT contains OMI. NfL, a non specific neurodegeneration marker, is null on all')
print('  three scores in every fluid and both cohorts, so this is not a generic sickness effect.\n')

# --- (b) robustness of the primary outcome, CSF pTau181/ABeta42 in prodromals --------------
print('(b) Robustness, CSF pTau181/ABeta42 in prodromals')
dP = attach22(PROD22, list(analyte22.keys()))
rob22 = {}
for lab, extra in [('age+sex+educ', []),
                   ('+ enrichment subgroup', ENR_SG),
                   ('+ MoCA', ['moca']),
                   ('+ APOE e4 dose', ['APOE_e4_n']),
                   ('+ subgroup + APOE e4', ENR_SG + ['APOE_e4_n'])]:
    r = three_way22(dP, 'pTau181_over_ABeta42_CSF', extra=extra)
    rob22[lab] = dict(n=r['OMI']['n'], b_UPSIT=r['UPSIT']['beta'], p_UPSIT=r['UPSIT']['p'],
                      b_OIS=r['OIS']['beta'], p_OIS=r['OIS']['p'],
                      b_OMI=r['OMI']['beta'], p_OMI=r['OMI']['p'])
    print(f"    {lab:<24}n={r['OMI']['n']:<5} UPSIT {r['UPSIT']['beta']:+.4f} (p={r['UPSIT']['p']:.3g})"
          f"  OIS {r['OIS']['beta']:+.4f} (p={r['OIS']['p']:.3g})"
          f"  OMI {r['OMI']['beta']:+.4f} (p={r['OMI']['p']:.3g})")

print('\n    within a single enrichment stratum (age+sex+educ):')
strat22 = {}
for k, gg in dP.groupby('enr_class'):
    if gg['pTau181_over_ABeta42_CSF'].notna().sum() < 40:
        continue
    r = three_way22(gg, 'pTau181_over_ABeta42_CSF')
    strat22[k] = dict(n=r['OMI']['n'], b_OMI=r['OMI']['beta'], p_OMI=r['OMI']['p'],
                      b_OIS=r['OIS']['beta'], p_OIS=r['OIS']['p'],
                      b_UPSIT=r['UPSIT']['beta'], p_UPSIT=r['UPSIT']['p'])
    print(f"      {k:<10}n={r['OMI']['n']:<5} UPSIT {r['UPSIT']['beta']:+.3f} (p={r['UPSIT']['p']:.3g})"
          f"  OIS {r['OIS']['beta']:+.3f} (p={r['OIS']['p']:.3g})"
          f"  OMI {r['OMI']['beta']:+.3f} (p={r['OMI']['p']:.3g})")

dD = attach22(PDC22, list(analyte22.keys()))
rpd = three_way22(dD, 'pTau181_over_ABeta42_CSF')
print(f"\n    PD replication (independent sample) n={rpd['OMI']['n']}: "
      f"UPSIT {rpd['UPSIT']['beta']:+.4f} (p={rpd['UPSIT']['p']:.3g})  "
      f"OIS {rpd['OIS']['beta']:+.4f} (p={rpd['OIS']['p']:.3g})  "
      f"OMI {rpd['OMI']['beta']:+.4f} (p={rpd['OMI']['p']:.3g})")

# --- (c) cross platform consistency, same subjects, two assays -----------------------------
print('\n(c) Cross platform consistency of the tau signal, prodromal baseline')
_nulf = NULISA_CNS_FILE
_nul = pd.read_excel(_nulf, sheet_name='NPQ Values')
_nul = _nul[_nul['SampleType'].astype(str).str.lower() == 'sample']
_nul = _nul.dropna(subset=['PATNO', 'Target', 'NPQ'])
_nul['PATNO'] = _nul['PATNO'].astype(int)
_nul = _nul[_nul['CLINICAL_EVENT'].astype(str) == 'BL']
NULW22 = _nul.pivot_table(index='PATNO', columns='Target', values='NPQ', aggfunc='mean')
dX = dP.copy()
dX['CSF_pTau181_p283'] = dX['pTau181_CSF']
dX['NULISA_pTau181'] = dX['PATNO'].map(NULW22['pTau-181'])
dX['NULISA_pTau217'] = dX['PATNO'].map(NULW22['pTau-217'])
_cc = dX.dropna(subset=['CSF_pTau181_p283', 'NULISA_pTau181'])
rho_pl, p_pl = stats.spearmanr(_cc['CSF_pTau181_p283'], _cc['NULISA_pTau181'])
r_pl = stats.pearsonr(_cc['CSF_pTau181_p283'], _cc['NULISA_pTau181'])[0]
_all = pd.DataFrame({'a': analyte22['pTau181_CSF']}).join(NULW22[['pTau-181']], how='inner').dropna()
rho_all = stats.spearmanr(_all['a'], _all['pTau-181'])[0]
print(f'    immunoassay pTau181 (project 283) versus NULISA pTau-181 (project 282):')
print(f'      prodromal subsample n={len(_cc)}: Spearman rho={rho_pl:.3f}, Pearson r={r_pl:.3f}')
print(f'      all subjects assayed on both n={len(_all)}: Spearman rho={rho_all:.3f}')
print('    The two assays are run on the same subjects, so this is platform consistency of one')
print('    measurement, not an independent replication.')
plat22 = {}
for c in ['CSF_pTau181_p283', 'NULISA_pTau181', 'NULISA_pTau217']:
    r = three_way22(dX, c)
    plat22[c] = dict(n=r['OMI']['n'], b_UPSIT=r['UPSIT']['beta'], p_UPSIT=r['UPSIT']['p'],
                     b_OIS=r['OIS']['beta'], p_OIS=r['OIS']['p'],
                     b_OMI=r['OMI']['beta'], p_OMI=r['OMI']['p'])
    print(f"    {c:<20}n={r['OMI']['n']:<5} UPSIT {r['UPSIT']['beta']:+.4f} (p={r['UPSIT']['p']:.3g})"
          f"  OIS {r['OIS']['beta']:+.4f} (p={r['OIS']['p']:.3g})"
          f"  OMI {r['OMI']['beta']:+.4f} (p={r['OMI']['p']:.3g})")

# --- save ----------------------------------------------------------------------------------
CSF22.to_csv(OUT / 'section22_csf_threeway.csv', index=False)
CSF_TAU_STATS = dict(
    # 标准误一并存入,正文与摘要按 beta +- 1.96*se 给 95% CI(生物医学期刊要求系数带区间)
    main={f'{r.cohort}|{r.key}': dict(n=int(r.n),
                                      b_UPSIT=r.b_UPSIT, se_UPSIT=r.se_UPSIT, p_UPSIT=r.p_UPSIT,
                                      b_OIS=r.b_OIS, se_OIS=r.se_OIS, p_OIS=r.p_OIS,
                                      b_OMI=r.b_OMI, se_OMI=r.se_OMI, p_OMI=r.p_OMI,
                                      q_OMI=getattr(r, 'q_OMI'))
          for r in CSF22.itertuples()},
    robustness_ratio_prodromal=rob22,
    strata_ratio_prodromal=strat22,
    pd_replication_ratio=dict(n=int(rpd['OMI']['n']), b_UPSIT=rpd['UPSIT']['beta'],
                              p_UPSIT=rpd['UPSIT']['p'], b_OIS=rpd['OIS']['beta'],
                              p_OIS=rpd['OIS']['p'], b_OMI=rpd['OMI']['beta'], p_OMI=rpd['OMI']['p']),
    cross_platform=plat22,
    platform_concordance=dict(n_prodromal=int(len(_cc)), spearman_prodromal=float(rho_pl),
                              pearson_prodromal=float(r_pl), n_all=int(len(_all)),
                              spearman_all=float(rho_all)))
with open(OUT / 'section16_stats.json') as f:
    _s = json.load(f)
_s['csf_tau'] = CSF_TAU_STATS
with open(OUT / 'section16_stats.json', 'w') as f:
    json.dump(_s, f, indent=2, default=float)
print('\n-> section22_csf_threeway.csv written; section16_stats.json updated with csf_tau')

# ============================================================================
# Section 13 (stage table). Residual by biological stage in diagnosed patients
# ============================================================================
# Section 13 (R3 clinical utility): OMI as NSD-stage readout
# Panel a: PD internal, OMI by NSD_STAGE (with Not-NSD highlighted)
# Panel b: prodromal Non-Hyposmia, OMI vs REM behavior cross-sectional

# --- A. PD internal: OMI vs NSD_STAGE ---
ord_map = {'Not NSD': 0, '1a': 1, '1b': 1.5, '2a': 2, '2b': 2.5, '3': 3, '4': 4, '5': 5, '6': 6}
bl = excel[excel['EVENT_ID'] == 'BL'].drop_duplicates('PATNO')
bl_nsd = bl[['PATNO', 'NSD_STAGE']].copy()
bl_nsd['nsd_ord'] = bl_nsd['NSD_STAGE'].map(ord_map)
pd_nsd = df[df['COHORT'] == 1][['PATNO', 'OMI', 'OIS']].merge(bl_nsd, on='PATNO', how='inner').dropna(subset=['OMI', 'nsd_ord'])
rho_pd, p_pd = stats.spearmanr(pd_nsd['OMI'], pd_nsd['nsd_ord'])
rho_pd_ois, p_pd_ois = stats.spearmanr(pd_nsd['OIS'], pd_nsd['nsd_ord'])
print(f'A. PD OMI vs NSD_STAGE: Spearman r={rho_pd:+.3f}, p={p_pd:.1e} | OIS r={rho_pd_ois:+.3f}, p={p_pd_ois:.1e}')
print('   Mean OMI by NSD_STAGE:')
for st, sub in pd_nsd.groupby('NSD_STAGE'):
    if len(sub) >= 5:
        print(f'     {st:8s}: N={len(sub):4d}, mean OMI={sub["OMI"].mean():+6.2f} (sd {sub["OMI"].std():.1f})')

_nsd_rows = [dict(stage=str(st_), n=int(len(sub)), mean_omi=float(sub['OMI'].mean()), sd_omi=float(sub['OMI'].std()))
             for st_, sub in pd_nsd.groupby('NSD_STAGE') if len(sub) >= 5]
sec16['nsd_stage_omi'] = dict(n=int(len(pd_nsd)), spearman_omi=float(rho_pd), p_omi=float(p_pd),
                              spearman_ois=float(rho_pd_ois), p_ois=float(p_pd_ois), rows=_nsd_rows)
with open('results/repro/section16_stats.json') as f:
    _m13 = json.load(f)
_m13.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m13, f, indent=2, default=float)
print('-> ledger key nsd_stage_omi (values previously copied by hand into the Supplementary generator)')

# ============================================================================
# Section 26. Choice of training cohort
# ============================================================================
# Section 26: why the training cohort must be PD + HC.
# A reviewer will ask why we do not train a normative model on controls alone. We answer with data.
from sklearn.model_selection import cross_val_predict as _cvp26
from sklearn.pipeline import make_pipeline as _mkpipe26
def _z(v):
    v = np.asarray(v, dtype=float)
    return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)
print('Section 26  Training-cohort comparison\n')
pool26 = df[df['COHORT'].isin([1, 2])].dropna(subset=FEATURES + ['upsit']).copy()
hyp26 = prod[prod['subgroup'] == 'Hyposmia'].dropna(subset=FEATURES + ['upsit', 'time_years', 'converted']).copy()
hyp26 = hyp26[(hyp26['converted'] == 1) | (hyp26['time_years'] > 0)]
hyp26['time_years'] = hyp26['time_years']  # already floored at 0.25 for converters upstream; do NOT re-clip censored rows
t26, e26 = hyp26['time_years'].values, hyp26['converted'].values.astype(int)
print(f'  hyposmia survival stratum: n={len(hyp26)}, events={int(e26.sum())}')

def _fit26(tr):
    sc = StandardScaler().fit(tr[FEATURES].values)
    m = Ridge(alpha=1.0).fit(sc.transform(tr[FEATURES].values), tr['upsit'].values)
    return m.predict(sc.transform(hyp26[FEATURES].values))

SETS26 = [('PD + HC (main)', pool26), ('HC only', pool26[pool26['COHORT'] == 2]),
          ('PD only', pool26[pool26['COHORT'] == 1])]
rows26, ois26 = [], {}
print(f"\n  {'training set':<16}{'n':>6}{'5-fold r':>10}{'C(OIS) hyposmia':>18}")
for nm, tr in SETS26:
    cvp = _cvp26(_mkpipe26(StandardScaler(), Ridge(alpha=1.0)), tr[FEATURES].values,
                 tr['upsit'].values, cv=KFold(5, shuffle=True, random_state=41))
    cvr = float(stats.pearsonr(cvp, tr['upsit'].values)[0])
    o = _fit26(tr); ois26[nm] = o
    c = float(concordance_index(t26, o, e26))
    rows26.append(dict(training_set=nm, n=int(len(tr)), cv_r=cvr, c_hyposmia=c))
    print(f'  {nm:<16}{len(tr):>6}{cvr:>10.3f}{c:>18.3f}')

# Downsample PD to the HC sample size: separates range restriction from sample size.
_rng26 = np.random.default_rng(42)
_pd26 = pool26[pool26['COHORT'] == 1]
_n_hc = int((pool26['COHORT'] == 2).sum())
_cs = [float(concordance_index(t26, _fit26(_pd26.sample(_n_hc, random_state=int(_rng26.integers(1e6)))), e26))
       for _ in range(200)]
_cs = np.array(_cs)
print(f'\n  PD downsampled to n={_n_hc}, 200 repeats: C = {_cs.mean():.3f} +/- {_cs.std():.3f}'
      f'  range [{_cs.min():.3f}, {_cs.max():.3f}]')
print(f'    {(_cs < 0.55).sum()}/200 below 0.55; the HC-only value {rows26[1]["c_hyposmia"]:.3f} sits at '
      f'percentile {(_cs < rows26[1]["c_hyposmia"]).mean()*100:.1f}')
print(f'  r(OIS_PDonly, OIS_main) on the hyposmia stratum = '
      f'{stats.pearsonr(ois26["PD only"], ois26["PD + HC (main)"])[0]:.3f}')

# Within-cohort smell-to-imaging association, linear and rank.
print('\n  Within-cohort association between UPSIT and striatal binding (why HC-only fails):')
print(f"    {'cohort':<10}{'n':>6}{'putamen r':>11}{'putamen rho':>13}{'caudate rho':>13}{'UPSIT SD':>10}")
for coh, nm in [(2, 'HC'), (1, 'PD'), (4, 'Prodromal')]:
    g = df[df['COHORT'] == coh].dropna(subset=['upsit', 'PUTAMEN_REF_CWM', 'CAUDATE_REF_CWM'])
    print(f"    {nm:<10}{len(g):>6}{stats.pearsonr(g['upsit'], g['PUTAMEN_REF_CWM'])[0]:>11.3f}"
          f"{stats.spearmanr(g['upsit'], g['PUTAMEN_REF_CWM'])[0]:>13.3f}"
          f"{stats.spearmanr(g['upsit'], g['CAUDATE_REF_CWM'])[0]:>13.3f}{g['upsit'].std():>10.2f}")
_hc26 = pool26[pool26['COHORT'] == 2]
print(f'\n  HC range restriction: UPSIT >= 36 in {(_hc26["upsit"] >= 36).mean()*100:.1f}% of controls;'
      f' putamen range [{_hc26["PUTAMEN_REF_CWM"].min():.2f}, {_hc26["PUTAMEN_REF_CWM"].max():.2f}]')
_out26 = ((hyp26['PUTAMEN_REF_CWM'] < _hc26['PUTAMEN_REF_CWM'].min()) |
          (hyp26['PUTAMEN_REF_CWM'] > _hc26['PUTAMEN_REF_CWM'].max())).mean()
print(f'  {_out26*100:.1f}% of the hyposmia stratum falls outside the HC imaging range (extrapolation)')

sec16['training_cohort'] = dict(models=rows26, downsample_mean=float(_cs.mean()), downsample_sd=float(_cs.std()),
    downsample_below_055=int((_cs < 0.55).sum()), n_downsample=_n_hc,
    r_pdonly_vs_main=float(stats.pearsonr(ois26['PD only'], ois26['PD + HC (main)'])[0]),
    hc_extrapolation_fraction=float(_out26),
    verdict='HC-only training yields a chance-level OIS; the failure is range restriction, not sample size')
print('\n-> sec16["training_cohort"] recorded')

# ============================================================================
# Section 27. Prodromal seed amplification
# ============================================================================
# Section 27: prodromal synuclein seed amplification. The data exist; both assays are unusable.
print('Section 27  Prodromal alpha-synuclein seed amplification\n')
_bio = pd.read_csv(BIO_FILE, low_memory=False)
_base27 = df.loc[df['COHORT'] == 4, ['PATNO', 'OIS', 'OMI', 'upsit', 'age', 'sex_num', 'educyrs', 'subgroup']].copy()
_base27['arm'] = np.where(_base27['subgroup'] == 'Hyposmia', 'Hyposmia',
                  np.where(_base27['subgroup'].astype(str).str.contains('RBD'), 'RBD', 'Genetic'))
_base27 = _base27.dropna(subset=['OIS', 'OMI', 'upsit', 'age', 'sex_num', 'educyrs'])
saa27 = {}

# (a) Assay quality check first. An outcome variable this noisy cannot support a negative conclusion.
print('  (a) Assay sensitivity check, skin synSAA (project 259)')
print(f"      {'assay':<26}{'cohort':<12}{'n':>5}{'positive':>10}{'rate':>8}")
_qual = []
for t in ['skin_synSAA R&D-v1 24h', 'skin_synSAA R&D-v2 24h']:
    g = _bio[_bio['TESTNAME'] == t].drop_duplicates('PATNO')
    g = g[g['TESTVALUE'].isin(['Positive', 'Negative'])]
    for coh in ['PD', 'Control', 'Prodromal']:
        s = g[g['COHORT'] == coh]
        if not len(s):
            continue
        pos = int((s['TESTVALUE'] == 'Positive').sum())
        _qual.append(dict(assay=t, cohort=coh, n=int(len(s)), positive=pos, rate=pos / len(s)))
        print(f'      {t:<26}{coh:<12}{len(s):>5}{pos:>10}{pos/len(s)*100:>7.0f}%')
print('      Reference: skin alpha-synuclein detection is positive in 92.7% of PD and 3.3% of controls')
print('      (Gibbons CH, et al. JAMA. 2024;331:1298-1306). The PPMI research-grade assay reaches')
print('      roughly half that sensitivity, and prodromal positivity exceeds the PD cohort, which')
print('      is the wrong ordering. A null result on this outcome is uninformative.')
saa27['assay_quality'] = _qual

# (b) Skin assay association, reported despite the quality problem.
print('\n  (b) Skin synSAA association (project 259)')
g = _bio[_bio['TESTNAME'] == 'skin_synSAA R&D-v2 24h'].drop_duplicates('PATNO')
g = g[g['TESTVALUE'].isin(['Positive', 'Negative'])].copy()
g['saa'] = (g['TESTVALUE'] == 'Positive').astype(int)
m27 = _base27.merge(g[['PATNO', 'saa']], on='PATNO', how='inner')
print(f"      n={len(m27)}, positive={int(m27['saa'].sum())} ({m27['saa'].mean()*100:.1f}%), "
      f"arms={dict(m27['arm'].value_counts())}")
_sk = []
for sc in ['upsit', 'OIS', 'OMI']:
    a = float(roc_auc_score(m27['saa'], -m27[sc]))
    X = sm.add_constant(pd.DataFrame({'s': _z(m27[sc]), 'age': _z(m27['age']),
                                      'sex': m27['sex_num'].values, 'edu': _z(m27['educyrs'])}))
    r = sm.Logit(m27['saa'].values, X).fit(disp=0)
    _sk.append(dict(score=sc, auc=a, odds_ratio=float(np.exp(r.params['s'])), p=float(r.pvalues['s'])))
    print(f'      {sc:<7} AUC={a:.3f}  adjusted OR={np.exp(r.params["s"]):.3f}  p={r.pvalues["s"]:.4g}')
saa27['skin'] = dict(n=int(len(m27)), n_positive=int(m27['saa'].sum()), results=_sk)

# (c) CSF dilution series, project 262.
print('\n  (c) CSF SAA dilution series (project 262)')
_st = _bio[(_bio['PROJECTID'] == 262) & (_bio['TESTNAME'].astype(str).str.endswith('_status'))
           & (_bio['CLINICAL_EVENT'] == 'BL')].copy()
_st['v'] = pd.to_numeric(_st['TESTVALUE'], errors='coerce')
_piv = _st.pivot_table(index='PATNO', columns='TESTNAME', values='v', aggfunc='max')
_piv['saa'] = (_piv.max(axis=1) > 0).astype(int)
m27b = df.merge(_piv[['saa']].reset_index(), on='PATNO', how='inner').dropna(subset=['OIS', 'OMI', 'upsit'])
_pd27 = m27b[m27b['COHORT'] == 1]
print(f"      PD positivity under this assay: {_pd27['saa'].mean()*100:.0f}% ({int(_pd27['saa'].sum())}/{len(_pd27)});"
      f" the main PPMI assay gives 86%. Different sensitivity setting, not comparable.")
p27 = m27b[m27b['COHORT'] == 4].copy()
p27['arm'] = np.where(p27['subgroup'] == 'Hyposmia', 'Hyposmia',
             np.where(p27['subgroup'].astype(str).str.contains('RBD'), 'RBD', 'Genetic'))
p27 = p27.dropna(subset=['age', 'sex_num', 'educyrs'])
print(f"      prodromal n={len(p27)}, positive={int(p27['saa'].sum())}, arms={dict(p27['arm'].value_counts())}")
_csf = []
for sc in ['upsit', 'OIS', 'OMI']:
    a = float(roc_auc_score(p27['saa'], -p27[sc]))
    base_cols = {'s': _z(p27[sc]), 'age': _z(p27['age']),
                 'sex': p27['sex_num'].values, 'edu': _z(p27['educyrs'])}
    r1 = sm.Logit(p27['saa'].values, sm.add_constant(pd.DataFrame(base_cols))).fit(disp=0)
    adj = dict(base_cols); adj['rbd'] = (p27['arm'] == 'RBD').astype(float).values
    r2 = sm.Logit(p27['saa'].values, sm.add_constant(pd.DataFrame(adj))).fit(disp=0)
    _csf.append(dict(score=sc, auc=a, p_unadjusted=float(r1.pvalues['s']), p_arm_adjusted=float(r2.pvalues['s'])))
    print(f'      {sc:<7} AUC={a:.3f}  p={r1.pvalues["s"]:.4g}  ->  arm-adjusted p={r2.pvalues["s"]:.4g}')
_hyp27 = p27[p27['arm'] == 'Hyposmia']
print(f"      within the hyposmia arm alone: n={len(_hyp27)}, positive={int(_hyp27['saa'].sum())}")

# Power for the single clean channel, simulated.
def _power27(n, rate, auc, nsim=2000, alpha=0.05):
    d = stats.norm.ppf(auc) * np.sqrt(2); hit = 0
    rr = np.random.default_rng(7)
    for _ in range(nsim):
        npos = rr.binomial(n, rate)
        if npos < 3 or n - npos < 3:
            continue
        x = np.r_[rr.normal(d, 1, npos), rr.normal(0, 1, n - npos)]
        yy = np.r_[np.ones(npos), np.zeros(n - npos)]
        a = roc_auc_score(yy, x)
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        n1, n2 = npos, n - npos
        se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2))
        if abs(a - 0.5) / se > stats.norm.ppf(1 - alpha / 2):
            hit += 1
    return hit / nsim
_pw = _power27(len(_hyp27), _hyp27['saa'].mean(), 0.70)
print(f'      power to detect AUC 0.70 in that channel: {_pw*100:.0f}%')
saa27['csf'] = dict(n=int(len(p27)), n_positive=int(p27['saa'].sum()), results=_csf,
                    pd_positivity=float(_pd27['saa'].mean()), hyposmia_arm_n=int(len(_hyp27)),
                    hyposmia_arm_positive=int(_hyp27['saa'].sum()), hyposmia_arm_power_auc070=float(_pw))
saa27['verdict'] = ('prodromal seed-amplification data exist in the biospecimen file but neither assay '
                    'can support a conclusion: the skin assay lacks sensitivity, the CSF dilution series '
                    'lacks power once the enrolment channel is controlled')
sec16['prodromal_saa'] = saa27
print('\n-> sec16["prodromal_saa"] recorded')

# ============================================================================
# Section 28. Genotype premise
# ============================================================================
# Section 28: the premise behind the genotype test, established inside this cohort.
from scipy.stats import fisher_exact as _fisher
print('Section 28  Synuclein burden differs by genotype (premise for section 20)\n')
_bl28 = excel[excel['EVENT_ID'] == 'BL'][['PATNO', 'Stage_S']].drop_duplicates('PATNO')
m28 = df.merge(_bl28, on='PATNO', how='inner')
m28 = m28[m28['Stage_S'].isin([0, 1])].copy()
m28['saa'] = (m28['Stage_S'] == 1).astype(int)
_pd28 = m28[m28['COHORT'] == 1]
print(f"  {'group':<16}{'n':>6}{'SAA positive':>14}{'rate':>8}")
rows28 = []
for sg, lab in [('GBA', 'GBA'), ('LRRK2', 'LRRK2'), ('Sporadic PD', 'Sporadic')]:
    g = _pd28[_pd28['subgroup'] == sg]
    if not len(g):
        continue
    pos = int(g['saa'].sum())
    rows28.append(dict(group=lab, n=int(len(g)), positive=pos, rate=float(pos / len(g))))
    print(f'  {lab:<16}{len(g):>6}{pos:>14}{pos/len(g)*100:>7.0f}%')
_g, _l = [r for r in rows28 if r['group'] == 'GBA'][0], [r for r in rows28 if r['group'] == 'LRRK2'][0]
_p28 = _fisher([[_g['positive'], _g['n'] - _g['positive']], [_l['positive'], _l['n'] - _l['positive']]])[1]
print(f'\n  GBA versus LRRK2, Fisher exact p = {_p28:.3g}')
print('  Direction agrees with the literature: LRRK2 carriers show relatively preserved olfaction')
print('  (36% hyposmic versus 75% in non-carriers, Ruiz-Martinez et al. Mov Disord. 2011;26:2026-2031).')
print('  Section 20 therefore tests a prediction whose direction was fixed in advance.')
sec16['genotype_saa_premise'] = dict(groups=rows28, fisher_p=float(_p28))
print('\n-> sec16["genotype_saa_premise"] recorded')

# ============================================================================
# Section 30. Cognitive endpoints
# ============================================================================
# Section 30: cognitive endpoints. Univariate looks promising; multivariable shows no increment.
print('Section 30  Cognitive endpoints, univariate and multivariable\n')
_L30 = excel[excel['COHORT'] == 4][['PATNO', 'age_at_visit', 'cogstate', 'MCI_testscores']].copy()
_L30['t_yr'] = _L30['age_at_visit'] - _L30.groupby('PATNO')['age_at_visit'].transform('min')
_b30 = df.loc[df['COHORT'] == 4, ['PATNO', 'OIS', 'OMI', 'upsit', 'age', 'sex_num', 'educyrs', 'subgroup']].copy()
_b30['arm'] = np.where(_b30['subgroup'] == 'Hyposmia', 'Hyposmia',
              np.where(_b30['subgroup'].astype(str).str.contains('RBD'), 'RBD', 'Genetic'))
for c in ['age', 'sex_num', 'educyrs']:
    _b30[c] = _b30[c].fillna(_b30[c].median())

def _incident30(var, bl_ok, ev_ok):
    d = _L30[['PATNO', 't_yr', var]].dropna().sort_values(['PATNO', 't_yr'])
    rows = []
    for p, g in d.groupby('PATNO'):
        g = g.sort_values('t_yr')
        if not bl_ok(g.iloc[0][var]):
            continue
        fu = g.iloc[1:]
        if not len(fu):
            continue
        ev = fu[fu[var].map(ev_ok)]
        if len(ev):
            rows.append(dict(PATNO=p, time=max(ev.iloc[0]['t_yr'], 0.25), event=1))
        elif fu['t_yr'].max() > 0:
            rows.append(dict(PATNO=p, time=fu['t_yr'].max(), event=0))
    return pd.DataFrame(rows).merge(_b30, on='PATNO').dropna(subset=['OIS', 'OMI', 'upsit'])

def _cox30(s, cols, arms=True):
    X = pd.DataFrame({'time': s['time'].values, 'event': s['event'].values, 'age': _z(s['age']),
                      'sex': s['sex_num'].values, 'edu': _z(s['educyrs'])})
    if arms:
        for lev in ['RBD', 'Genetic']:
            X[lev] = (s['arm'] == lev).astype(float).values
    for c in cols:
        X[c] = _z(s[c])
    m = CoxPHFitter().fit(X, 'time', 'event')
    return m

rows30 = []
for var, bl, ev, lab in [('MCI_testscores', lambda v: v == 0, lambda v: v == 1, 'incident MCI (neuropsych)'),
                         ('cogstate', lambda v: v == 1, lambda v: v >= 2, 'incident cognitive impairment (clinician)')]:
    s30 = _incident30(var, bl, ev)
    print(f'  {lab}:  at risk n={len(s30)}, events={int(s30["event"].sum())}')
    print(f"    {'model':<26}{'result'}")
    for cols, nm in [(['upsit'], 'univariate UPSIT'), (['OIS'], 'univariate OIS'), (['OMI'], 'univariate OMI'),
                     (['OIS', 'OMI'], 'OIS + OMI'), (['upsit', 'OMI'], 'UPSIT + OMI (increment)')]:
        m = _cox30(s30, cols)
        parts = []
        for c in cols:
            hr, p = float(np.exp(m.params_[c])), float(m.summary.loc[c, 'p'])
            parts.append(f'{c} HR={hr:.3f} p={p:.4g}')
            rows30.append(dict(endpoint=lab, model=nm, term=c, hr=hr, p=p,
                               concordance=float(m.concordance_index_), n=int(len(s30)),
                               events=int(s30['event'].sum())))
        print(f'    {nm:<26}' + '  |  '.join(parts) + f'   C={m.concordance_index_:.3f}')
    print()
print('  Reading the two endpoints together. On the neuropsychology-based MCI endpoint the two')
print('  components have indistinguishable effects (residual HR 1.045, p=0.81 once the total score is')
print('  in the model), so the equal-weight sum is already optimal and the decomposition adds nothing.')
print('  On the clinician-assessed endpoint the two components differ significantly (p<0.001), but the')
print('  signal sits in the imaging-explained component (OIS HR 0.670, p<0.001) rather than in the')
print('  residual (HR 0.873, p=0.054). Neither endpoint supports an independent claim for the residual.')
print('  Note that in a model containing the total score, the coefficient on the residual estimates the')
print('  difference between the two component effects rather than an incremental effect of its own.')
sec16['cognitive_endpoints'] = rows30
print('-> sec16["cognitive_endpoints"] recorded')

# ============================================================================
# Section 32. Age in the prodromal residual
# ============================================================================
# Section 32: age absorption transfers imperfectly to the application cohort.
print('Section 32  Residual age dependence in the application cohort\n')
print(f"  {'cohort':<20}{'n':>6}{'r(age, UPSIT)':>16}{'r(age, OIS)':>14}{'r(age, OMI)':>14}{'p(OMI)':>10}")
rows32 = []
for coh, nm in [(1, 'PD (training)'), (2, 'HC (training)'), (4, 'Prodromal (applied)')]:
    g = df[df['COHORT'] == coh].dropna(subset=['upsit', 'OIS', 'OMI', 'age'])
    ru = stats.pearsonr(g['age'], g['upsit'])[0]
    ro = stats.pearsonr(g['age'], g['OIS'])[0]
    rm, pm = stats.pearsonr(g['age'], g['OMI'])
    rows32.append(dict(cohort=nm, n=int(len(g)), r_age_upsit=float(ru), r_age_ois=float(ro),
                       r_age_omi=float(rm), p_age_omi=float(pm)))
    print(f'  {nm:<20}{len(g):>6}{ru:>16.3f}{ro:>14.3f}{rm:>14.3f}{pm:>10.3g}')
_h32 = df[(df['COHORT'] == 4) & (df['subgroup'] == 'Hyposmia')].dropna(subset=['OMI', 'age'])
_rh, _ph = stats.pearsonr(_h32['age'], _h32['OMI'])
rows32.append(dict(cohort='Hyposmia stratum', n=int(len(_h32)), r_age_upsit=None, r_age_ois=None,
                   r_age_omi=float(_rh), p_age_omi=float(_ph)))
print(f'  {"Hyposmia stratum":<20}{len(_h32):>6}{"":>16}{"":>14}{_rh:>14.3f}{_ph:>10.3g}')
print('\n  Age is a model input, so the residual is age-independent inside the training cohort by')
print('  construction. That adjustment does not transfer fully to the prodromal cohort. Every analysis')
print('  involving the residual adjusts for age, but the quantity itself retains an age component in')
print('  the application cohort, and this belongs in the limitations.')

# Confirm the absence of confounder fields, field by field.
_pats = {'nasal': ['nasal', 'sinus', 'rhin', 'polyp', 'septum'], 'head trauma': ['trauma', 'injury', 'concus', 'tbi'],
         'smoking': ['smok', 'tobacco', 'cigar', 'nicotin'], 'post-infectious': ['infect', 'covid', 'viral', 'urti']}
print('\n  Confounder fields present in the curated cut:')
_absent32 = []
for k, pats in _pats.items():
    hit = [c for c in excel.columns if any(p in c.lower() for p in pats)]
    print(f'    {k:<18}{hit if hit else "none"}')
    if not hit:
        _absent32.append(k)
sec16['age_absorption'] = dict(correlations=rows32, unrecorded_confounders=_absent32)
print('\n-> sec16["age_absorption"] recorded')
# Merge rather than overwrite: several earlier sections write straight to the file and never
# touch the in-memory sec16, so dumping sec16 alone would silently drop those keys.
with open('results/repro/section16_stats.json') as f:
    _merged32 = json.load(f)
_before32 = set(_merged32)
_merged32.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_merged32, f, indent=2, default=float)
print(f'-> section16_stats.json merged: {len(_before32)} keys before, {len(_merged32)} after,'
      f' {len(set(_merged32) - _before32)} added, 0 dropped')

# ============================================================================
# Section 42. Partial correlations
# ============================================================================
from scipy import stats as _st42

def _partial_r(y, x, covars):
    # partial correlation between y and x given covars, via residualisation
    C = sm.add_constant(covars.astype(float))
    ry = sm.OLS(np.asarray(y, float), C).fit().resid
    rx = sm.OLS(np.asarray(x, float), C).fit().resid
    r, p = _st42.pearsonr(ry, rx)
    return float(r), float(p), int(len(ry))

rows42 = []

# --- (a) prodromal, CSF pTau181/ABeta42 -----------------------------------------------------
_d42 = attach22(PROD22, ['pTau181_over_ABeta42_CSF'])
_d42 = _d42[['pTau181_over_ABeta42_CSF', 'OIS', 'OMI', 'UPSIT',
             'age', 'sex_m', 'educ']].apply(pd.to_numeric, errors='coerce').dropna()
print(f'(a) Prodromal, CSF pTau181/ABeta42, n = {len(_d42)}')
for sc in ['OMI', 'OIS', 'UPSIT']:
    r, pv, n = _partial_r(_d42[sc], _d42['pTau181_over_ABeta42_CSF'], _d42[['age', 'sex_m', 'educ']])
    print(f'    partial r({sc:5s}, tau ratio) = {r:+.3f}   P = {pv:.3g}   n = {n}')
    rows42.append(dict(test='CSF pTau181/ABeta42', cohort='Prodromal', score=sc,
                       r=r, abs_r=abs(r), p=pv, n=n))

# --- (b) manifest PD, GBA vs LRRK2, point-biserial ------------------------------------------
_g42 = PDC22[PDC22['subgroup'].isin(['GBA', 'LRRK2'])].copy()
_g42['carrier'] = (_g42['subgroup'] == 'GBA').astype(int)
_g42 = _g42[['carrier', 'OIS', 'OMI', 'UPSIT', 'age', 'sex_m', 'educ']].apply(
    pd.to_numeric, errors='coerce').dropna()
print(f'\n(b) Manifest PD, GBA vs LRRK2, n = {len(_g42)} '
      f'(GBA {int(_g42.carrier.sum())}, LRRK2 {int((1 - _g42.carrier).sum())})')
for sc in ['OMI', 'OIS', 'UPSIT']:
    r, pv, n = _partial_r(_g42[sc], _g42['carrier'], _g42[['age', 'sex_m', 'educ']])
    print(f'    partial r({sc:5s}, genotype) = {r:+.3f}   P = {pv:.3g}   n = {n}')
    rows42.append(dict(test='GBA vs LRRK2 genotype', cohort='PD', score=sc,
                       r=r, abs_r=abs(r), p=pv, n=n))

print()
print('  Reading: on a single common scale the residual carries the association in both tests')
print('  while the fitted value sits near zero, which is the point of the decomposition. The')
print('  genotype contrast is the same model as section 23, expressed as a correlation so that')
print('  the two tests can be plotted on one axis.')

pr42 = pd.DataFrame(rows42)
pr42.to_csv('results/repro/section42_partial_r.csv', index=False)
sec16['partial_r_two_tests'] = pr42.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m42 = json.load(f)
_b42 = set(_m42); _m42.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m42, f, indent=2, default=float)
print(f'-> results/repro/section42_partial_r.csv written')
print(f'-> section16_stats.json: {len(_b42)} -> {len(_m42)} keys')

# ============================================================================
# Section 45. Tau tertile means (Fig. 1d)
# ============================================================================
from scipy import stats as _st45
print('Section 45  Residual and fitted value across CSF tau burden, as drawn in Fig 1d\n')
_d45 = attach22(PROD22, ['pTau181_over_ABeta42_CSF'])
_d45 = _d45[['pTau181_over_ABeta42_CSF', 'OIS', 'OMI', 'UPSIT', 'age', 'sex_m', 'educ']].apply(
    pd.to_numeric, errors='coerce').dropna()
print(f'  prodromal participants with a CSF pTau181/ABeta42 ratio: n = {len(_d45)}')
_C45 = sm.add_constant(_d45[['age', 'sex_m', 'educ']].astype(float))
for c in ['OMI', 'OIS', 'UPSIT']:
    _d45[c + '_adj'] = sm.OLS(_d45[c].values, _C45).fit().resid + _d45[c].mean()
rows45 = []
for k, lab in [(3, 'tertile'), (2, 'median split')]:
    _d45['g'] = pd.qcut(_d45['pTau181_over_ABeta42_CSF'], k, labels=[f'G{i+1}' for i in range(k)])
    groups = [(str(nm), sub) for nm, sub in _d45.groupby('g', observed=True)]
    print(f'\n  --- CSF tau {lab}, G1 = lowest tau burden (means adjusted for age, sex and education)')
    for sc in ['OMI', 'OIS']:
        pv = _st45.f_oneway(*[sub[sc + '_adj'] for _, sub in groups]).pvalue
        cells = '  '.join(f'{nm}: {sub[sc + "_adj"].mean():+6.2f} (SE {sub[sc + "_adj"].sem():.2f}, n={len(sub)})'
                          for nm, sub in groups)
        print(f'    {sc:5s} {cells}   across-group P = {pv:.3g}')
        for nm, sub in groups:
            rows45.append(dict(split=lab, n_groups=k, group=nm, score=sc, n=int(len(sub)),
                               adj_mean=float(sub[sc + '_adj'].mean()), adj_se=float(sub[sc + '_adj'].sem()),
                               p_across_groups=float(pv),
                               tau_mean=float(sub['pTau181_over_ABeta42_CSF'].mean())))
    d_omi = groups[-1][1]['OMI_adj'].mean() - groups[0][1]['OMI_adj'].mean()
    d_ois = groups[-1][1]['OIS_adj'].mean() - groups[0][1]['OIS_adj'].mean()
    print(f'    highest minus lowest tau group: OMI {d_omi:+.2f} points, OIS {d_ois:+.2f} points')
print()
print('  Reading: the residual falls monotonically as tau burden rises while the fitted value does not move.')
print('  This is the grouped form of the continuous association reported in section 22 and is shown in Fig 1d;')
print('  three groups are used there because the claim is a gradient rather than a split.')
T45 = pd.DataFrame(rows45)
T45.to_csv('results/repro/section45_tau_groups.csv', index=False)
sec16['tau_gradient_fig1d'] = T45.to_dict('records')
with open('results/repro/section16_stats.json') as f:
    _m45 = json.load(f)
_b45 = set(_m45); _m45.update(sec16)
with open('results/repro/section16_stats.json', 'w') as f:
    json.dump(_m45, f, indent=2, default=float)
print(f'-> results/repro/section45_tau_groups.csv written')
print(f'-> section16_stats.json: {len(_b45)} -> {len(_m45)} keys')
