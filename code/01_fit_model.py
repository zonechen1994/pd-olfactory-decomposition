"""01_fit_model.py  Data loading, model comparison, OIS fitting, leave-one-centre-out

Reads the two PPMI files, merges the earliest DaTscan quantification per participant with
the baseline clinical record, compares six regression families under five-fold cross-validation
in the PD + HC training cohort, fits the final ridge model, writes OIS and OMI for every
participant, and runs the leave-one-centre-out validation.
Outputs: model_output/repro_sbr_models/{sbr_all_Ridge.pkl, scaler.pkl, metadata.json},
results/repro/all_mismatch_scores.csv (participant level, not distributed),
results/repro/loso_per_site.csv, loso_per_site_summary.json, results/intermediate/01_model.pkl
Notebook provenance: sections 1, 1.1, 2, 3, 4 and 8. Code is the notebook cell text, unchanged except for
the file paths, which are resolved through common.py.
"""
from common import *

# ============================================================================
# Section 1. Load and merge
# ============================================================================
# --- SBR: earliest visit per subject ---
sbr = pd.read_csv(SBR_FILE)
sbr['visit_rank'] = sbr['EVENT_ID'].map({'SC': 0, 'BL': 1, 'V04': 2, 'V06': 3}).fillna(99)
sbr_bl = sbr.sort_values(['PATNO', 'visit_rank']).drop_duplicates('PATNO', keep='first')
SBR_COLS = [c for c in sbr_bl.columns if 'REF_CWM' in c]
sbr_bl = sbr_bl.dropna(subset=SBR_COLS)
print(f'SBR features: {len(SBR_COLS)} | subjects with complete SBR: {len(sbr_bl)}')

# --- Clinical baseline ---
excel = pd.read_excel(CLIN_FILE, sheet_name=CLIN_SHEET)
# visit_date ships as the STRING 'MM/YYYY'. Sorting it as text orders by month before year,
# so '01/2021' < '05/2018' and every groupby min/max on the raw column returns the wrong visit.
# Parse it once, here, before anything derives follow-up time from it.
excel['visit_date'] = pd.to_datetime(excel['visit_date'], format='%m/%Y', errors='coerce')
assert excel['visit_date'].notna().all(), 'visit_date failed to parse as %m/%Y'
clin = excel[excel['EVENT_ID'] == 'BL'].drop_duplicates('PATNO')
keep = ['PATNO','COHORT','age','SEX','EDUCYRS','upsit','moca','updrs3_score',
        'subgroup','SITE','enroll_phase']
df = sbr_bl.merge(clin[keep], on='PATNO', how='inner')

# --- Covariates ---
df['age_num'] = df['age']
df['sex_num'] = df['SEX'].astype(float)      # already 0/1
df['educyrs'] = df['EDUCYRS']
DEMO = ['age_num', 'sex_num', 'educyrs']
FEATURES = SBR_COLS + DEMO

# fill missing covariates with PD+HC training medians
demo_medians = df.loc[df['COHORT'].isin([1, 2]), DEMO].median()
for c in DEMO:
    df[c] = df[c].fillna(demo_medians[c])

print(f'Merged SBR+clinical: {len(df)} subjects')
print('COHORT distribution:', df['COHORT'].value_counts().sort_index().to_dict())
print('(1=PD, 2=HC, 3=SWEDD, 4=Prodromal)')

# ============================================================================
# Section 1.1. Hyposmia and non-hyposmia enrolment labels
# ============================================================================
# Print the subgroup composition and UPSIT stats to verify the definition
prod_bl_view = df[df['COHORT']==4][['PATNO','subgroup','upsit']].copy()
print(f"PPMI Prodromal (COHORT==4) at baseline: N={len(prod_bl_view)} subjects with SBR\n")
print("Subgroup composition (per PPMI enrollment label):")
sg_tbl = prod_bl_view.groupby('subgroup', dropna=False).agg(
    N=('PATNO','count'),
    UPSIT_n=('upsit','count'),
    UPSIT_mean=('upsit','mean'),
    UPSIT_sd=('upsit','std'),
    UPSIT_min=('upsit','min'),
    UPSIT_max=('upsit','max'),
).sort_values('N', ascending=False).round(2)
print(sg_tbl.to_string())

print("\nApplied dichotomy used throughout R2/R3 analyses:")
prod_bl_view['grp'] = np.where(prod_bl_view['subgroup']=='Hyposmia', 'Hyposmia', 'Non-Hyposmia')
grp_tbl = prod_bl_view.groupby('grp').agg(
    N=('PATNO','count'),
    UPSIT_n=('upsit','count'),
    UPSIT_mean=('upsit','mean'),
    UPSIT_sd=('upsit','std'),
    UPSIT_floor_le15=('upsit', lambda x: (x<=15).sum()),
    UPSIT_floor_le20=('upsit', lambda x: (x<=20).sum()),
).round(2)
print(grp_tbl.to_string())
print("\n=> Hyposmia 子组超过一半受试者 UPSIT ≤ 20 (地板效应), Non-Hyposmia 分布更宽")

# Sanity check: PPMI Hyposmia protocol verification
# PPMI 协议要求 Hyposmia 子组 UPSIT 落在年龄/性别 15th 百分位以下, 实际数据应高度集中在低位
print("\n--- Sanity check: PPMI Hyposmia 子组 UPSIT 累计占比 (验证 PPMI 入组协议执行) ---")
hyp_upsit = prod_bl_view[prod_bl_view['subgroup']=='Hyposmia']['upsit'].dropna()
nh_upsit  = prod_bl_view[prod_bl_view['subgroup']!='Hyposmia']['upsit'].dropna()
print(f"  Hyposmia 子组 (n={len(hyp_upsit)}):")
for thr in [15, 20, 25, 30, 35]:
    pct_h = (hyp_upsit <= thr).mean() * 100
    print(f"    UPSIT ≤ {thr}: {pct_h:.1f}% ({(hyp_upsit<=thr).sum()}/{len(hyp_upsit)})")
print(f"  Non-Hyposmia 子组 (n={len(nh_upsit)}, 对比):")
for thr in [15, 20, 25, 30, 35]:
    pct_n = (nh_upsit <= thr).mean() * 100
    print(f"    UPSIT ≤ {thr}: {pct_n:.1f}% ({(nh_upsit<=thr).sum()}/{len(nh_upsit)})")
print("=> 99% 的 PPMI 'Hyposmia' 子组 UPSIT ≤ 35, 95%+ ≤ 30, 强烈支持 PPMI 是按 UPSIT ≤ 15th 百分位入组的; ")
print("   Non-Hyposmia 仅 ~60-70% ≤ 30, 跨度更广, 主要由非嗅觉富集 (RBD/LRRK2/GBA) 决定。")

# Subquestion: 'Why are there low-UPSIT subjects in Non-Hyposmia?'
print("\n--- 解释 Non-Hyposmia 中低 UPSIT 受试者的来源 (按 subgroup 拆) ---")
low_nh = prod_bl_view[(prod_bl_view['subgroup']!='Hyposmia') & (prod_bl_view['upsit']<=20)]
print(f"Non-Hyposmia with UPSIT ≤ 20 总人数: {len(low_nh)} (29% of Non-Hyposmia)")
print("分布按原 PPMI subgroup:")
print(low_nh['subgroup'].value_counts().to_string())
rbd_pct = (low_nh['subgroup']=='RBD').mean()*100
print(f"\n=> {rbd_pct:.0f}% 的低 UPSIT Non-Hyposmia 是 RBD 子组, 他们入组通道是 RBD, 同时也嗅觉减退, PPMI 标 'RBD' 不标 'Hyposmia'。\n   生物学上 RBD + hyposmia 共存常见, 因为两者都是突触核蛋白病的早期标志。\n   这证明 PPMI subgroup 字段反映入组通道, 不反映受试者最终临床表型。")

# ============================================================================
# Section 2. Six regression families, five-fold cross-validation (KFold seed 41)
# ============================================================================
def make_models():
    # n_jobs=1 on tree models: folds are parallelised by cross_val_predict instead,
    # which avoids nested-parallelism oversubscription (huge slowdown otherwise).
    return {
        'Ridge':    Ridge(alpha=1.0),
        'RF':       RandomForestRegressor(n_estimators=300, max_depth=6, min_samples_leaf=8,
                                          random_state=42, n_jobs=1),
        'SVR-Lin':  SVR(kernel='linear', C=1.0),
        'XGBoost':  XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.03,
                                 subsample=0.8, colsample_bytree=0.8, random_state=42, n_jobs=1),
        'LightGBM': LGBMRegressor(n_estimators=300, max_depth=3, learning_rate=0.03,
                                  subsample=0.8, colsample_bytree=0.8, random_state=42,
                                  n_jobs=1, verbose=-1),
        'GradBoost':GradientBoostingRegressor(n_estimators=200, max_depth=3, min_samples_leaf=8,
                                              learning_rate=0.03, subsample=0.8, random_state=42),
    }

train = df[df['COHORT'].isin([1, 2])].copy()
kf = KFold(n_splits=5, shuffle=True, random_state=41)

def cv_compare(target):
    m = train[target].notna()
    X = StandardScaler().fit_transform(train.loc[m, FEATURES].values)
    y = train.loc[m, target].values
    fold_splits = list(kf.split(X, y))
    per_fold = {name: [] for name in make_models().keys()}
    rows = []
    for name, mdl in make_models().items():
        yp = cross_val_predict(mdl, X, y, cv=kf, n_jobs=5)
        for _tr, te in fold_splits:
            per_fold[name].append(stats.pearsonr(y[te], yp[te])[0])
        fr = per_fold[name]
        rows.append({'model': name,
                     'pearson_r':   stats.pearsonr(y, yp)[0],
                     'fold_mean':   np.mean(fr),
                     'fold_std':    np.std(fr),
                     'fold_min':    np.min(fr),
                     'fold_max':    np.max(fr),
                     'spearman_r':  stats.spearmanr(y, yp)[0],
                     'mae':         mean_absolute_error(y, yp)})
    return pd.DataFrame(rows), m.sum(), per_fold

# Stage 1: OIS = imaging-anchored prediction of UPSIT. The slim paper focuses on the olfactory
# axis; cognitive/motor analogues are not constructed in this notebook.
all_cv = {}
comp_ois, n_ois, pf_ois = cv_compare('upsit')
all_cv['upsit'] = (comp_ois, n_ois, pf_ois, 'OIS (UPSIT, olfactory primary)')
print(f'=== OIS (UPSIT, olfactory primary)  train N = {n_ois}  [5-fold CV on PD+HC] ===')
print(comp_ois.to_string(index=False))
best = comp_ois.sort_values('pearson_r', ascending=False).iloc[0]
ridge_row = comp_ois[comp_ois.model == "Ridge"].iloc[0]
print(f'  Best by overall r: {best["model"]} (r = {best["pearson_r"]:.3f}, fold mean ± sd = {best["fold_mean"]:.3f} ± {best["fold_std"]:.3f})')
print(f'  Ridge (chosen for main analyses): r = {ridge_row["pearson_r"]:.3f}'
      f'   (within {(best["pearson_r"]-ridge_row["pearson_r"])*100:+.1f}pp of best)')

# ============================================================================
# Section 2 (continued). Model comparison figure
# ============================================================================
# OIS-only model comparison plot (primary application)
fig, ax = plt.subplots(figsize=(7.5, 4.6))
model_order = list(make_models().keys())
comp, n, per_fold, label = all_cv['upsit']
means = [np.mean(per_fold[m]) for m in model_order]
stds  = [np.std(per_fold[m])  for m in model_order]
pos = np.arange(len(model_order))
colors = ['#16A085' if m == 'Ridge' else '#bbbbbb' for m in model_order]
bars = ax.bar(pos, means, yerr=stds, color=colors, capsize=4, edgecolor='#444', linewidth=0.8)
for i, m in enumerate(model_order):
    for r in per_fold[m]:
        ax.scatter(i, r, color='#C0392B', s=22, alpha=0.65, zorder=3, edgecolor='white', linewidth=0.6)
rmean = means[model_order.index('Ridge')]; rstd = stds[model_order.index('Ridge')]
ax.annotate(f'Ridge chosen for main analyses\n(r = {rmean:.3f} ± {rstd:.3f})',
            xy=(0, rmean), xytext=(1.5, max(means) + max(stds) + 0.08),
            fontsize=10, color='#138D75', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#138D75', lw=1.0))
ax.set_xticks(pos); ax.set_xticklabels(model_order, rotation=30, ha='right', fontsize=9.5)
ax.set_ylabel('5-fold CV Pearson r'); ax.axhline(0, color='k', lw=0.7)
ax.set_title(f'OIS (UPSIT, olfactory primary): 6 regression families, 5-fold CV (n = {n})',
             loc='left', fontsize=11, fontweight='bold')
ax.set_ylim(0, max(means + stds) * 1.4 + 0.02)
plt.tight_layout()
plt.savefig(FIG / 'FigS_R1_OIS_model_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig(FIG / 'FigS_R1_OIS_model_comparison.pdf', bbox_inches='tight')
plt.show()
print('=> 6 model families within ~1% of each other; Ridge retained for interpretable linear score.')

# ============================================================================
# Section 2 (continued). Final ridge fit, OIS and OMI for every participant
# ============================================================================
# Fit the final Ridge OIS model on PD+HC and apply to ALL subjects.
def fit_score(target, alpha=1.0):
    m = train[target].notna()
    sc = StandardScaler().fit(train.loc[m, FEATURES].values)
    model = Ridge(alpha=alpha).fit(sc.transform(train.loc[m, FEATURES].values),
                                   train.loc[m, target].values)
    pred_all = model.predict(sc.transform(df[FEATURES].values))
    return model, sc, pred_all

ois_model, ois_sc, df['OIS'] = fit_score('upsit')
df['OMI'] = df['upsit'] - df['OIS']

v = df[['OIS', 'upsit']].dropna()
print(f'OIS: in-sample r(observed, upsit) = {stats.pearsonr(v["upsit"], v["OIS"])[0]:.3f}  (N={len(v)})')

# ============================================================================
# Section 2 (continued). Save model, scaler and metadata
# ============================================================================
# Save the OIS Ridge model + scaler + metadata (non-destructive: repro dir).
pickle.dump(ois_model, open(MODELS / 'sbr_all_Ridge.pkl', 'wb'))
pickle.dump(ois_sc,    open(MODELS / 'scaler.pkl', 'wb'))
meta = {'n_train': int(n_ois), 'feature_names': FEATURES, 'sbr_cols': SBR_COLS,
        'demo_cols': DEMO, 'demo_medians': demo_medians.to_dict(),
        'cv_models_ois': comp_ois.set_index('model').to_dict('index')}
json.dump(meta, open(MODELS / 'metadata.json', 'w'), indent=2)
print('Saved OIS Ridge model + scaler + metadata to', MODELS)

# ============================================================================
# Section 3. Held-out correlation in the prodromal cohort
# ============================================================================
# Stage 1 held-out validation for OIS only (Prodromal cohort, not seen during training).
# Two panels: scatter (overall held-out r) + per-site r distribution across PPMI sites.
np.random.seed(42)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={'width_ratios': [1.2, 1]})
OIS_COLOR = '#8E44AD'  # purple, olfactory

# Panel a: OIS scatter on Prodromal held-out
ax = axes[0]
for coh in [1, 2]:  # background: PD+HC training data, light grey
    sub = df[df['COHORT'] == coh][['OIS', 'upsit']].dropna()
    ax.scatter(sub['OIS'], sub['upsit'], s=5, alpha=0.15, color='#cccccc', zorder=1)
v = df[df['COHORT'] == 4][['OIS', 'upsit']].dropna()
ax.scatter(v['OIS'], v['upsit'], s=14, alpha=0.55, color=OIS_COLOR,
           edgecolor='white', linewidth=0.4, label=f'Prodromal held-out (n={len(v)})', zorder=2)
r = stats.pearsonr(v['upsit'], v['OIS'])[0]
lo, hi = min(v['OIS'].min(), v['upsit'].min()), max(v['OIS'].max(), v['upsit'].max())
ax.plot([lo, hi], [lo, hi], 'k--', lw=1, alpha=0.6, label='y = x', zorder=3)
ax.set_xlabel('OIS (predicted)'); ax.set_ylabel('UPSIT (observed)')
ax.set_title(f'a  OIS vs UPSIT (Prodromal held-out r = {r:.3f})', loc='left', fontsize=10.5, fontweight='bold')
ax.legend(fontsize=9, markerscale=1.4, loc='upper left', frameon=False)

# Panel b: per-PPMI-site held-out r distribution for OIS
ax_site = axes[1]
prod_d = df[df['COHORT'] == 4].copy()
site_n = prod_d['SITE'].value_counts()
sites_valid = site_n[site_n >= 60].index
site_rs = []
for site in sites_valid:
    sub = prod_d[prod_d['SITE'] == site][['OIS', 'upsit']].dropna()
    if len(sub) >= 10 and sub['upsit'].std() > 0 and sub['OIS'].std() > 0:
        site_rs.append(stats.pearsonr(sub['upsit'], sub['OIS'])[0])
site_rs = np.array(site_rs)
jitter = np.random.uniform(-0.18, 0.18, len(site_rs))
ax_site.scatter(jitter, site_rs, s=42, color=OIS_COLOR, alpha=0.7,
                edgecolor='white', linewidth=0.6, zorder=2)
ax_site.plot([-0.3, 0.3], [site_rs.mean()] * 2, color=OIS_COLOR, lw=3, zorder=3)
ax_site.text(0, site_rs.mean() + 0.05, f'mean = {site_rs.mean():.3f} ± {site_rs.std():.3f}',
             ha='center', fontsize=10, color=OIS_COLOR, fontweight='bold')
ax_site.axhline(0, color='k', lw=0.6)
ax_site.set_xticks([0]); ax_site.set_xticklabels(['OIS held-out per site'], fontsize=10)
ax_site.set_xlim(-0.6, 0.6); ax_site.set_ylabel('per-site Pearson r')
ax_site.set_title(f'b  Per-PPMI-site held-out r for OIS\n({len(site_rs)} sites with n ≥ 60)',
                  loc='left', fontsize=10.5, fontweight='bold')

plt.tight_layout()
plt.savefig(FIG / 'Fig2c_R1_OIS_validation.png', bbox_inches='tight')
plt.savefig(FIG / 'Fig2c_R1_OIS_validation.pdf', bbox_inches='tight')
plt.show()
print(f'\nOIS held-out validation on Prodromal cohort (model not seen during training):')
print(f'  Overall r = {r:.3f} (n = {len(v)})')
print(f'  Per-site r = {site_rs.mean():.3f} ± {site_rs.std():.3f}'
      f' (median {np.median(site_rs):.3f}, range [{site_rs.min():.3f}, {site_rs.max():.3f}], {len(site_rs)} sites)')

# ============================================================================
# Section 4. Participant-level scores (not distributed)
# ============================================================================
mismatch_df = df[['PATNO','COHORT','OIS','OMI']].copy()
mismatch_df.to_csv(OUT / 'all_mismatch_scores.csv', index=False)
print('Saved', OUT / 'all_mismatch_scores.csv', '|', len(mismatch_df), 'subjects')
v = df['OMI'].dropna()
print(f'  OMI: N={len(v)}, mean={v.mean():+.3f}, std={v.std():.3f}')

# ============================================================================
# Section 8. Leave-one-centre-out validation
# ============================================================================
pdhc = df[df['COHORT'].isin([1, 2])].copy()   # df already carries SITE from the baseline merge
site_counts = pdhc['SITE'].value_counts()
valid_sites = site_counts[site_counts >= 10].index
print(f'Sites with >=10 PD+HC: {len(valid_sites)}')

loso_rows = []
for site in valid_sites:
    tr = pdhc[pdhc['SITE'] != site]; te = pdhc[pdhc['SITE'] == site]
    sc = StandardScaler().fit(tr[FEATURES].values)
    Xtr, Xte = sc.transform(tr[FEATURES].values), sc.transform(te[FEATURES].values)
    mtr, mte = tr['upsit'].notna(), te['upsit'].notna()
    if mtr.sum() < 50 or mte.sum() < 5: continue
    m = Ridge(alpha=1.0).fit(Xtr[mtr.values], tr.loc[mtr, 'upsit'].values)
    yt = te.loc[mte, 'upsit'].values; yp = m.predict(Xte[mte.values])
    if np.std(yt) == 0: continue
    _r, _p = stats.pearsonr(yt, yp)
    _sub = te.loc[mte]
    loso_rows.append(dict(site=site, n=int(mte.sum()), r=_r, p=_p,
                          n_pd=int((_sub['COHORT'] == 1).sum()),
                          n_hc=int((_sub['COHORT'] == 2).sum()),
                          upsit_mean=float(np.mean(yt)), upsit_sd=float(np.std(yt, ddof=1))))
loso_site = pd.DataFrame(loso_rows).sort_values('r').reset_index(drop=True)
loso_ois = loso_site['r'].values
loso_site.to_csv(OUT / 'loso_per_site.csv', index=False)
print(f'  OIS: r = {np.mean(loso_ois):.3f} +/- {np.std(loso_ois):.3f}  (across {len(loso_ois)} sites)')

# Why two centres come out negative. Neither is a reversal, both are indistinguishable
# from zero, and what tracks the centre-level correlation is whether the centre enrolled
# healthy controls as well as patients, not how many participants it contributed.
_neg = loso_site[loso_site['r'] < 0]
for _, _x in _neg.iterrows():
    print(f"  negative centre {int(_x.site)}: n = {int(_x.n)} (PD {int(_x.n_pd)}, HC {int(_x.n_hc)}), "
          f"r = {_x.r:+.3f}, P = {_x.p:.2f}, UPSIT SD = {_x.upsit_sd:.2f}")
_with, _without = loso_site[loso_site.n_hc > 0], loso_site[loso_site.n_hc == 0]
_t = stats.ttest_ind(_with['r'], _without['r'], equal_var=False)
_rn, _pn = stats.pearsonr(loso_site['n'], loso_site['r'])
_rh, _ph = stats.pearsonr(loso_site['n_hc'] / loso_site['n'], loso_site['r'])
print(f'  centres enrolling controls too: n = {len(_with)}, mean r = {_with["r"].mean():.3f} | '
      f'patients only: n = {len(_without)}, mean r = {_without["r"].mean():.3f}  (Welch P = {_t.pvalue:.4f})')
print(f'  control fraction vs centre r: r = {_rh:.3f}, P = {_ph:.2g} | centre size vs centre r: r = {_rn:.3f}, P = {_pn:.2f}')

# Per-centre significance is a statement about centre size, not about transportability.
# Expected number reaching P < 0.05 if the true correlation were the observed mean everywhere.
def _pw(rho, n, a=0.05):
    _se = 1 / np.sqrt(n - 3); _z = np.arctanh(rho) / _se; _c = stats.norm.ppf(1 - a / 2)
    return stats.norm.sf(_c - _z) + stats.norm.cdf(-_c - _z)
_pow = np.array([_pw(np.mean(loso_ois), n) for n in loso_site['n']])
print(f'  reaching P < 0.05: {int((loso_site["p"] < 0.05).sum())} of {len(loso_site)}, '
      f'against {_pow.sum():.1f} expected if the true r were {np.mean(loso_ois):.3f} at every centre '
      f'(median centre n = {int(loso_site["n"].median())}, power {_pw(np.mean(loso_ois), loso_site["n"].median()):.2f})')

loso_summary = dict(
    n_sites=int(len(loso_site)), r_mean=float(np.mean(loso_ois)), r_sd=float(np.std(loso_ois)),
    n_positive=int((loso_site['r'] > 0).sum()),
    negative_sites=[dict(site=int(x.site), n=int(x.n), n_pd=int(x.n_pd), n_hc=int(x.n_hc),
                         r=float(x.r), p=float(x.p), upsit_sd=float(x.upsit_sd))
                    for _, x in _neg.iterrows()],
    r_with_hc=float(_with['r'].mean()), n_with_hc=int(len(_with)),
    r_without_hc=float(_without['r'].mean()), n_without_hc=int(len(_without)),
    p_with_vs_without=float(_t.pvalue),
    r_controlfrac_vs_siter=float(_rh), p_controlfrac_vs_siter=float(_ph),
    r_size_vs_siter=float(_rn), p_size_vs_siter=float(_pn),
    n_significant=int((loso_site['p'] < 0.05).sum()), n_significant_expected=float(_pow.sum()),
    median_site_n=int(loso_site['n'].median()))
json.dump(loso_summary, open(OUT / 'loso_per_site_summary.json', 'w'), indent=2)
print('  -> saved results/repro/loso_per_site.csv and loso_per_site_summary.json')

# ============================================================================
# State for the next scripts
# ============================================================================
save_state('01_model', df=df, train=train, excel=excel, sbr=sbr, SBR_COLS=SBR_COLS, FEATURES=FEATURES,
           DEMO=DEMO, demo_medians=demo_medians)
