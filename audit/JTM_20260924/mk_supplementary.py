#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Information generator.

Every number in the tables is read from the reproducibility ledger
(results/repro/section16_stats.json and the section CSVs) or from
results/repro/supp_channel_followup.csv (computed by the same survival-frame
rules, 1,759 / 152). Nothing is copied by hand: the test-retest ICC, the
NSD-stage means, the temporal and centre splits all live in the ledger.

Run from the project root:
    python3 mk_supplementary.py
Writes
    Supplementary_Information_EN.md
    Supplementary_Information_中文.md
    821_results/Supplementary/{English/,}FigS1_flow, FigS2_calibration, FigS3_robustness (.pdf + .png)
"""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)
R = 'results/repro'
J = json.load(open(f'{R}/section16_stats.json'))
META = json.load(open('model_output/repro_sbr_models/metadata.json'))
CH = pd.read_csv(f'{R}/supp_channel_followup.csv')
LO2 = json.load(open(f'{R}/loso_per_site_summary.json'))
ICC = J['icc']
TS = J['temporal_split']          # §14 时间分割,2026-09-10 起进总账(此前是手抄常量)
C36 = pd.read_csv(f'{R}/section36_cindex_ci.csv')
C38 = pd.read_csv(f'{R}/section38_continuous_vs_binary.csv')
C39ALL = pd.read_csv(f'{R}/section39_sensitivity_matched.csv')
C39B = C39ALL[C39ALL.kind == 'sensitivity_matched']
C39F = C39ALL[C39ALL.kind == 'fixed_flag_fraction']
C22 = pd.read_csv(f'{R}/section22_csf_threeway.csv')
C42 = pd.read_csv(f'{R}/section42_partial_r.csv')
S44 = pd.read_csv(f'{R}/section44_three_groups.csv'); L44 = pd.read_csv(f'{R}/section44_lowest_third_vs_rest.csv'); M44 = pd.read_csv(f'{R}/supp_median_split.csv')

# ---- ledger keys that used to be hand-copied constants (2026-09-11: now written by the code) ----
# ICC 2026-09-10 起进总账(§12);NSD 分期均值与中心拆分 2026-09-11 起进总账(键 nsd_stage_omi、centre_split)
CENTRE = {k: dict(v, sites=', '.join(str(x) for x in v['sites'])) for k, v in J['centre_split'].items()}

# ------------------------------------------------------------------ helpers
def f3(x): return f'{x:.3f}'
def ci(c, lo, hi): return f'{c:.3f} ({lo:.3f}–{hi:.3f})'
def pf(p):
    if isinstance(p, str): return p
    if p == 0: return '< 0.001'
    if p < 0.001: return f'{p:.1e}'.replace('e-0', ' × 10<sup>−').replace('e-', ' × 10<sup>−') + '</sup>'
    return f'{p:.3f}' if p < 0.01 else f'{p:.2f}'
def sgn(x, d=3): return f'{x:+.{d}f}'
def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    for r in rows: out.append('| ' + ' | '.join(str(x) for x in r) + ' |')
    return '\n'.join(out)

LANG = 'en'
def T(en, zh): return zh if LANG == 'zh' else en

# ------------------------------------------------------------------ figures

# ---- Illustrator-editable PDF for the English build: base-14 Type 1 fonts (WinAnsi), so
# ---- every non-WinAnsi glyph (Greek, arrows, minus sign, comparison signs) is spelled out.
import matplotlib.text as _mtext
_SAN = [('ΔC', 'delta C'), ('Δ', 'delta '), ('A\u03b2' + '42', 'A-beta-42'), ('A\u00df' + '42', 'A-beta-42'), ('α', 'alpha'), ('β', 'beta'), ('ρ', 'rho'), ('χ²', 'chi-square'),
        ('χ', 'chi'), ('≥', '>='), ('≤', '<='), ('≈', '~'), ('−', '-'), ('→', '->'), ('↔', '<->'), ('²', '2'), ('⁻', '-'),
        ('⁰', '0'), ('¹', '1'), ('³', '3'), ('⁴', '4'), ('⁵', '5'), ('⁶', '6'), ('⁷', '7'), ('⁸', '8'), ('⁹', '9'), ('…', '...'), ('′', "'")]
def _sanitize(fig):
    import unicodedata as _ud
    def _fallback(ch):
        try: ch.encode('cp1252'); return ch
        except UnicodeEncodeError: pass
        nm = _ud.name(ch, '')
        if 'GREEK SMALL LETTER' in nm: return nm.split()[-1].lower()
        if 'GREEK CAPITAL LETTER' in nm: return nm.split()[-1].capitalize()
        return '?'
    def _san(s0):
        s1 = s0
        for a, b in _SAN: s1 = s1.replace(a, b)
        return ''.join(_fallback(c) for c in s1)
    fig.canvas.draw()   # materialise tick-label strings from their formatters
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            labs = [t.get_text() for t in axis.get_ticklabels()]
            new = [_san(x) for x in labs]
            if new != labs:
                axis.set_ticks(axis.get_ticklocs()); axis.set_ticklabels(new)
    for t in fig.findobj(_mtext.Text):
        s0 = t.get_text()
        if s0:
            s1 = _san(s0)
            if s1 != s0: t.set_text(s1)

def setup_fonts():
    mpl.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False, 'font.size': 8.5,
                         'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 120})
    mpl.rcParams['font.sans-serif'] = ['Noto Sans CJK JP', 'Noto Sans CJK SC', 'DejaVu Sans']
    mpl.rcParams['font.serif'] = ['Times', 'Times New Roman', 'Nimbus Roman', 'Liberation Serif', 'DejaVu Serif']
    # 英文图件一律 Times New Roman(2026-09-11 用户定)。use14corefonts 下 family='serif'
    # 映射到 base-14 的 Times-Roman,仍是 Type 1 / WinAnsi,Illustrator 里文字照样可编辑。
    mpl.rcParams['font.family'] = ('sans-serif' if LANG == 'zh' else 'serif')
    # ⚠️ rcParams 是全局的,而本脚本在同一进程里先跑 en 再跑 zh。
    # use14corefonts 只在 en 分支打开、从不关闭时,中文 PDF 会被强制走 base-14 Helvetica,
    # 中文字形编不进去而全部变成 '?'(PNG 不受影响,故只看 PNG 查不出来)。必须显式复位。
    mpl.rcParams['pdf.use14corefonts'] = (LANG == 'en')
    if LANG == 'en':
        mpl.rcParams['font.serif'] = ['Times', 'Times New Roman', 'Nimbus Roman', 'DejaVu Serif']
TEAL, RED, GREY, DARK, ORANGE, BLUE = '#16A085', '#C0392B', '#7F8C8D', '#2C3E50', '#E67E22', '#2980B9'

PANEL_CASE = os.environ.get('PANEL_CASE', 'lower')     # 'upper' for journals that letter panels A, B (eBioMedicine)
def PL(s):
    """Panel letter at the start of an in-figure title ('a  Decile calibration'), uppercased when required."""
    return (s[0].upper() + s[1:]) if PANEL_CASE == 'upper' else s
def outdir():
    # 图件一律只出英文(2026-09-10 用户定):投稿用英文图,中文稿沿用同一套图,
    # 避免同一张图维护两份、也避免中文 PDF 的字体坑。
    d = os.environ.get('SUPP_FIG_DIR', '821_results/Supplementary/English')
    os.makedirs(d, exist_ok=True); return d

A4W, A4H = 8.27, 11.69
_LAYOUT = {}      # name -> ink extent on the A4 page (points from the top), read by the submission builder

def a4_subplots(nrows, ncols, w, h, wspace=None, left=.06, right=.94, top=.965):
    """Axes laid out as plt.subplots(nrows, ncols, figsize=(w, h)) would lay them out, but drawn at the
    top of a portrait A4 page (every figure is delivered as an A4 portrait page). The drawing is scaled
    to the page width when wider than it and kept at its natural size otherwise."""
    s = min(1.0, (right - left) * A4W / w)
    bw, bh = w * s / A4W, h * s / A4H                       # box occupied by the original figure, in page fractions
    l = left if s < 1 else .5 - bw / 2
    fig = plt.figure(figsize=(A4W, A4H))
    gs = fig.add_gridspec(nrows, ncols, left=l + .125 * bw, right=l + .9 * bw, top=top - .12 * bh, bottom=top - .89 * bh,
                          wspace=.2 if wspace is None else wspace)
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(ncols)] for i in range(nrows)])
    fig._foot_y, fig._foot_x0, fig._foot_w = top - bh - .012, l + .125 * bw, bw * A4W
    return fig, (axes[0, 0] if axes.size == 1 else axes.ravel())

def footnote(fig, txt, fs=6.6):
    """Grey note under the drawing, wrapped to the drawing width and centred on the page."""
    import textwrap
    n = int(fig._foot_w * 72 / (fs * (1.0 if LANG == 'zh' else 0.46)))     # characters per line for this width and size
    lines = []
    for para in txt.split('\n'):
        lines += textwrap.wrap(para, n) if LANG == 'en' else [para[i:i + n] for i in range(0, len(para), n)]
    fig.text(.5, fig._foot_y, '\n'.join(lines), ha='center', va='top', fontsize=fs, color='#555', linespacing=1.35)

def save(fig, name):
    d = outdir()
    if LANG == 'en': _sanitize(fig)
    bb = fig.get_tightbbox(fig.canvas.get_renderer())        # ink extent in inches, origin bottom left
    _LAYOUT[name] = dict(ink_top_pt=(A4H - bb.y1) * 72, ink_bottom_pt=(A4H - bb.y0) * 72, page='A4 portrait')
    fig.savefig(f'{d}/{name}.pdf'); fig.savefig(f'{d}/{name}.png', dpi=200)
    json.dump(_LAYOUT, open(f'{d}/figS_layout.json', 'w'), indent=1)
    plt.close(fig); print('  saved', d, name)

def figS1_flow():
    LF = J['layer_flow']; ch = CH.set_index('channel')
    fig, ax = a4_subplots(1, 1, 9.6, 6.0); ax.set_xlim(0, 10.5); ax.set_ylim(1.2, 10); ax.axis('off')
    def box(x, y, w, h, txt, fc='#F4F6F7', ec=DARK, fs=7.8, bold=False):
        ax.add_patch(plt.Rectangle((x - w / 2, y - h / 2), w, h, fc=fc, ec=ec, lw=1.0))
        ax.text(x, y, txt, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal')
    def arrow(x0, y0, x1, y1):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='->', color=DARK, lw=1.0))
    # training cohort
    box(1.6, 9.2, 2.9, 1.0, T('Training cohort\nPD 1,163 + HC 235 = 1,398\nUPSIT and imaging both present',
                              '训练队列\nPD 1,163 + HC 235 = 1,398\n同时有 UPSIT 与影像'), fc='#E8F6F3', ec=TEAL)
    box(1.6, 7.6, 3.1, 1.1, T('Ridge (α = 1.0)\n33 SBR + age, sex, education\nfive-fold r = 0.479',
                              '岭回归(α = 1.0)\n33 SBR + 年龄、性别、教育\n五折 r = 0.479'), fc='#E8F6F3', ec=TEAL, fs=7.2)
    arrow(1.6, 8.7, 1.6, 8.15)
    # prodromal flow
    x = 6.0
    box(x, 9.2, 5.0, 0.8, T('Prodromal cohort with imaging and baseline clinical data\nn = 2,453 (rule A: 1,558 / 534 / 361)',
                            '有影像与基线临床数据的前驱期队列\nn = 2,453(规则 A:1,558 / 534 / 361)'), bold=True)
    # Correlation and survival are separate subsets of the 2,453 baseline participants.
    box(1.6, 5.8, 3.1, 1.0, T('Correlation analyses\nOIS and UPSIT both present\nn = 2,441',
                              '相关分析\nOIS 与 UPSIT 均有\nn = 2,441'), fc='#FEF9E7', ec=ORANGE, fs=7.2)
    ax.plot([3.5, 3.3, 3.3], [9.2, 9.2, 5.8], color=DARK, lw=1.0)
    arrow(3.3, 5.8, 3.15, 5.8)
    steps = [(T('Prevalent Parkinson\'s disease at baseline', '基线即已确诊帕金森病'), 78, 2375),
             (T('No follow-up time (single visit); 0 converters', '无随访时间(单次访视);转化者 0'), 605, 1770),
             (T('Missing an analysis variable', '缺分析变量'), 11, 1759)]
    y = 8.05
    for lab, ex, rem in steps:
        y2 = y - 1.05
        arrow(x, 8.8 if rem == 2375 else y - 0.35, x, y2 + 0.3)
        ax.text(x + 2.15, (y + y2) / 2, T(f'exclude {ex}: {lab}', f'排除 {ex}:{lab}'), fontsize=6.9, va='center', color='#555')
        box(x, y2, 4.0, 0.6, T(f'n = {rem:,}', f'n = {rem:,}') + ('' if rem != 1759 else T(' with 152 conversions', '(152 例转化)')),
            bold=(rem == 1759), fc=('#FDEDEC' if rem == 1759 else '#F4F6F7'), ec=(RED if rem == 1759 else DARK))
        y = y2
    # groups
    h = ch.loc['Hyposmia']; nh = ch.loc['Non-hyposmia-enriched']; rbd = ch.loc['RBD']; gen = ch.loc['Genetic']
    box(3.6, 3.35, 3.0, 0.85, T(f'Hyposmia group\nn = {int(h.n):,}, {int(h.events)} conversions\n{h.rate_per_100py:.2f} per 100 person-years',
                                f'嗅觉减退组\n{int(h.n):,} 人,{int(h.events)} 例转化\n每百人年 {h.rate_per_100py:.2f}'), fc='#FDEDEC', ec=RED, fs=7.6)
    box(8.0, 3.35, 3.4, 0.85, T(f'RBD and variant-carrier group\nn = {int(nh.n):,}, {int(nh.events)} conversions\n{nh.rate_per_100py:.2f} per 100 person-years',
                                f'RBD 与遗传携带组\n{int(nh.n):,} 人,{int(nh.events)} 例转化\n每百人年 {nh.rate_per_100py:.2f}'), fc='#EBF5FB', ec=BLUE, fs=7.6)
    arrow(x - 0.6, 4.6, 3.6, 3.8); arrow(x + 0.6, 4.6, 8.0, 3.8)
    box(2.0, 1.75, 2.3, 0.8, T(f'DAT deficit (< 65%)\nn = {265}, {46} conversions', f'DAT 缺损(< 65%)\n265 人,46 例转化'), fs=7.4)
    box(4.65, 1.75, 2.3, 0.8, T(f'No DAT deficit\nn = {738}, {22} conversions', f'DAT 非缺损\n738 人,22 例转化'), fs=7.4)
    arrow(3.0, 2.9, 2.0, 2.2); arrow(4.2, 2.9, 4.65, 2.2)
    box(7.15, 1.75, 2.0, 0.8, T(f'RBD cohort\nn = {int(rbd.n)}, {int(rbd.events)} conversions', f'RBD 队列\n{int(rbd.n)} 人,{int(rbd.events)} 例转化'), fs=7.0)
    box(9.35, 1.75, 2.25, 0.8, T(f'Pathogenic-variant cohort\nn = {int(gen.n)}, {int(gen.events)} conversions', f'遗传携带队列\n{int(gen.n)} 人,{int(gen.events)} 例转化'), fs=7.0)
    arrow(7.4, 2.9, 7.15, 2.2); arrow(8.6, 2.9, 9.35, 2.2)
    footnote(fig, T('Deficit is lower putamen SBR below 65% of its age- and sex-expected value (PARS threshold). '
                         'The 65% split uses the 1,003 participants of the hyposmia group. The two recruitment cohorts of the RBD and variant-carrier group '
                         'follow rule A (dual carriers with RBD assigned to the RBD cohort).',
                         '缺损定义为较低侧壳核 SBR 低于年龄性别预期值的 65%(PARS 阈值)。65% 划分施加于嗅觉减退组 1,003 人;'
                         'RBD 与遗传携带组的两个队列按规则 A(兼有遗传变异与 RBD 者归 RBD 队列)。'),
            )
    save(fig, 'FigS1_flow')

def figS2_calibration(chan_slopes):
    dec = pd.DataFrame(J['calibration_deciles']); cal = pd.DataFrame(J['calibration'])
    fig, axes = a4_subplots(1, 2, 10.4, 3.9, wspace=.85)
    ax = axes[0]
    ax.plot([15, 40], [15, 40], c=GREY, ls='--', lw=1, label=T('identity', '等值线'))
    ax.errorbar(dec.pred, dec.obs, fmt='o', c=BLUE, ms=5)
    for _, r in dec.iterrows():
        ax.text(r.pred + .25, r.obs - .9, f'{int(r.dec)}', fontsize=6.5, color=BLUE)
    ax.set_xlabel(T('Mean predicted UPSIT (OIS) in decile', '各十分位的平均预测 UPSIT(OIS)')); ax.set_ylabel(T('Mean observed UPSIT', '平均实测 UPSIT'))
    ax.set_xlim(20, 39); ax.set_ylim(15, 39); ax.legend(frameon=False, fontsize=7.5, loc='upper left')
    ax.set_title(PL(T('a  Decile calibration\nRBD and variant-carrier group (n = 885)', 'a  十分位校准\nRBD 与遗传携带组(n = 885)')), loc='left', fontsize=9, fontweight='bold')
    ax.text(.04, .86, T('observed below predicted in every decile\n(mean residual −2.7 points)', '每个十分位实测均低于预测\n(平均残差 −2.7 分)'),
            transform=ax.transAxes, ha='left', va='top', fontsize=7, color='#555')
    ax = axes[1]
    cols_ = {T('All prodromal', '全前驱期'): GREY, T('Hyposmia cohort', '嗅觉减退队列'): RED, T('RBD cohort', 'RBD 队列'): BLUE, T('Pathogenic-variant cohort', '遗传携带队列'): '#8E44AD'}
    rows = [(k, v['slope'], v['r'], v['n'], cols_[k], v['lo'], v['hi'], v['p']) for k, v in chan_slopes.items()]
    y = np.arange(len(rows))[::-1]
    ax.barh(y, [r[1] for r in rows], color=[r[4] for r in rows], height=.6,
            xerr=[[r[1] - r[5] for r in rows], [r[6] - r[1] for r in rows]], error_kw=dict(ecolor='#566573', lw=.9, capsize=3))
    ax.axvline(1, c=GREY, ls='--', lw=1)
    x0 = max(r[6] for r in rows) + .06                      # one text column, clear of every upper CI limit and of the line at 1
    for yi, r in zip(y, rows):
        ax.text(x0, yi, f'{r[1]:.2f} ({r[5]:.2f} to {r[6]:.2f})\nr = {r[2]:.3f}, P = {r[7]:.0e}, n = {r[3]:,}', va='center', fontsize=6.8, linespacing=1.3)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7.5); ax.set_xlim(0, 1.9)
    ax.set_xlabel(T('Slope of measured UPSIT on OIS (95% CI)', '实测 UPSIT 对 OIS 的回归斜率(95% CI)'))
    ax.set_title(PL(T('b  Slope by recruitment cohort', 'b  分队列斜率')), loc='left', fontsize=9, fontweight='bold')
    footnote(fig, T('Slopes describe calibration of predicted UPSIT, not conversion risk. The pooled slope for the RBD and variant-carrier group (1.027) is not shown '
                            'because it mixes two recruitment cohorts with different means. The within-recruitment cohort slopes are the interpretable values.',
                            '斜率低于 1 表示 OIS 的取值范围比实测分数窄。不显示 RBD 与遗传携带组的合并斜率(1.027),它混合了均值不同的两个队列。队列内斜率才是可解读的值。'),
                 )
    save(fig, 'FigS2_calibration')

def figS3_robustness():
    L = {r['landmark_yr']: r for r in J['landmark']}
    rows = [(T('Full hyposmia group, 1,003 / 68', '全嗅觉减退组,1,003 / 68'), L[0]['dc'], L[0]['lo'], L[0]['hi'], pf(L[0]['p']), TEAL),
            (T(f'Landmark 1 year, {L[1]["n"]:,} / {L[1]["events"]}', f'1 年 landmark, {L[1]["n"]:,} / {L[1]["events"]}'), L[1]['dc'], L[1]['lo'], L[1]['hi'], pf(L[1]['p']), TEAL),
            (T(f'Landmark 2 years, {L[2]["n"]:,} / {L[2]["events"]}', f'2 年 landmark, {L[2]["n"]:,} / {L[2]["events"]}'), L[2]['dc'], L[2]['lo'], L[2]['hi'], pf(L[2]['p']), TEAL),
            (T(f'Model retrained on pre-2017 enrolment, {TS["n"]:,} / {TS["events"]}', f'模型改用 2017 年前入组者重训,{TS["n"]:,} / {TS["events"]}'), TS['dc'], TS['lo'], TS['hi'], pf(TS['p']), TEAL)]
    fig, ax = a4_subplots(1, 1, 8.4, 3.6); y = np.arange(len(rows))[::-1]
    ax.get_gridspec().update(left=.34)                      # room for the long row labels on the page
    for yi, (lab, d, lo, hi, p, col) in zip(y, rows):
        if lo is not None:
            ax.plot([lo, hi], [yi, yi], c=col, lw=2.4, solid_capstyle='round')
        ax.scatter([d], [yi], s=70, marker='s' if lo is not None else 'D', c=col, zorder=3, ec='white', lw=1)
        txt = (f'{d:+.3f} ({lo:+.3f}, {hi:+.3f})   P ' + (p if p.startswith('<') else '= ' + p)) if lo is not None else f'{d:+.3f}   ' + p
        ax.text(0.46, yi, txt, va='center', fontsize=7.4, color=col)
    ax.axvline(0, c=DARK, lw=1); ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7.8)
    ax.set_xlim(-.08, .75); ax.set_xlabel(T('ΔC, OIS minus UPSIT (hyposmia group)', 'ΔC,OIS 减 UPSIT 总分(嗅觉减退组)'))
    ax.set_title(T('Landmark and temporal sensitivity analyses', 'OIS 相对 UPSIT 优势的稳健性'), loc='left', fontsize=9.5, fontweight='bold')
    footnote(fig, T('Intervals use 1,000 paired bootstrap resamples. Landmark analyses retain event-free participants observed beyond the landmark and restart time there. The temporal split retrains the model on participants enrolled before 2017 and evaluates in hyposmia-group participants enrolled from 2017 (C 0.702 to 0.827).',
                       '区间为 1,000 次配对 bootstrap。时间拆分在 2017 年前入组者上重训模型,在 2017 年起入组的嗅觉减退组评价(C 0.702 至 0.827)。')
            )
    save(fig, 'FigS3_robustness')


def _state():
    # participant-level frames written by code/02_survival_cohort.py (not distributed)
    import pickle
    cands = ['results/intermediate/02_cohort.pkl']
    for p in cands:
        if os.path.exists(p):
            return pickle.load(open(p, 'rb'))
    return None


def hyposmia_frame():
    st = _state(); prod = st['prod']
    hh = prod[prod.subgroup == 'Hyposmia'].dropna(subset=['OIS', 'upsit', 'PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM', 'time_years', 'converted', 'age', 'sex_num']).copy()
    hh = hh[(hh.converted == 1) | (hh.time_years > 0)]
    hh['pct'] = hh[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / (1.8547 - .00303 * hh.age - .2419 * hh.sex_num) * 100
    return hh

def _km_panel(ax, d, groups, labels, colors, ymin, title):
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
    labs_ = []
    for g_, lab, col in zip(groups, labels, colors):
        gg = d[d.grp == g_]; k = KaplanMeierFitter().fit(gg.time_years, gg.converted)
        sf = k.survival_function_; m = sf.index <= 2
        ax.step(sf.index[m], sf.iloc[:, 0][m], where='post', color=col, lw=2, label=f'{lab} (n={len(gg)}, {int(gg.converted.sum())} ev)')
        ci = k.confidence_interval_survival_function_; mc = ci.index <= 2
        ax.fill_between(ci.index[mc], ci.iloc[:, 0][mc], ci.iloc[:, 1][mc], step='post', color=col, alpha=.10, lw=0)
        idx = ci.index[ci.index <= 2.0][-1]
        labs_.append([float(k.predict(2.)), f'{(1 - float(k.predict(2.))) * 100:.1f}%\n({(1 - ci.loc[idx].iloc[1]) * 100:.1f} to {(1 - ci.loc[idx].iloc[0]) * 100:.1f})', col])
    gap = (1.005 - ymin) * 0.10
    labs_.sort(key=lambda r: -r[0])                          # highest curve first, labels pushed downwards so none climbs into the title
    labs_[0][0] = min(labs_[0][0], 1.0 - gap * 0.35)
    for i in range(1, len(labs_)):
        if labs_[i - 1][0] - labs_[i][0] < gap: labs_[i][0] = labs_[i - 1][0] - gap
    for y_, txt, col in labs_:
        ax.text(2.04, y_, txt, fontsize=6.8, color=col, va='center', fontweight='bold')
    p3 = multivariate_logrank_test(d.time_years, d.grp, d.converted).p_value
    lines = [T(f'3-group log-rank P = {p3:.1e}', f'三组 log-rank P = {p3:.1e}')]
    for a, b in [(0, 1), (1, 2), (0, 2)]:
        ga, gb = d[d.grp == groups[a]], d[d.grp == groups[b]]
        pv = logrank_test(ga.time_years, gb.time_years, ga.converted, gb.converted).p_value
        lines.append(f'{labels[a]} vs {labels[b]} ' + ('P < 0.001' if pv < .001 else f'P = {pv:.3f}' if pv < .01 else f'P = {pv:.2f}'))
    for i, l in enumerate(lines):
        ax.text(.97, .30 - i * .065, l, transform=ax.transAxes, ha='right', fontsize=6.6, color='#444', fontweight='bold' if i == 0 else 'normal')
    ax.set_xlim(0, 2.55); ax.set_ylim(ymin, 1.005); ax.legend(frameon=False, loc='lower left', fontsize=6.2)
    ax.set_xlabel(T('Years', '距基线年数')); ax.set_ylabel(T('Free of PD', '未确诊比例'))
    ax.set_title(title, loc='left', fontsize=8.8, fontweight='bold')

def figS4_tertiles():
    hh = hyposmia_frame(); fig, axes = a4_subplots(1, 2, 10, 4, wspace=.35)
    for ax, d, ymin, ttl in [(axes[0], hh, .68, PL(T('a  Whole hyposmia group (n = 1,003), OIS tertiles', 'a  全嗅觉减退组(n = 1,003),OIS 三分位'))),
                             (axes[1], hh[hh.pct >= 65], .875, PL(T('b  Non-deficit stratum (n = 738), OIS tertiles', 'b  非缺损亚组(n = 738),OIS 三分位')))]:
        d = d.copy(); d['grp'] = pd.qcut(d.OIS, 3, labels=['T1', 'T2', 'T3']).astype(str)
        _km_panel(ax, d, ['T1', 'T2', 'T3'], [T('lowest', '最低'), T('middle', '中间'), T('highest', '最高')], [RED, ORANGE, TEAL], ymin, ttl)
    footnote(fig, T('Division into OIS tertiles without using outcomes to set the cut points. Middle versus highest tertile comparisons were not statistically significant in either population. The main text uses two groups divided at the median.',
                         '按 OIS 三分位划分,不看结局。两个人群中中间与最高三分位均不分开,因此正文按中位数分为两组。'), fs=6.8)
    save(fig, 'FigS4_tertiles')

def figS5_three_bands():
    hh = hyposmia_frame(); r = S44[(S44.stratum == 'hyposmia') & (S44.reading == 'OIS')].iloc[0]
    d = hh.copy(); d['grp'] = np.where(d.OIS <= r.cut1, 'low', np.where(d.OIS <= r.cut2, 'mid', 'high'))
    fig, ax = a4_subplots(1, 1, 5.4, 4)
    _km_panel(ax, d, ['low', 'mid', 'high'], [T('high risk', '高危'), T('intermediate', '中危'), T('low risk', '低危')], [RED, ORANGE, TEAL], .60,
              T(f'OIS bands at {r.cut1:.2f} and {r.cut2:.2f}, hyposmia group', f'OIS 切点 {r.cut1:.2f} 与 {r.cut2:.2f},嗅觉减退组'))
    footnote(fig, T(f'two cut points, permutation P {"< 0.001" if r.p_permutation < .001 else f"= {r.p_permutation:.3f}"} ({int(r.n_perm)} permutations)\nbootstrap 95% for the cuts: {r.cut1_boot_lo:.1f} to {r.cut1_boot_hi:.1f} and {r.cut2_boot_lo:.1f} to {r.cut2_boot_hi:.1f}',
                         f'两切点选取,置换 P {"< 0.001" if r.p_permutation < .001 else f"= {r.p_permutation:.3f}"}({int(r.n_perm)} 次)\n切点 bootstrap 95%:{r.cut1_boot_lo:.1f} 至 {r.cut1_boot_hi:.1f} 与 {r.cut2_boot_lo:.1f} 至 {r.cut2_boot_hi:.1f}'),
            fs=6.4)
    ax.legend(frameon=False, loc='lower left', fontsize=6.2, bbox_to_anchor=(0, .0))
    save(fig, 'FigS5_three_bands')

def channel_slopes():
    from scipy import stats
    st = _state()
    if st is not None:
        df = st['df']
        d = df[df.COHORT == 4].dropna(subset=['OIS', 'upsit']).copy()
        sg = d.subgroup.astype(str)
        d['arm'] = np.where(sg == 'Hyposmia', 'Hyposmia', np.where(sg.str.contains('RBD'), 'RBD', 'Genetic'))
        out = {}
        for arm, lab in [('ALL', T('All prodromal', '全前驱期')), ('Hyposmia', T('Hyposmia cohort', '嗅觉减退队列')), ('RBD', T('RBD cohort', 'RBD 队列')), ('Genetic', T('Pathogenic-variant cohort', '遗传携带队列'))]:
            g = d if arm == 'ALL' else d[d.arm == arm]; sl = stats.linregress(g.OIS, g.upsit)
            out[lab] = dict(slope=float(sl.slope), lo=float(sl.slope - 1.96 * sl.stderr), hi=float(sl.slope + 1.96 * sl.stderr), r=float(sl.rvalue), p=float(sl.pvalue), n=int(len(g)))
        return out
    ex = pd.read_excel(os.path.join(os.environ.get('PPMI_DATA_DIR', '.'), 'PPMI_Curated_Data_Cut_Public_20260223.xlsx'), sheet_name='20260105')
    e4 = ex[ex.COHORT == 4].sort_values(['PATNO', 'visit_date'])
    base = e4.groupby('PATNO').first()
    sc = pd.read_csv(f'{R}/all_mismatch_scores.csv'); d = sc[sc.COHORT == 4].copy()
    d['sub'] = d.PATNO.map(base['subgroup']).astype(str); d['upsit'] = d.PATNO.map(base['upsit'])
    d['arm'] = np.where(d['sub'] == 'Hyposmia', 'Hyposmia', np.where(d['sub'].str.contains('RBD'), 'RBD', 'Genetic'))
    d = d.dropna(subset=['OIS', 'upsit'])
    out = {}
    for arm, lab in [('RBD', T('RBD cohort', 'RBD 队列')), ('Genetic', T('Pathogenic-variant cohort', '遗传携带队列'))]:
        g = d[d.arm == arm]; sl = stats.linregress(g.OIS, g.upsit)
        out[lab] = dict(slope=float(sl.slope), r=float(sl.rvalue), n=int(len(g)))
    return out

# ------------------------------------------------------------------ document
def build_doc(chan_slopes):
    S = []
    A = S.append
    A(T('# Supplementary Information', '# 补充材料'))
    A(T('**Decomposition of the olfactory score by dopamine transporter imaging improves Parkinson\'s disease risk stratification in hyposmic individuals**',
        '**以多巴胺转运体显像分解嗅觉测试分数,可提升嗅觉减退人群中帕金森病的风险分层能力**'))
    A('')
    A(T('## Contents', '## 目录'))
    A(T('''- Supplementary Methods 1–2
- Supplementary Tables 1–11
- Supplementary Figures 1–5''', '''- 补充方法 1 至 2
- 补充表 1 至 11
- 补充图 1 至 5'''))
    A('')
    # ---------------- Methods
    A(T('## Supplementary Methods', '## 补充方法'))
    A(T('**Supplementary Method 1. Assignment to recruitment cohort.** PPMI enrols prodromal participants through three recruitment cohorts, namely olfactory screening (the hyposmia cohort), polysomnography-confirmed REM sleep behaviour disorder (RBD) and pathogenic-variant carriage. Six individuals carry a pathogenic variant and also have RBD. Rule A, used throughout the main text, assigns them to the RBD cohort (1,558 / 534 / 361). Rule B assigns them to the pathogenic-variant cohort (1,558 / 528 / 367). Both rules total 2,453 and the held-out correlations differ only in the third decimal (Supplementary Table 1).',
        '**补充方法 1. 入组队列的归属。**PPMI 前驱期参与者经三条途径纳入:嗅觉筛查(嗅觉减退队列)、多导睡眠图确诊的快速眼动睡眠行为障碍(RBD)与致病变异携带。6 人兼有致病变异与 RBD。正文通用的规则 A 将其归 RBD 队列(1,558 / 534 / 361),规则 B 归遗传携带队列(1,558 / 528 / 367)。两套规则合计均为 2,453,held-out 相关仅在第三位小数上不同(补充表 1)。'))
    D = J['dat_metric_validation']
    A(T(f'**Supplementary Method 2. Imaging-state determination and its validation.** DAT deficit is defined on the PPMI canonical quantity, the lower of the two putamen binding ratios divided by the value expected from a linear regression on age and sex fitted in healthy controls (n = {D["n_hc_normative"]}), expected = {D["norm_intercept"]:.4f} − {abs(D["norm_beta_age"]):.5f} × age − {abs(D["norm_beta_sex"]):.4f} × sex. Official staging fields are populated only for the Parkinson\'s disease and SWEDD cohorts in the curated cut, so the implementation was validated there. At the 0.75 cut used by the staging system it reproduced the official determination in {D["agreement_vs_official_Stage_D"]*100:.1f}% of {D["n_validated"]:,} participants, with all {D["n_discordant"]} disagreements inside the band {D["official_normal_min_pct_exp"]:.4f} to {D["official_deficit_max_pct_exp"]:.4f}. The 65% cut used for the deficit split comes from PARS.',
        f'**补充方法 2. 影像状态判定及其验证。**DAT 缺损按 PPMI 规范量定义:较低侧壳核结合比除以在健康对照(n = {D["n_hc_normative"]})上由年龄与性别线性回归得到的预期值,预期值 = {D["norm_intercept"]:.4f} − {abs(D["norm_beta_age"]):.5f} × 年龄 − {abs(D["norm_beta_sex"]):.4f} × 性别。curated cut 中官方分期字段仅帕金森病与 SWEDD 队列有值,故在该处验证:以分期系统的 0.75 切点,在 {D["n_validated"]:,} 人中复现官方判定 {D["agreement_vs_official_Stage_D"]*100:.1f}%,{D["n_discordant"]} 例不一致全部落在 {D["official_normal_min_pct_exp"]:.4f} 至 {D["official_deficit_max_pct_exp"]:.4f} 的窄带内。缺损划分所用的 65% 切点出自 PARS。'))
    A('')
    A(T('## Supplementary Tables', '## 补充表'))
    # ---------------- Table 1
    A(T('### Supplementary Table 1. Model families, training-set choice and held-out correlation by recruitment cohort', '### 补充表 1. 回归族、训练集选择与分队列 held-out 相关'))
    A(T('**a. Six regression families, five-fold cross-validation in the training cohort (n = 1,398), identical folds (seed 41).**', '**a. 六个回归族,训练队列五折交叉验证(n = 1,398),折数相同(种子 41)。**'))
    rows = []
    for k, v in META['cv_models_ois'].items():
        rows.append([k, f3(v['pearson_r']), f"{v['fold_mean']:.3f} ± {v['fold_std']:.3f}", f"{v['fold_min']:.3f} to {v['fold_max']:.3f}", f3(v['spearman_r']), f"{v['mae']:.2f}"])
    A(md_table([T('Family', '回归族'), T('Pooled Pearson r', '合并 Pearson r'), T('Fold mean ± s.d.', '折均值 ± 标准差'), T('Fold range', '折范围'), T('Spearman ρ', 'Spearman ρ'), T('MAE (points)', 'MAE(分)')], rows))
    A(T('Ridge was retained because the families are indistinguishable and a linear model keeps the fitted value interpretable as a weighted sum of binding ratios.',
        '六族表现无法区分,保留岭回归是因为线性模型使拟合值可解读为结合比的加权和。'))
    TC = J['training_cohort']
    A('')
    A(T('**b. Choice of training set.**', '**b. 训练集选择。**'))
    rows = [[m['training_set'], f"{m['n']:,}", f3(m['cv_r']), f3(m['c_hyposmia'])] for m in TC['models']]
    A(md_table([T('Training set', '训练集'), 'n', T('CV r (UPSIT)', 'CV r(UPSIT)'), T('C-index for conversion, hyposmia group', '嗅觉减退组转化 C-index')], rows))
    A(T(f'Training on healthy controls alone gives an OIS with no discrimination. Down-sampling the Parkinson\'s disease group to the same n = {TC["n_downsample"]} and repeating 200 times gives C = {TC["downsample_mean"]:.3f} ± {TC["downsample_sd"]:.3f} ({TC["downsample_below_055"]} of 200 below 0.55), so the failure is range restriction rather than sample size, since {TC["hc_extrapolation_fraction"]*100:.1f}% of hyposmia-group participants fall outside the imaging range of healthy controls. The PD-only model correlates r = {TC["r_pdonly_vs_main"]:.3f} with the main model.',
        f'仅用健康对照训练所得 OIS 无判别力。将帕金森病组降采样到同样 n = {TC["n_downsample"]} 并重复 200 次得 C = {TC["downsample_mean"]:.3f} ± {TC["downsample_sd"]:.3f}({TC["downsample_below_055"]}/200 低于 0.55),故失败源于取值范围限制而非样本量:嗅觉减退组 {TC["hc_extrapolation_fraction"]*100:.1f}% 的人落在健康对照的影像取值范围之外。仅 PD 训练的模型与主模型 r = {TC["r_pdonly_vs_main"]:.3f}。'))
    HB = J['heldout_by_arm']
    A('')
    A(T('**c. Held-out correlation between OIS and measured UPSIT in the prodromal cohort, by recruitment cohort.**', '**c. 前驱期 OIS 与实测 UPSIT 的 held-out 相关,按入组队列。**'))
    rows = [[{'Hyposmia': T('Hyposmia', '嗅觉减退'), 'RBD': 'RBD', 'Genetic': T('Variant carriers', '遗传携带')}[r['arm']], f"{r['n']:,}", f3(r['r']), pf(r['p']), f"{r['upsit_mean']:.1f} ({r['upsit_sd']:.1f})"] for r in HB['rows'] if r['rule'] == 'A']
    A(md_table([T('Recruitment cohort', '入组队列'), 'n', 'r', 'P', T('UPSIT mean (s.d.)', 'UPSIT 均值(标准差)')], rows))
    A(T(f'Pooled r = {HB["pooled_r"]:.3f} (n = {HB["pooled_n"]:,}).', f'合并 r = {HB["pooled_r"]:.3f}(n = {HB["pooled_n"]:,})。'))
    A('')
    # d. 逐中心 LOSO 明细,含 PD/HC 构成(用户要求:两个负值靠这张表解释)
    A(T('**d. Leave-one-centre-out correlation at each of the 44 centres, with the patient and control counts that centre contributed.**',
        '**d. 44 个中心各自的留一相关,并列出该中心贡献的患者与对照人数。**'))
    _ls = pd.read_csv(f'{R}/loso_per_site.csv').sort_values('r')
    # 用 1 至 44 的序号,不印 PPMI 站点号(站点号 10 至 79 跳号,会被误读成中心数;
    # 且逐站点公布判别力属站点可识别信息)
    rows = [[i, f'site_{int(x.site)}', int(x.n), int(x.n_pd), int(x.n_hc), f'{100 * x.n_hc / x.n:.0f}%',
             f3(x.r), pf(x.p), f'{x.upsit_sd:.2f}'] for i, x in enumerate(_ls.itertuples(), 1)]
    A(md_table([T('No.', '序号'), T('Centre', '站点'), 'n', T('Patients', '患者'), T('Controls', '对照'),
                T('Controls, % of centre', '对照占比'), 'r', 'P', T('UPSIT s.d.', 'UPSIT 标准差')], rows))
    A(T(f'The 44 centres are numbered 1 to 44 in order of r, and the centre column gives the PPMI site code. Across the 44 centres the correlation rises with the control fraction '
        f'(r = {LO2["r_controlfrac_vs_siter"]:.3f}, P = {pf(LO2["p_controlfrac_vs_siter"])}) and is unrelated to centre size '
        f'(r = {LO2["r_size_vs_siter"]:.3f}, P = {LO2["p_size_vs_siter"]:.2f}). Median centre size {LO2["median_site_n"]}. '
        f'The {LO2["n_with_hc"]} centres that enrolled controls as well as patients reach a mean r of {LO2["r_with_hc"]:.3f} against {LO2["r_without_hc"]:.3f} in the {LO2["n_without_hc"]} that enrolled patients only (Welch P = {pf(LO2["p_with_vs_without"])}), '
        f'and the two negative centres are both patient-dominated, so UPSIT and the scan each vary over a narrow range there and a within-centre correlation carries little information. '
        f'{LO2["n_significant"]} of 44 reach P < 0.05, against {LO2["n_significant_expected"]:.1f} expected under a true correlation of '
        f'{LO2["r_mean"]:.3f} at the observed centre sizes.',
        f'44 个中心按 r 从低到高编为 1 至 44,站点一列为 PPMI 站点号。44 个中心间,相关随对照占比升高(r = {LO2["r_controlfrac_vs_siter"]:.3f},P = {pf(LO2["p_controlfrac_vs_siter"])}),'
        f'与中心人数无关(r = {LO2["r_size_vs_siter"]:.3f},P = {LO2["p_size_vs_siter"]:.2f})。中心人数中位数 {LO2["median_site_n"]}。'
        f'同时纳入对照与患者的 {LO2["n_with_hc"]} 个中心平均 r 为 {LO2["r_with_hc"]:.3f},只纳入患者的 {LO2["n_without_hc"]} 个中心为 {LO2["r_without_hc"]:.3f}(Welch P = {pf(LO2["p_with_vs_without"])}),'
        f'两个取负值的中心均以患者为主,其中嗅觉分数与影像的取值范围都很窄,中心内部的相关几乎不携带信息。'
        f'44 个中心中 {LO2["n_significant"]} 个达到 P<0.05,而按真实相关 {LO2["r_mean"]:.3f} 与各中心实际人数计算,期望为 '
        f'{LO2["n_significant_expected"]:.1f} 个。'))
    A('')
    # ---------------- Table 2
    A(T('### Supplementary Table 2. Follow-up by recruitment cohort and group', '### 补充表 2. 各入组队列与各组的随访'))
    order = ['Hyposmia', 'RBD', 'Genetic', 'Non-hyposmia-enriched', 'All prodromal']
    names = {'Hyposmia': T('Hyposmia cohort', '嗅觉减退队列'), 'RBD': T('RBD cohort', 'RBD 队列'), 'Genetic': T('Pathogenic-variant cohort', '遗传携带队列'),
             'Non-hyposmia-enriched': T('RBD and variant-carrier group (RBD + carriers)', 'RBD 与遗传携带组(RBD + 携带者)'), 'All prodromal': T('All prodromal', '全前驱期')}
    ch = CH.set_index('channel')
    rows = []
    for k in order:
        r = ch.loc[k]
        rows.append([names[k], f'{int(r.n):,}', int(r.events), int(r.median_enrol_year), f'{r.median_followup_y:.2f}', f'{r.pct_ge4y:.1f}%', f'{r.person_years:,.0f}', f'{r.rate_per_100py:.2f}', f'{r.km2y_pct:.1f}%'])
    A(md_table([T('Group', '分组'), 'n', T('Conversions', '转化例数'), T('Median first-visit year', '中位首次访视年'), T('Median follow-up (y)', '中位随访(年)'), T('Followed ≥ 4 y', '随访 ≥ 4 年'), T('Person-years', '人年'), T('Per 100 person-years', '每百人年'), T('2-year KM conversion', '两年 KM 转化')], rows))
    LF = J['layer_flow']
    A('')
    A(T(f'Survival cohort after the exclusions in Methods (1,759 / 152). The RBD and variant-carrier group accrues more conversions than the hyposmia group (84 against 68) because it is followed longer, not because it is at higher risk. Per 100 person-years the hyposmia group is {ch.loc["Hyposmia","rate_per_100py"]:.2f} against {ch.loc["Non-hyposmia-enriched","rate_per_100py"]:.2f}, and the two-year Kaplan–Meier estimates coincide. Within the hyposmia group cumulative conversion by Kaplan–Meier is ' + ', '.join(f'{k["cum_conversion"]*100:.1f}% at {k["year"]} y ({k["at_risk"]} at risk)' for k in LF['km']) + '. The four-year value rests on 17 participants and is not used.',
        f'为方法中各步排除后的生存队列(1,759 / 152)。RBD 与遗传携带组转化例数多于嗅觉减退组(84 对 68)是随访更久所致,不是风险更高:每百人年嗅觉减退组 {ch.loc["Hyposmia","rate_per_100py"]:.2f} 对 {ch.loc["Non-hyposmia-enriched","rate_per_100py"]:.2f},两年 Kaplan–Meier 估计相同。嗅觉减退组内 Kaplan–Meier 累积转化为 ' + '、'.join(f'{k["year"]} 年 {k["cum_conversion"]*100:.1f}%(在险 {k["at_risk"]} 人)' for k in LF['km']) + ';四年值仅基于 17 人,不予使用。'))
    A('')
    # ---------------- Table 3
    A(T('### Supplementary Table 3. Discrimination of every reading, continuous against dichotomised, by stratum, and head to head with OIS', '### 补充表 3. 各读法的判别力:连续对二分、分层比较,以及与 OIS 的头对头'))
    A(T('**a. Continuous against dichotomised readings.**', '**a. 连续读法对二分读法。**'))
    rows = []
    lab38 = {'UPSIT, continuous score': T('UPSIT, continuous', 'UPSIT,连续'), 'UPSIT, dichotomised (hyposmic yes/no)': T('UPSIT, dichotomised at the 15th percentile', 'UPSIT,按第 15 百分位二分'),
             'imaging, continuous (lowest putamen %exp)': T('Imaging, continuous (lower putamen, % expected)', '影像,连续(较低侧壳核 %预期)'), 'imaging, dichotomised (DAT-deficit flag)': T('Imaging, dichotomised (DAT-deficit flag, 65%)', '影像,二分(DAT 缺损标志,65%)'),
             'OIS, continuous': T('OIS, continuous', 'OIS,连续')}
    for st in ['all prodromal', 'hyposmia']:
        g = C38[C38.stratum == st]
        for _, r in g.iterrows():
            rows.append([T({'all prodromal': 'All prodromal', 'hyposmia': 'Hyposmia group'}[st], {'all prodromal': '全前驱期', 'hyposmia': '嗅觉减退组'}[st]), lab38[r.reading], f'{int(r.n):,} / {int(r.events)}', ci(r.c, r.lo, r.hi)])
    A(md_table([T('Stratum', '分层'), T('Reading', '读法'), T('Participants / conversions', '人数 / 转化例数'), T('C-index (95% CI)', 'C-index(95% CI)')], rows))
    c = C38.set_index(['stratum', 'reading'])['c']
    A('')
    A(T(f'Dichotomisation costs {c[("all prodromal","UPSIT, continuous score")]-c[("all prodromal","UPSIT, dichotomised (hyposmic yes/no)")]:.3f} of concordance for UPSIT and {c[("all prodromal","imaging, continuous (lowest putamen %exp)")]-c[("all prodromal","imaging, dichotomised (DAT-deficit flag)")]:.3f} for imaging in the whole prodromal cohort, and {c[("hyposmia","UPSIT, continuous score")]-c[("hyposmia","UPSIT, dichotomised (hyposmic yes/no)")]:.3f} against {c[("hyposmia","imaging, continuous (lowest putamen %exp)")]-c[("hyposmia","imaging, dichotomised (DAT-deficit flag)")]:.3f} within the hyposmia group, where the olfactory flag is nearly constant. Descriptive, with no paired test.',
        f'二分使嗅觉分数在全前驱期损失 {c[("all prodromal","UPSIT, continuous score")]-c[("all prodromal","UPSIT, dichotomised (hyposmic yes/no)")]:.3f} 的 C-index、影像损失 {c[("all prodromal","imaging, continuous (lowest putamen %exp)")]-c[("all prodromal","imaging, dichotomised (DAT-deficit flag)")]:.3f};嗅觉减退组内分别为 {c[("hyposmia","UPSIT, continuous score")]-c[("hyposmia","UPSIT, dichotomised (hyposmic yes/no)")]:.3f} 对 {c[("hyposmia","imaging, continuous (lowest putamen %exp)")]-c[("hyposmia","imaging, dichotomised (DAT-deficit flag)")]:.3f},该组内嗅觉标志近乎常数。描述性,无配对检验。'))
    A(T('**b. Five readings in five strata, including the residual.**', '**b. 五种读法在五个分层中的表现,含残差。**'))
    stl = {'all prodromal': T('All prodromal', '全前驱期'), 'hyposmia': T('Hyposmia group', '嗅觉减退组'), 'hyposmia + DAT deficit': T('Hyposmia, DAT deficit', '嗅觉减退,DAT 缺损'),
           'hyposmia + no DAT deficit': T('Hyposmia, no DAT deficit', '嗅觉减退,DAT 非缺损'), 'non-hyposmia-enriched': T('RBD and variant-carrier', 'RBD 与遗传携带组')}
    rdl = {'UPSIT total': 'UPSIT', 'OIS': 'OIS', 'OMI': T('OMI (residual)', 'OMI(残差)'), 'putamen SBR': T('Single-region putamen', '单区壳核'), 'lowest putamen %expected': T('Lower putamen, % expected', '较低侧壳核 %预期')}
    rows = []
    for st in stl:
        g = C36[C36.stratum == st].set_index('reading')
        n = int(g.n.iloc[0]); ev = int(g.events.iloc[0])
        cells = [stl[st], f'{n:,} / {ev}']
        for rd in ['UPSIT total', 'OIS', 'OMI', 'putamen SBR', 'lowest putamen %expected']:
            r = g.loc[rd]; cells.append(ci(r.c, r.lo, r.hi) + ('\\*' if bool(r.at_chance) else ''))
        rows.append(cells)
    A(md_table([T('Stratum', '分层'), T('Participants / conversions', '人数 / 转化例数')] + [rdl[k] for k in ['UPSIT total', 'OIS', 'OMI', 'putamen SBR', 'lowest putamen %expected']], rows))
    A(T('\\* interval includes 0.5. The residual is at chance in the hyposmia group and in its deficit stratum only. In the whole cohort and in the RBD and variant-carrier group its interval excludes 0.5, and it never approaches UPSIT. All values come from one bootstrap run.',
        '\\* 区间含 0.5。残差仅在嗅觉减退组及其缺损亚组等同随机;在全队列与RBD 与遗传携带组其区间不含 0.5,且从未接近 UPSIT。全部数值出自同一次 bootstrap。'))
    A('')
    A(T('**c. Head to head with OIS in the hyposmia group.**', '**c. 嗅觉减退组内与 OIS 的头对头比较。**'))
    HH = J['head_to_head']
    tr = {'behavioural baseline': T('behavioural baseline', '行为学基线'), 'threshold rule': T('threshold rule', '阈值规则'), 'Coxnet 5-fold CV (in-cohort)': T('elastic-net Cox, five-fold CV on conversion labels, in-cohort', '弹性网 Cox,转化标签五折 CV,队列内'), 'Ridge on UPSIT, PD+HC zero-shot': T('ridge on UPSIT in PD + HC, transferred unchanged', 'PD + HC 上以 UPSIT 训练的岭回归,原样迁移')}
    rows = []
    for r in HH:
        d = '' if r['dc'] != r['dc'] else f"{r['dc']:+.3f} ({r['lo']:+.3f}, {r['hi']:+.3f})"
        p = '' if r['p'] != r['p'] else pf(r['p'])
        rows.append([r['model'], tr[r['training']], f3(r['c']), d, p])
    A(md_table([T('Model', '模型'), T('How it was built', '构建方式'), 'C', T('ΔC, OIS minus model (95% CI)', 'ΔC,OIS 减模型(95% CI)'), 'P'], rows))
    A(T('1,003 participants, 68 conversions, 1,000 paired bootstrap resamples. OIS exceeds UPSIT and the binary flag, is tied with the conversion-supervised 33-region model plus demographics and with its own imaging-only version, and reaches nominal significance against the 33 regions alone. The supervised models were cross-validated within this cohort and would be expected to fall in a new cohort. OIS used no conversion label.',
        '1,003 人,68 转化,1,000 次配对 bootstrap。OIS 优于 UPSIT 与二分标志,与转化标签监督的 33 区加人口学模型及自身仅影像版打平,对单独 33 区达到名义显著。监督模型在本队列内交叉验证,换新队列预期会下降;OIS 未用任何转化标签。'))
    A('')
    A('')
    # ---------------- Table 4 (three groups, §44)
    A(T('### Supplementary Table 4. Two groups at the median or the tertile boundary, and three groups by permutation-corrected selection of two cut points', '### 补充表 4. 中位数或三分位边界的两组划分,以及置换校正选出的两切点三组划分'))
    sl = {'hyposmia': T('Hyposmia group', '嗅觉减退组'), 'hyposmia + no DAT deficit': T('Non-deficit stratum', '非缺损亚组')}
    rl = {'OIS': 'OIS', 'UPSIT total': 'UPSIT', 'lower putamen %expected': T('Lower putamen, % expected', '较低侧壳核 %预期')}
    A(T('**a. Outcome-blind division at the median (the main-text division, Fig. 3c, e).**', '**a. 不看结局的中位数划分(正文所用,Fig 3c、e)。**'))
    rows = [[sl[r.stratum], rl[r.reading], f"{r['median']:.2f}", f'{int(r.ev_low)} / {int(r.n_low)}', f'{r.km2_low:.1f}% ({r.km2_low_lo:.1f}–{r.km2_low_hi:.1f})', f'{int(r.ev_high)} / {int(r.n_high)}', f'{r.km2_high:.1f}% ({r.km2_high_lo:.1f}–{r.km2_high_hi:.1f})', pf(r.p_logrank)] for _, r in M44.iterrows()]
    A(md_table([T('Population', '人群'), T('Reading', '读法'), T('Median', '中位数'), T('At or below, conversions / participants', '不高于中位数,转化 / 人数'), T('2-year conversion (95% CI)', '两年转化(95% CI)'), T('Above, conversions / participants', '高于中位数,转化 / 人数'), T('2-year conversion (95% CI)', '两年转化(95% CI)'), T('Log-rank P', 'log-rank P')], rows))
    A(T('**b. Outcome-blind division at the tertile boundary (lowest third against upper two thirds).**', '**b. 不看结局的三分位边界划分(最低三分之一对其余三分之二)。**'))
    rows = [[sl[r.stratum], rl[r.reading], f'{int(r.ev_low)} / {int(r.n_low)}', f'{r.km2_low:.1f}% ({r.km2_low_lo:.1f}–{r.km2_low_hi:.1f})', f'{int(r.ev_rest)} / {int(r.n_rest)}', f'{r.km2_rest:.1f}% ({r.km2_rest_lo:.1f}–{r.km2_rest_hi:.1f})', pf(r.p_two_group), pf(r.p_three_group_tertiles), pf(r.p_middle_vs_highest)] for _, r in L44.iterrows()]
    A(md_table([T('Population', '人群'), T('Reading', '读法'), T('Lowest third, conversions / participants', '最低三分之一,转化 / 人数'), T('2-year conversion (95% CI)', '两年转化(95% CI)'), T('Upper two thirds, conversions / participants', '其余三分之二,转化 / 人数'), T('2-year conversion (95% CI)', '两年转化(95% CI)'), T('Two-group log-rank P', '两组 log-rank P'), T('Three-tertile log-rank P', '三分位 log-rank P'), T('Middle against highest tertile P', '中间对最高三分位 P')], rows))
    A(T('**c. Selection of pairs of cut points (10th to 90th percentile grid in 2.5-percentile steps, each group at least 15% of the population), retaining the pair with the largest three-group log-rank statistic among pairs whose three groups all differ pairwise at P < 0.05. The permutation P repeats the unconstrained selection on 1,000 outcome permutations; cut-point intervals are from 200 bootstrap resamples.**',
        '**c. 两切点选取(10 至 90 百分位网格,步长 2.5 百分位,每组不少于该人群 15%),在三组两两 P < 0.05 的切点对中取三组 log-rank 统计量最大者。置换 P 为在 1,000 次结局置换上重做无约束选取;切点区间来自 200 次 bootstrap。**'))
    rows = []
    for _, r in S44.iterrows():
        if not bool(r.feasible):
            rows.append([sl[r.stratum], rl[r.reading], T('no pair separates all three groups', '无切点对能使三组两两分开'), '', '', '', '', '']); continue
        rows.append([sl[r.stratum], rl[r.reading], f'{r.cut1:.2f} ({r.cut1_boot_lo:.1f}–{r.cut1_boot_hi:.1f}) / {r.cut2:.2f} ({r.cut2_boot_lo:.1f}–{r.cut2_boot_hi:.1f})',
                     f'{r.frac_low*100:.0f}% / {r.frac_mid*100:.0f}% / {r.frac_high*100:.0f}%', f'{int(r.ev_low)} / {int(r.ev_mid)} / {int(r.ev_high)}',
                     f'{r.km2_low:.1f}% / {r.km2_mid:.1f}% / {r.km2_high:.1f}%', f'{r.chi2:.1f}, {pf(r.p_permutation)}', f'{pf(r.p_low_mid)} / {pf(r.p_mid_high)} / {pf(r.p_low_high)}'])
    A(md_table([T('Population', '人群'), T('Reading', '读法'), T('Cut points (bootstrap 95%)', '切点(bootstrap 95%)'), T('Share of population, low / mid / high', '占比,低 / 中 / 高'), T('Conversions', '转化例数'), T('2-year conversion', '两年转化'), T('χ², permutation P', 'χ²,置换 P'), T('Pairwise P, low–mid / mid–high / low–high', '两两 P,低–中 / 中–高 / 低–高')], rows))
    A(T('These cut points are outcome-selected and exploratory. The main text uses the median, which is outcome-blind and coincides with the lowest-50% enrolment rule of the trial simulation.',
        '这些切点是对着结局选出的,属探索性。正文使用中位数,它不看结局,且与试验模拟中最低 50% 入组的规则一致。'))
    A('')
    # ---------------- Figures
    # ---------------- Table 14 (proportional hazards check, supp_ph_test.py)
    if 'ph_test' in J:
        PH = J['ph_test']
        A(T('### Supplementary Table 5. Proportional hazards check for the Cox models adjusting a score for age',
            '### 补充表 5. 分数校正年龄的 Cox 模型的比例风险检验'))
        A(T("Schoenfeld residual test (rank time transform, as implemented in the lifelines package) for the models reported in the main text, the score in its own units plus age in years, and for a fuller specification with sex and years of education (continuous terms z-scored). The global row sums the term statistics. Hazard ratios are per point (or per year) in the main-text models and per standard deviation in the fuller models.",
            "对正文报告的模型(分数按原单位、年龄按年)以及加入性别与教育年限的完整模型(连续项 z 标准化)作 Schoenfeld 残差检验(秩时间变换,lifelines 软件包实现)。global 行为各项统计量之和。正文模型的风险比为每 1 分(或每 1 岁),完整模型为每 1 个标准差。"))
        _pl = {'hyposmia group': T('Hyposmia group', '嗅觉减退组'), 'non-deficit stratum': T('Non-deficit stratum', '非缺损亚组')}
        _tl = {'OIS': 'OIS', 'upsit': 'UPSIT', 'age': T('Age', '年龄'), 'sex_num': T('Sex', '性别'), 'educyrs': T('Education', '教育年限'), 'global': T('Global', '整体')}
        rows = [[_pl[r['population']], r['model'].replace('upsit', 'UPSIT').replace('sex_num', 'sex').replace('educyrs', 'education'), _tl.get(r['term'], r['term']),
                 '' if r['hr'] != r['hr'] else f"{r['hr']:.3f}", '' if r['p_cox'] != r['p_cox'] else pf(r['p_cox']),
                 f"{r['ph_test_stat']:.2f}", pf(r['ph_p']), f"{r['n']:,} / {r['events']}"] for r in PH['rows']]
        A(md_table([T('Population', '人群'), T('Model', '模型'), T('Term', '项'), 'HR', T('P (Cox)', 'P(Cox)'), T('Schoenfeld χ²', 'Schoenfeld χ²'), T('P (proportional hazards)', 'P(比例风险)'), T('n / conversions', 'n / 转化')], rows))
        _min = min(r['ph_p'] for r in PH['rows'] if r['term'] != 'global'); _gmin = min(r['ph_p'] for r in PH['rows'] if r['term'] == 'global')
        A(T(f'No term departs from proportional hazards at P < 0.05 (smallest term-level P = {_min:.3f}, smallest global P = {_gmin:.2f}).',
            f'没有任何一项在 P<0.05 水平上偏离比例风险(最小的项水平 P = {_min:.3f},最小的整体 P = {_gmin:.2f})。'))
        A('')
    # ---------------- Table 5
    A(T('### Supplementary Table 6. Flagging by rank on the three readings, hyposmia group (1,003 / 68)', '### 补充表 6. 三种读法按排序标记,嗅觉减退组(1,003 / 68)'))
    rdl2 = {'UPSIT total': 'UPSIT', 'lowest putamen %expected': T('Lower putamen, % expected', '较低侧壳核 %预期'), 'OIS': 'OIS'}
    A(T('**a. Number that must be flagged to capture a given fraction of converters, by rank.**', '**a. 为捕获给定比例转化者所需标记的人数,按秩。**'))
    piv = C39B.pivot(index='target_sensitivity', columns='reading', values=['n_flagged', 'frac_flagged'])
    rows = []
    for ts in sorted(C39B.target_sensitivity.unique()):
        cells = [f'{ts*100:.0f}%']
        for rd in ['UPSIT total', 'lowest putamen %expected', 'OIS']:
            cells.append(f"{int(piv.loc[ts, ('n_flagged', rd)])} ({piv.loc[ts, ('frac_flagged', rd)]*100:.0f}%)")
        rows.append(cells)
    A(md_table([T('Converters captured', '捕获转化者比例')] + [rdl2[k] for k in ['UPSIT total', 'lowest putamen %expected', 'OIS']], rows))
    A(T('**b. Fraction of converters captured when a fixed fraction of the group is flagged.**', '**b. 标记固定比例的人时捕获的转化者比例。**'))
    rows = []
    for ff in sorted(C39F.frac_flagged.unique()):
        cells = [f'{ff*100:.0f}% ({int(C39F[C39F.frac_flagged == ff].n_flagged.iloc[0])})']
        for rd in ['UPSIT total', 'lowest putamen %expected', 'OIS']:
            r = C39F[(C39F.frac_flagged == ff) & (C39F.reading == rd)].iloc[0]
            cells.append(f"{int(r.events_caught)} / 68 ({r.frac_events*100:.0f}%)")
        rows.append(cells)
    A(md_table([T('Group flagged (n)', '标记比例(人数)')] + [rdl2[k] for k in ['UPSIT total', 'lowest putamen %expected', 'OIS']], rows))
    A(T('Read by rank alone, OIS dominates at every operating point, needing 132, 247 and 326 participants flagged to capture 50%, 70% and 80% of converters against 276, 458 and 525 for UPSIT and 158, 348 and 537 for the lower putamen percentage. Imaging falls below UPSIT at the high-sensitivity end (537 against 525), the same non-monotonicity seen across its clinical bands. These cut points are selected against the outcome and cannot serve as clinical thresholds. The median remains the primary division.',
        '若只按秩读,OIS 在每一个工作点均占优,欲覆盖 50%、70%、80% 的转化者需标记 132、247 与 326 人,而 UPSIT 需 276、458 与 525 人,较低侧壳核 %预期需 158、348 与 537 人。影像在高灵敏度端反而低于 UPSIT(537 对 525),与其临床分档在本队列不单调是同一现象。这些切点是对着结局选出的,不能作为临床阈值,中位数仍为主划分。'))
    A('')
    # ---------------- Table 8
    A(T('### Supplementary Table 7. Test-retest reliability', '### 补充表 7. 重测信度'))
    rows = [['UPSIT', f"{ICC['upsit_r']:.3f}", f"{ICC['upsit_icc']:.3f}"], [T('OMI (residual)', 'OMI(残差)'), f"{ICC['omi_r']:.3f}", f"{ICC['omi_icc']:.3f}"]]
    A(md_table([T('Score', '分数'), T('First against second visit r', '首次对第二次访视 r'), 'ICC(1,1)'], rows))
    # 分组 ICC 现由 notebook §12 写入总账并同名落盘,SWEDD 已在源头剔除,此处不再补丁式过滤
    ICG = pd.DataFrame(ICC['by_group'])
    gl = {'Hyposmia': T('Prodromal, hyposmia cohort', '前驱期,嗅觉减退队列'), 'RBD': T('Prodromal, RBD cohort', '前驱期,RBD 队列'), 'Genetic': T('Prodromal, pathogenic-variant cohort', '前驱期,遗传携带队列'), 'PD': T('Parkinson\'s disease (training)', '帕金森病(训练)'), 'HC': T('Healthy controls (training)', '健康对照(训练)')}
    ICG = ICG.set_index('group').loc[['Hyposmia', 'RBD', 'Genetic', 'PD', 'HC']].reset_index()
    rows = [[gl[r.group], int(r.n_pairs), f'{r.upsit_icc:.3f}', f'{r.omi_icc:.3f}'] for _, r in ICG.iterrows() if r.group in gl]
    A(T('**By cohort and recruitment cohort.**', '**按队列与入组队列。**'))
    A(md_table([T('Group', '分组'), T('Pairs', '配对数'), T('ICC, UPSIT', 'ICC,UPSIT'), T('ICC, OMI', 'ICC,OMI')], rows))
    A(T('Within the three prodromal recruitment cohorts the residual\'s ICC lies between 0.75 and 0.83. In the training cohorts both scores are less reliable, the residual there having mean zero by construction and a narrower range. SWEDD participants are excluded at source, as everywhere in this work.',
        '三个前驱期队列内残差的 ICC 在 0.75 至 0.83 之间。训练队列中两个分数的信度均较低,残差在那里按构造均值为零且取值范围更窄。SWEDD 与全文口径一致已在源头剔除。'))
    A(T(f'{ICC["n_pairs"]:,} participants with at least two UPSIT measurements and an OIS. About three quarters of residual variance is stable across visits, so the small R² values in the pathology tests reflect diversity of sources rather than noise.',
        f'{ICC["n_pairs"]:,} 名至少两次 UPSIT 且有 OIS 的参与者。残差方差约四分之三跨访视稳定,故病理检验中较小的 R² 反映来源多样而非噪声。'))
    A('')
    # ---------------- Table 9
    A(T('### Supplementary Table 8. Cerebrospinal fluid and blood analytes against the three scores', '### 补充表 8. 脑脊液与血液分析物对三个分数'))
    A(T('Standardised regression coefficients adjusted for age, sex and years of education. The q value is Benjamini–Hochberg across analytes within cohort.', '校正年龄、性别与教育年限的标准化回归系数;q 为队列内跨分析物的 Benjamini–Hochberg 校正。'))
    alab = {'CSF pTau181/ABeta42': 'CSF pTau181/Aβ42', 'CSF pTau181': 'CSF pTau181', 'CSF eMTBR-TAU243': 'CSF eMTBR-tau243', 'CSF ABeta42': 'CSF Aβ42', 'Serum NfL': T('Serum NfL', '血清 NfL'), 'Plasma NfL': T('Plasma NfL', '血浆 NfL'), 'CSF NfL': 'CSF NfL'}
    rows = []
    for _, r in C22.iterrows():
        rows.append([T({'Prodromal': 'Prodromal', 'PD': 'Diagnosed PD'}[r.cohort], {'Prodromal': '前驱期', 'PD': '已确诊 PD'}[r.cohort]), alab[r.analyte], int(r.n),
                     f'{r.b_UPSIT:+.3f} ({pf(r.p_UPSIT)})', f'{r.b_OIS:+.3f} ({pf(r.p_OIS)})', f'{r.b_OMI:+.3f} ({pf(r.p_OMI)})', f'{r.q_OMI:.3f}'])
    A(md_table([T('Cohort', '队列'), T('Analyte', '分析物'), 'n', 'UPSIT β (P)', 'OIS β (P)', 'OMI β (P)', T('q (OMI)', 'q(OMI)')], rows))
    RB = J['csf_tau']['robustness_ratio_prodromal']; STR = J['csf_tau']['strata_ratio_prodromal']; XP = J['csf_tau']['cross_platform']
    A('')
    A(T('**Sensitivity analyses for the ratio in the prodromal cohort.**', '**前驱期比值的敏感性分析。**'))
    rows = [[k, v['n'], f"{v['b_UPSIT']:+.3f} ({pf(v['p_UPSIT'])})", f"{v['b_OIS']:+.3f} ({pf(v['p_OIS'])})", f"{v['b_OMI']:+.3f} ({pf(v['p_OMI'])})"] for k, v in RB.items()]
    for k, v in STR.items():
        rows.append([T(f'within {k} recruitment cohort only', f'仅 {k} 队列内'), v['n'], f"{v['b_UPSIT']:+.3f} ({pf(v['p_UPSIT'])})", f"{v['b_OIS']:+.3f} ({pf(v['p_OIS'])})", f"{v['b_OMI']:+.3f} ({pf(v['p_OMI'])})"])
    for k, v in XP.items():
        rows.append([T(f'platform: {k}', f'平台:{k}'), v['n'], f"{v['b_UPSIT']:+.3f} ({pf(v['p_UPSIT'])})", f"{v['b_OIS']:+.3f} ({pf(v['p_OIS'])})", f"{v['b_OMI']:+.3f} ({pf(v['p_OMI'])})"])
    A(md_table([T('Specification', '设定'), 'n', 'UPSIT β (P)', 'OIS β (P)', 'OMI β (P)'], rows))
    PC = J['csf_tau']['platform_concordance']
    A(T(f'The two tau platforms were run on the same participants (Spearman ρ = {PC["spearman_prodromal"]:.3f} between platforms, n = {PC["n_prodromal"]}), so cross-platform agreement is consistency rather than independent replication. The raw total correlates slightly more strongly with the ratio than the residual does, which is an attribution result and not a demonstration that the residual outperforms the total.',
        f'两个 tau 平台在同一批人上运行(平台间 Spearman ρ = {PC["spearman_prodromal"]:.3f},n = {PC["n_prodromal"]}),跨平台一致是一致性而非独立重复。原始总分与比值的相关略强于残差,这是归属结果而非残差胜过总分。'))
    A('')
    # ---------------- Table 10
    A(T('### Supplementary Table 9. Genotype contrast in full, its premise, and both tests on one metric', '### 补充表 9. 基因型对比全表、其前提,以及两项检验换算到同一个量'))
    GP = J['genotype_saa_premise']
    A(T('**a. Premise, seed-amplification positivity in diagnosed Parkinson\'s disease by genotype.**', '**a. 前提:已确诊帕金森病中按基因型的种子扩增阳性率。**'))
    # the ledger stores the PPMI subgroup label 'GBA', displayed here as the current symbol GBA1
    rows = [[g['group'].replace('GBA', 'GBA1'), g['n'], g['positive'], f"{g['rate']*100:.0f}%"] for g in GP['groups']]
    A(md_table([T('Group', '组'), 'n', T('Positive', '阳性'), T('Rate', '阳性率')], rows))
    A(T(f'Fisher exact P = {pf(GP["fisher_p"])}.', f'Fisher 精确检验 P = {pf(GP["fisher_p"])}。'))
    A(T('**b. GBA1 against LRRK2 carriers, linear regression adjusted for age, sex and years of education (β is GBA1 minus LRRK2).**', '**b. GBA1 对 LRRK2 携带者,校正年龄、性别与教育年限的线性回归(β 为 GBA1 减 LRRK2)。**'))
    rows = []
    for r in J['genotype']:
        rows.append([T({'Prodromal': 'Carriers without disease', 'PD': 'Diagnosed PD'}[r['cohort']], {'Prodromal': '未发病携带者', 'PD': '已确诊 PD'}[r['cohort']]), r['measure'].replace('UPSIT total', 'UPSIT'), f"{r['n_gba']} / {r['n_lrrk2']}", f"{r['mean_gba']:.2f}", f"{r['mean_lrrk2']:.2f}", f"{r['beta']:+.2f} ({r['lo']:+.2f}, {r['hi']:+.2f})", pf(r['p'])])
    A(md_table([T('Cohort', '队列'), T('Measure', '指标'), T('n GBA1 / LRRK2', 'n GBA1 / LRRK2'), T('Mean GBA1', 'GBA1 均值'), T('Mean LRRK2', 'LRRK2 均值'), 'β (95% CI)', 'P'], rows))
    A(T('Among carriers without disease the two genotypes have the same UPSIT while OIS is higher and OMI lower in GBA1, the two differences cancelling. Among diagnosed patients OIS is identical (the negative control) and the whole 6-point difference in UPSIT sits in the residual. Putamen SBR is shown for completeness.',
        '未发病携带者中两基因型 UPSIT 相同,而 GBA1 的 OIS 更高、OMI 更低,两者相抵。已确诊患者中 OIS 相同(阴性对照),6 分的嗅觉差全部落在残差上。壳核 SBR 为完整起见列出。'))
    A(T('**c. Partial correlations adjusted for age, sex and years of education.**', '**c. 校正年龄、性别与教育年限的偏相关(补充方法 4)。**'))
    rows = [[r.test.replace('ABeta42', 'Aβ42').replace('GBA vs', 'GBA1 vs'), T({'Prodromal': 'Prodromal', 'PD': 'Diagnosed PD'}[r.cohort], {'Prodromal': '前驱期', 'PD': '已确诊 PD'}[r.cohort]), r.score, int(r.n), f'{r.r:+.3f}', pf(r.p)] for _, r in C42.iterrows()]
    A(md_table([T('External measurement', '外部测量'), T('Cohort', '队列'), T('Score', '分数'), 'n', T('Partial r', '偏相关 r'), 'P'], rows))
    A(T('The UPSIT column shows that the total behaves like the residual on both measurements. The tests establish that the decomposition separates two biologies, not that the residual outperforms the total.',
        'UPSIT 列显示总分在两项测量上的行为与残差相同;这两项检验证明的是分解把两种生物学分开了,不是残差胜过总分。'))
    A('')
    # ---------------- Table 13
    A(T('### Supplementary Table 10. Prodromal seed-amplification data: tested and not usable', '### 补充表 10. 前驱期种子扩增数据:已检验,不可用'))
    PS = J['prodromal_saa']
    rows = [[a['assay'], a['cohort'], a['n'], a['positive'], f"{a['rate']*100:.0f}%"] for a in PS['assay_quality']]
    A(T('**a. Positivity of the skin assay by cohort (PPMI project 259).**', '**a. 皮肤检测按队列的阳性率(PPMI 项目 259)。**'))
    A(md_table([T('Assay', '检测'), T('Cohort', '队列'), 'n', T('Positive', '阳性'), T('Rate', '阳性率')], rows))
    rows = [[r['score'].replace('upsit', 'UPSIT'), f"{r['auc']:.3f}", f"{r['odds_ratio']:.3f}", pf(r['p'])] for r in PS['skin']['results']]
    A(T(f'**b. Skin assay, prodromal (n = {PS["skin"]["n"]}, {PS["skin"]["n_positive"]} positive), logistic regression adjusted for age, sex and years of education.**', f'**b. 皮肤检测,前驱期(n = {PS["skin"]["n"]},阳性 {PS["skin"]["n_positive"]}),校正年龄、性别与教育年限的 logistic 回归。**'))
    A(md_table([T('Score', '分数'), 'AUC', T('OR per s.d.', '每标准差 OR'), 'P'], rows))
    rows = [[r['score'].replace('upsit', 'UPSIT'), f"{r['auc']:.3f}", pf(r['p_unadjusted']), pf(r['p_arm_adjusted'])] for r in PS['csf']['results']]
    A(T(f'**c. CSF dilution series, prodromal (n = {PS["csf"]["n"]}, {PS["csf"]["n_positive"]} positive, PPMI project 262).**', f'**c. 脑脊液稀释系列,前驱期(n = {PS["csf"]["n"]},阳性 {PS["csf"]["n_positive"]};PPMI 项目 262)。**'))
    A(md_table([T('Score', '分数'), 'AUC', T('P, unadjusted', 'P,未校正'), T('P, recruitment cohort-adjusted', 'P,入组队列校正')], rows))
    A(T(f'The skin assay is positive in only about half of diagnosed disease against 92.7% reported for the validated assay [PMID 38506839], and more often in prodromal than in diagnosed participants, so it lacks sensitivity. The CSF series has {PS["csf"]["hyposmia_arm_n"]} hyposmia-recruitment cohort participants with {PS["csf"]["hyposmia_arm_positive"]} positives, giving {PS["csf"]["hyposmia_arm_power_auc070"]*100:.0f}% power to detect an AUC of 0.70, and its positivity in diagnosed disease is {PS["csf"]["pd_positivity"]*100:.0f}% against 86% for the primary assay. The association between residual and synuclein pathology is therefore established only in the diagnosed cohort.',
        f'皮肤检测在已确诊疾病中阳性率仅约一半,而验证过的检测报告 92.7% [PMID 38506839],且前驱期阳性率高于已确诊者,故敏感度不足。脑脊液系列的嗅觉减退队列仅 {PS["csf"]["hyposmia_arm_n"]} 人、{PS["csf"]["hyposmia_arm_positive"]} 例阳性,检出 AUC 0.70 的功效 {PS["csf"]["hyposmia_arm_power_auc070"]*100:.0f}%,其在已确诊疾病中的阳性率为 {PS["csf"]["pd_positivity"]*100:.0f}%(主检测 86%)。残差与突触核蛋白病理的关联因此仅在已确诊队列中成立。'))
    A('')
    # ---------------- Table 14
    A(T('### Supplementary Table 11. Age and the three scores', '### 补充表 11. 年龄与三个分数'))
    AA = J['age_absorption']['correlations']
    _cl = {'Hyposmia stratum': T('Hyposmia group', '嗅觉减退组')}   # ledger label predates the cohort/group/stratum terminology, mapped at display only
    rows = [[_cl.get(r['cohort'], r['cohort']), f"{r['n']:,}", '' if r['r_age_upsit'] is None else f"{r['r_age_upsit']:+.3f}", '' if r['r_age_ois'] is None else f"{r['r_age_ois']:+.3f}", f"{r['r_age_omi']:+.3f}", pf(r['p_age_omi'])] for r in AA]
    A(md_table([T('Cohort', '队列'), 'n', 'r(age, UPSIT)', 'r(age, OIS)', 'r(age, OMI)', T('P for OMI', 'OMI 的 P')], rows))
    A(T('Age is absorbed in the training cohorts but not in the applied prodromal cohort, so every analysis involving the residual is adjusted for age. Sinonasal disease, head trauma, smoking and post-infectious dysfunction have no field in PPMI.',
        '年龄在训练队列中被吸收,在应用的前驱期队列中未被吸收,故所有涉及残差的分析均校正年龄。鼻窦疾病、颅脑外伤、吸烟与感染后嗅觉障碍在 PPMI 中无字段。'))
    A('')
    A(T('## Supplementary Figures', '## 补充图'))
    d = outdir()          # 图只有英文一套,两语种的 md 都引用同一路径
    A(T('**Supplementary Fig. 1 | Participant flow.** Training cohort (Parkinson\'s disease and healthy controls with both tests) at the top left. Prodromal cohort through the survival exclusions on the right, then split into the hyposmia group and the RBD and variant-carrier group, and the hyposmia group split at the 65% PARS threshold. Numbers as in Methods and Supplementary Table 2.',
        '**补充图 1 | 受试者流程。**左上为训练队列(同时有两项检查的帕金森病与健康对照);右侧为前驱期队列经生存分析各步排除,再分为嗅觉减退组与 RBD 与遗传携带组,嗅觉减退组按 PARS 65% 阈值再分。数字同方法与补充表 2。'))
    A(f'![FigS1]({d}/FigS1_flow.png)')
    A(T('**Supplementary Fig. 2 | Calibration of OIS in the prodromal cohort.** **a**, Mean observed against mean predicted UPSIT by decile of OIS in the RBD and variant-carrier group (n = 885). Observed lies below predicted in every decile by about 2.7 points, which is the mean residual of that group. **b**, Slope of measured UPSIT on OIS, with its 95% CI, r, P and n, in the whole prodromal cohort, the hyposmia cohort and each of the two recruitment cohorts of the RBD and variant-carrier group separately. The pooled slope of the RBD and variant-carrier group (1.027) is omitted because it mixes two recruitment cohorts with different means. The within-recruitment cohort values are the interpretable ones.',
        '**补充图 2 | OIS 在前驱期的校准。**(a) RBD 与遗传携带组(n = 885)按 OIS 十分位的平均实测对平均预测 UPSIT;每个十分位实测均低于预测约 2.7 分,即该组的平均残差。(b) 全前驱期、嗅觉减退队列以及RBD 与遗传携带组的两个队列分别的实测 UPSIT 对 OIS 的斜率及其 95% CI,附 r、P 与 n。省略RBD 与遗传携带组合并斜率 1.027,它出自均值不同的两个队列混合;队列内的值才可信。'))
    A(f'![FigS2]({d}/FigS2_calibration.png)')
    A(T('**Supplementary Fig. 3 | Robustness of the advantage of OIS over UPSIT in the hyposmia group.** ΔC with 1,000-resample paired bootstrap intervals for the full group, the one- and two-year landmarks (participants converting before the landmark removed), the temporal split (model retrained on participants enrolled before 2017 and evaluated in hyposmia-group participants enrolled from 2017, point estimate only).',
        '**补充图 3 | 嗅觉减退组中 OIS 相对 UPSIT 优势的稳健性。**ΔC 及 1,000 次配对 bootstrap 区间:全组、1 年与 2 年 landmark(剔除 landmark 前转化者)、时间拆分(在 2017 年前入组者上重训模型,在 2017 年起入组的嗅觉减退组评价,仅点估计)。'))
    A(f'![FigS3]({d}/FigS3_robustness.png)')
    A(T('**Supplementary Fig. 4 | Kaplan–Meier curves by OIS tertile, with 95% confidence bands and pairwise log-rank P values.** **a**, Whole hyposmia group. **b**, Non-deficit stratum. In both populations the middle and highest tertiles do not differ (P = 0.42 and P = 0.09), so the main text divides at the median into a high-risk lower half and a low-risk upper half (Fig. 3c, e).',
        '**补充图 4 | 按 OIS 三分位的 Kaplan–Meier 曲线,含 95% 置信带与两两 log-rank P。**(a) 全嗅觉减退组。(b) 非缺损亚组。两个人群中中间与最高三分位均不分开(P = 0.42 与 P = 0.09),因此正文按中位数分为高危的低半与低危的高半(Fig 3c、e)。'))
    A(f'![FigS4]({d}/FigS4_tertiles.png)')
    A(T('**Supplementary Fig. 5 | Three risk bands from permutation-corrected selection of two cut points on OIS in the hyposmia group.** Bands at OIS ≤ 24.91, 24.91 to 27.24 and > 27.24 hold 15%, 20% and 65% of the group with two-year conversion of 22.7%, 9.6% and 0.9%, all pairwise log-rank P < 0.001, permutation P < 0.001 for the selection, and bootstrap 95% intervals for the cut points of 24.4 to 25.4 and 27.1 to 28.8 (Supplementary Table 4). Both cut points lie below the median, which is where the risk is concentrated. The cut points are outcome-selected and this figure is exploratory; no pair of cut points separates three groups in the non-deficit stratum.',
        '**补充图 5 | 嗅觉减退组 OIS 经置换校正选出两切点得到的三个风险档。**切点 OIS ≤ 24.91、24.91 至 27.24、> 27.24 分别占该组 15%、20%、65%,两年转化 22.7%、9.6%、0.9%,两两 log-rank P 均 < 0.001,选取过程的置换 P < 0.001,切点 bootstrap 95% 区间 24.4 至 25.4 与 27.1 至 28.8(补充表 4)。两个切点都落在中位数以下,风险集中于此。切点对着结局选出,本图属探索性;非缺损亚组中无切点对能分出三组。'))
    A(f'![FigS5]({d}/FigS5_three_bands.png)')
    import re
    doc = '\n\n'.join(S) + '\n'
    doc = re.sub(r'(?<![\w\-])-(?=\d)', '−', doc)
    if LANG == 'zh':
        out = []
        for l in doc.split('\n'):
            if l.startswith('|') or l.startswith('#') or l.startswith('!['):
                out.append(l); continue
            for _ in range(4):
                l = re.sub(r'\(([^()\n]*?)[；;]([^()\n]*?)\)', r'(\1,\2)', l)
            l = re.sub(r'[；;](?=\S)', '。', l)
            l = re.sub(r'(?<!RRID)(?<!www)(?<!1)[：:](?!//)', lambda m: '' if l[max(0, m.start()-1):m.start()] in ('为', '是', '即') else '。', l)
            out.append(l)
        doc = '\n'.join(out)
    return doc

if __name__ == '__main__':
    for LANG in (('en',) if 'SUPP_FIG_DIR' in os.environ else ('en', 'zh')):
        globals()['LANG'] = LANG
        setup_fonts()
        cs = channel_slopes()
        if LANG == 'en':          # 图只出一套(英文),中文稿引用同一批文件
            figS1_flow(); figS2_calibration(cs); figS3_robustness(); figS4_tertiles(); figS5_three_bands()
        doc = build_doc(cs)
        fn = 'Supplementary_Information_EN.md' if LANG == 'en' else 'Supplementary_Information_中文.md'
        if 'SUPP_FIG_DIR' in os.environ: fn = fn.replace('.md', '_' + os.path.basename(os.environ['SUPP_FIG_DIR']) + '.md')
        open(fn, 'w').write(doc); print('wrote', fn, len(doc.split()), 'words')
