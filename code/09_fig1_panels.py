#!/usr/bin/env python3
"""09_fig1_panels.py  Data panels of Fig. 1d, drawn from the ledger so the figure carries measured values.

LEFT   Kaplan-Meier schematic, two groups, matching the median split of Fig 3c.
MIDDLE Residual and fitted value across cerebrospinal fluid tau tertiles, measured values
       from section 45 with the across-group P printed.
RIGHT  Both external tests on one common metric, partial correlations from section 42.

Each panel is written as its own transparent, Illustrator-editable PDF so it can be
dropped into the Fig 1 artwork.

    python3 code/09_fig1_panels.py
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Illustrator opens matplotlib PDFs as uneditable when the font is embedded as a CID
# subset with Identity-H encoding. The base-14 core fonts use WinAnsiEncoding instead,
# which Illustrator treats as ordinary live type. This restricts text to WinAnsi.
plt.rcParams.update({
    'pdf.use14corefonts': True, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'font.family': 'serif',
    'font.serif': ['Times', 'Times New Roman', 'Nimbus Roman', 'Liberation Serif', 'DejaVu Serif'],
    'figure.facecolor': 'none', 'axes.facecolor': 'none',
    'savefig.facecolor': 'none', 'savefig.edgecolor': 'none',
    'axes.unicode_minus': False, 'font.size': 8,
})
NAVY, TERRA, GREY, TEAL = '#2e3a5c', '#d4624a', '#6b6b6b', '#16A085'
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); os.makedirs('figures/main', exist_ok=True)
J = json.load(open('results/repro/section16_stats.json'))


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(f'figures/main/{name}.{ext}', dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)
    print('  wrote', name)


def pfmt(p):
    return 'P < 0.001' if p < 0.001 else f'P = {p:.3f}' if p < 0.01 else f'P = {p:.2f}'


# ----------------------------------------------------------------- left, KM schematic
def km_two_groups():
    fig, ax = plt.subplots(figsize=(3.0, 2.5))
    t = np.linspace(0, 2, 200)
    # idealised shapes only; the measured two-year values are in Fig 3c
    # decay constants chosen so the two-year end points sit at the measured values of Fig. 3c
    for lab, col, k, lw in [('high risk\n(OIS at or below median)', TERRA, 0.0335, 2.2),
                            ('low risk\n(OIS above median)', TEAL, 0.0026, 2.2)]:
        y = 100 * np.exp(-k * t ** 1.8)
        ax.step(t, y, where='post', color=col, lw=lw)
        ax.annotate(lab, xy=(2.0, y[-1]), xytext=(6, 0), textcoords='offset points',
                    fontsize=6.6, color=col, va='center', fontweight='bold')
    ax.set_xlim(0, 2.6)
    ax.set_ylim(86, 100.5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('time', fontsize=7.5)
    ax.set_ylabel('free of Parkinson\'s disease', fontsize=7.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.text(0.0, -0.16, 'OIS separates who converts', transform=ax.transAxes,
            fontsize=7.4, fontweight='bold', color=NAVY)
    ax.text(0.0, -0.28, 'shapes only, measured values in Fig. 3c',
            transform=ax.transAxes, fontsize=6.2, color=GREY, style='italic')
    save(fig, 'Fig1d_left_KM')


# ----------------------------------------------------------------- middle, tau gradient
def tau_gradient():
    T = pd.DataFrame(J['tau_gradient_fig1d'])
    T = T[(T.split == 'tertile')]
    fig, axes = plt.subplots(2, 1, figsize=(3.1, 2.9), sharex=True,
                             gridspec_kw=dict(hspace=.35))
    labels = ['low', 'middle', 'high']
    for ax, score, col in [(axes[0], 'OMI', TERRA), (axes[1], 'OIS', NAVY)]:
        d = T[T.score == score].sort_values('group')
        y = d.adj_mean.values; e = d.adj_se.values
        ax.bar(np.arange(3), y, 0.6, color=col, yerr=e,
               error_kw=dict(ecolor='#555', lw=.8, capsize=2))
        for i, (v, s) in enumerate(zip(y, e)):
            ax.text(i, v - s - (0.28 if score == 'OMI' else -0.28), f'{v:.2f}',
                    ha='center', va='top' if score == 'OMI' else 'bottom',
                    fontsize=6.6, color=col, fontweight='bold')
        p = float(d.p_across_groups.iloc[0])
        ax.text(0.02, 0.06 if score == 'OMI' else 0.80, f'across groups {pfmt(p)}',
                transform=ax.transAxes, ha='left', fontsize=6.4,
                color=col if p < 0.05 else GREY, fontweight='bold' if p < 0.05 else 'normal')
        ax.set_ylabel(f'{score}\n(adjusted mean)', fontsize=7, color=col)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        if score == 'OMI':
            ax.set_ylim(min(y - e) - 2.0, 0)
            ax.axhline(0, color='#888', lw=.8)
        else:
            ax.set_ylim(28.0, 30.2)
    axes[1].set_xticks(np.arange(3)); axes[1].set_xticklabels(labels, fontsize=7.5)
    axes[1].set_xlabel('cerebrospinal fluid tau burden (tertile)', fontsize=7.5)
    axes[0].text(0.0, 1.16, 'the residual falls as tau burden rises,\nthe fitted value does not move',
                 transform=axes[0].transAxes, fontsize=7.2, fontweight='bold', color='#333a45')
    axes[1].text(0.0, -0.62, f'n = {int(T[T.score == "OMI"].n.sum())} prodromal participants, means adjusted\nfor age, sex and years of education',
                 transform=axes[1].transAxes, fontsize=6.2, color=GREY, style='italic')
    save(fig, 'Fig1d_mid_tau')


# ----------------------------------------------------------------- right, one metric
def one_metric():
    P = {(r['test'], r['score']): r for r in J['partial_r_two_tests']}
    # these are ledger keys and keep the PPMI label, the axis text below shows GBA1
    tests = ['CSF pTau181/ABeta42', 'GBA vs LRRK2 genotype']
    fig, ax = plt.subplots(figsize=(3.1, 2.6))
    x = np.arange(2); w = 0.3
    for k, (score, col, lab) in enumerate([('OMI', TERRA, 'OMI, not explained by the scan'),
                                           ('OIS', NAVY, 'OIS, explained by the scan')]):
        h = [abs(P[(t, score)]['r']) for t in tests]
        pos = x + (k - 0.5) * w
        ax.bar(pos, h, w * 0.92, color=col, label=lab)
        for xi, t in zip(pos, tests):
            r = P[(t, score)]
            ax.text(xi, abs(r['r']) + 0.012, f"{abs(r['r']):.2f}\n{pfmt(r['p'])}",
                    ha='center', fontsize=6.2, color=col, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['CSF pTau181 / A-beta-42\nprodromal, n = 827',
                        'GBA1 vs LRRK2 carriers\nmanifest PD, n = 69 vs 126'], fontsize=6.8)
    ax.set_ylim(0, 0.46); ax.set_yticks([0, 0.2, 0.4])
    ax.set_ylabel('partial correlation with the\nexternal measurement  |r|', fontsize=7.2)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=6.2, frameon=False, loc='upper left', handlelength=1.1,
              handleheight=0.9, borderpad=0.0, labelspacing=0.3)
    ax.text(0.0, -0.30, 'adjusted for age, sex and years of education',
            transform=ax.transAxes, fontsize=6.2, color=GREY, style='italic')
    save(fig, 'Fig1d_right_metric')


if __name__ == '__main__':
    km_two_groups()
    for key, fn, sec in [('tau_gradient_fig1d', tau_gradient, '45'), ('partial_r_two_tests', one_metric, '42')]:
        if key in J:
            fn()
        else:
            print(f'  section {sec} not yet in the ledger, skipping that panel')
