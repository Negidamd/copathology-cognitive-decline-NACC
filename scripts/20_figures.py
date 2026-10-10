#!/usr/bin/env python3
"""20_figures.py — Figures 1-4.
  Fig 1  counts among the 3,391 participants with all four proteinopathies classifiable (R1.13)
  Fig 2a observed means in a FIXED-COMPOSITION subsample followed >= 6 years, groups with >= 50
         participants, so survivorship cannot produce apparent improvement (R1.7)
  Fig 3  Model 3 coefficients (order as in the original)
  Fig 4  unadjusted pairwise rates (original definition; cells with n < 10 left blank)
PNG longest side <= 1800 px; PDF alongside."""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns

HERE = Path(__file__).resolve().parents[1]
FIG = HERE/'figures'; FIG.mkdir(exist_ok=True)
R = json.load(open(HERE/'tables'/'results.json'))
df = pd.read_csv(HERE/'data'/'analytical_dataset.csv', low_memory=False)
npm = pd.read_csv(HERE/'data'/'neuropath_master.csv', low_memory=False)
P = df[df.primary == 1].copy()
rates = pd.read_csv(HERE/'tables'/'individual_rates.csv')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8, 'legend.fontsize': 7, 'axes.linewidth': 0.6, 'lines.linewidth': 1.2,
                     'axes.spines.top': False, 'axes.spines.right': False})
BC = {0: '#7f7f7f', 1: '#4878d0', 2: '#ee854a', 3: '#6acc64', 4: '#d65f5f'}
LAB = lambda s: (s.replace('ADNC', 'AD').replace('TDP-43', 'TDP').replace('Other tauopathy', 'FTLD-tau')
                  .replace(' + ', '+').replace('None above threshold', 'None/minimal'))


def save(fig, name, dpi):
    fig.savefig(FIG/f'{name}.png', dpi=dpi, bbox_inches='tight'); fig.savefig(FIG/f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)

# Figure 1
comp = npm[npm.n_path_assessed == 4]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={'width_ratios': [1, 1.3]})
bc = comp.copath_count.value_counts().sort_index()
bars = ax[0].bar(bc.index, bc.values, color=[BC[int(i)] for i in bc.index], edgecolor='black', linewidth=0.3)
for bb, v in zip(bars, bc.values):
    ax[0].text(bb.get_x() + bb.get_width()/2, v + 15, f'{v:,}', ha='center', va='bottom', fontsize=6)
ax[0].set_xlabel('Number of co-occurring pathologies'); ax[0].set_ylabel('Number of participants')
ax[0].set_ylim(0, bc.max()*1.12)
cc = comp.path_combo.value_counts().head(12)
ax[1].barh(range(len(cc)), cc.values, color='#2166AC', edgecolor='black', linewidth=0.3, height=0.7)
ax[1].set_yticks(range(len(cc))); ax[1].set_yticklabels([LAB(s) for s in cc.index], fontsize=6.5); ax[1].invert_yaxis()
for i, v in enumerate(cc.values):
    ax[1].text(v + 6, i, f'{v:,}', va='center', fontsize=6)
ax[1].set_xlabel('Number of participants'); ax[1].set_xlim(0, cc.max()*1.14)
for a, l in zip(ax, 'ab'):
    a.text(-0.14, 1.04, l, transform=a.transAxes, fontsize=11, fontweight='bold', va='bottom')
fig.tight_layout(); save(fig, 'Figure_1', 240)

# Figure 2
fu = P.groupby('NACCID').years.max(); keep = fu[fu >= 6].index
S = P[P.NACCID.isin(keep)].assign(yb=lambda x: x.years.round().astype(int)); S = S[S.yb <= 6]
S['g'] = S.copath_count.astype(int)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.6))
for g in range(5):
    sub = S[S.g == g]; n = sub.NACCID.nunique()
    if n < 50:
        continue
    m = sub.groupby('yb').CDRSUM.agg(['mean', 'sem'])
    ax[0].plot(m.index, m['mean'], '-o', color=BC[g], ms=3.5, label=f'{g} patholog{"y" if g == 1 else "ies"} (n={n:,})', zorder=3)
    ax[0].fill_between(m.index, m['mean'] - 1.96*m['sem'], m['mean'] + 1.96*m['sem'], color=BC[g], alpha=0.18, lw=0)
ax[0].set_xlabel('Years from baseline'); ax[0].set_ylabel('CDR Sum of Boxes'); ax[0].set_xticks(range(7))
ax[0].legend(frameon=False, loc='upper left', fontsize=6.5)
ax[0].set_ylim(0, S.groupby(['g', 'yb']).CDRSUM.mean().max()*1.45)
rates['g'] = rates.g.astype(int)
bp = ax[1].boxplot([rates.loc[rates.g == n, 'rate'].dropna() for n in range(5)], tick_labels=[str(n) for n in range(5)],
                   patch_artist=True, widths=0.6, showfliers=False, medianprops=dict(color='black', linewidth=1),
                   boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5))
for i, pt in enumerate(bp['boxes']):
    pt.set_facecolor(BC[i]); pt.set_alpha(0.7)
ax[1].axhline(0, color='grey', ls='--', lw=0.5)
ax[1].set_xlabel('Number of co-occurring pathologies'); ax[1].set_ylabel('Annual CDR-SB change (points/year)')
sp = R['spearman']
ax[1].text(0.97, 0.97, f"ρ = {sp['rho']:.3f}\nP < 0.001", transform=ax[1].transAxes, ha='right', va='top', fontsize=7.5,
           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='grey', linewidth=0.5))
for a, l in zip(ax, 'ab'):
    a.text(-0.14, 1.04, l, transform=a.transAxes, fontsize=11, fontweight='bold', va='bottom')
fig.tight_layout(); save(fig, 'Figure_2', 240)
json.dump(dict(n=int(S.NACCID.nunique()), by_group={int(g): int(S[S.g == g].NACCID.nunique()) for g in range(5)}),
          open(HERE/'tables'/'figure2_fixed_composition.json', 'w'))

# Figure 3 (original order: TDP-43, ADNC, LBD, FTLD-tau)
items = [('years:p_TDP', 'TDP-43 proteinopathy'), ('years:p_ADNC', 'AD neuropathologic change'),
         ('years:p_LBD', 'Lewy body disease'), ('years:p_FTLDtau', 'FTLD-tauopathy')]
fig, ax = plt.subplots(figsize=(5.0, 3.0))
ax.axvline(0, color='grey', ls='--', lw=0.7, zorder=1)
his = []
for i, (k, lab) in enumerate(items):
    t = R['m3'][k]; y = len(items) - 1 - i; his.append(t['hi'])
    col = '#4878d0' if t['lo'] > 0 else '#999999'
    ax.errorbar(t['beta'], y, xerr=[[t['beta'] - t['lo']], [t['hi'] - t['beta']]], fmt='o', color=col, ms=6, capsize=4,
                lw=1.2, mec='black', mew=0.4, zorder=3)
    ax.text(t['hi'] + 0.02, y, f"{t['beta']:.2f} ({t['lo']:.2f}, {t['hi']:.2f})", va='center', fontsize=7.5)
ax.set_yticks(range(len(items))); ax.set_yticklabels([l for _, l in items][::-1], fontsize=8.5)
ax.set_xlabel('Additional CDR-SB points per year\n(interaction with time)')
ax.set_xlim(min(R['m3'][k]['lo'] for k, _ in items) - 0.05, max(his) + 0.32)
fig.tight_layout(); save(fig, 'Figure_3', 300)

# Figure 4
order = ['ADNC', 'LBD', 'TDP-43', 'FTLD-tau']
full = ['AD neuropathologic\nchange', 'Lewy body\ndisease', 'TDP-43\nproteinopathy', 'FTLD-\ntauopathy']
H = pd.DataFrame(R['heatmap'])
M = H.pivot(index='row', columns='col', values='mean').loc[order, order]
N = H.pivot(index='row', columns='col', values='n').loc[order, order]
M = M.where(N >= 10)
ann = np.array([[f'{M.iloc[i, j]:.2f}\n(n={int(N.iloc[i, j]):,})' if pd.notna(M.iloc[i, j]) else f'n={int(N.iloc[i, j])}'
                 for j in range(4)] for i in range(4)])
fig, ax = plt.subplots(figsize=(5.0, 4.0))
sns.heatmap(M.values, annot=ann, fmt='', cmap='YlOrRd', ax=ax, linewidths=0.8, linecolor='white',
            xticklabels=full, yticklabels=full, cbar_kws={'label': 'Annual CDR-SB change (points/year)', 'shrink': 0.8},
            annot_kws={'size': 8})
for k in range(4):
    ax.add_patch(plt.Rectangle((k, k), 1, 1, fill=False, edgecolor='black', lw=1.5))
ax.tick_params(labelsize=8)
fig.tight_layout(); save(fig, 'Figure_4', 300)
print('done')
