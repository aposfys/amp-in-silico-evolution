"""
Generate 4 side-by-side comparison plots: MACREL vs MACREL+PFES
Saves to comparison_plots/
"""
import os, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
warnings.filterwarnings('ignore')

os.makedirs('comparison_plots', exist_ok=True)

CA = '#1565C0'    # MACREL blue
CB = '#E65100'    # MACREL+PFES orange
CA_L = '#90CAF9'
CB_L = '#FFCC80'
LABEL_A = 'MACREL'
LABEL_B = 'MACREL + PFES'

DPI = 180
FONT = {'family': 'DejaVu Sans', 'size': 11}
plt.rc('font', **FONT)


def load(path, amp_col=None):
    df = pd.read_csv(path, sep='\t', comment='#', low_memory=False)
    # normalise amp column
    if 'amp_prob' in df.columns:
        df['amp'] = pd.to_numeric(df['amp_prob'], errors='coerce')
    elif 's_amp' in df.columns:
        df['amp'] = pd.to_numeric(df['s_amp'], errors='coerce')
    for col in ['score', 'hemo_prob', 'mean_plddt', 'ptm', 'seq_len',
                'prot_len_penalty', 'max_alpha_penalty', 'max_beta_penalty']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    df['gen'] = df['gndx'].str.extract(r'(\d+)').astype(float)
    return df.dropna(subset=['score', 'gen'])


A = load('results_macrel_prod/progress.log')
B = load('results_macrel_pfes_pro/progress.log')

LAST_GEN = int(min(A.gen.max(), B.gen.max()))


def smooth(s, w=7):
    return s.rolling(w, min_periods=1, center=True).mean()


def gen_stat(df, col):
    g = df.groupby('gen')[col]
    return g.max(), g.mean(), g.quantile(0.25), g.quantile(0.75)


def legend_handles():
    return [
        Line2D([0], [0], color=CA, lw=2.5, label=LABEL_A),
        Line2D([0], [0], color=CB, lw=2.5, label=LABEL_B),
    ]


# ── 1. Fitness over generations ─────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
fig.suptitle('Fitness over Generations', fontsize=15, fontweight='bold', y=0.98)

for ax, df, color, light, label in [
    (axes[0], A, CA, CA_L, LABEL_A),
    (axes[1], B, CB, CB_L, LABEL_B),
]:
    mx, mn, q25, q75 = gen_stat(df, 'score')
    gens = mx.index
    ax.fill_between(gens, smooth(q25), smooth(q75), color=light, alpha=0.4, label='IQR')
    ax.plot(gens, smooth(mx), color=color, lw=2.5, label='Best score/gen')
    ax.plot(gens, smooth(mn), color=color, lw=1.2, ls='--', alpha=0.7, label='Mean score/gen')
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title(label, fontsize=12, color=color, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='lower right')

axes[1].set_xlabel('Generation', fontsize=11)
plt.tight_layout()
plt.savefig('comparison_plots/1_fitness_over_generations.png', dpi=DPI, bbox_inches='tight')
plt.close()


# ── 2. Score components ─────────────────────────────────────────────────────
components = [
    ('score',    'Overall Score',  'black'),
    ('amp',      'AMP probability','#2E7D32'),
    ('hemo_prob','Hemolytic prob', '#C62828'),
    ('mean_plddt','pLDDT',         '#1565C0'),
    ('ptm',      'pTM',            '#7B1FA2'),
]

fig, axes = plt.subplots(len(components), 2, figsize=(13, 11), sharex=True)
fig.suptitle('Score Components Over Generations', fontsize=15, fontweight='bold', y=0.99)

axes[0, 0].set_title(LABEL_A, fontsize=13, color=CA, fontweight='bold', pad=8)
axes[0, 1].set_title(LABEL_B, fontsize=13, color=CB, fontweight='bold', pad=8)

for row, (col, label, color) in enumerate(components):
    for ax_col, (df, run_color) in enumerate([(A, CA), (B, CB)]):
        ax = axes[row, ax_col]
        if col not in df.columns:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes)
            continue
        mx, mn, q25, q75 = gen_stat(df, col)
        gens = mx.index
        ax.fill_between(gens, smooth(q25), smooth(q75), color=color, alpha=0.15)
        ax.plot(gens, smooth(mx), color=color, lw=2.0)
        ax.plot(gens, smooth(mn), color=color, lw=1.0, ls='--', alpha=0.6)
        ax.set_ylabel(label, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.25)

for ax in axes[-1]:
    ax.set_xlabel('Generation', fontsize=10)

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig('comparison_plots/2_score_components.png', dpi=DPI, bbox_inches='tight')
plt.close()


# ── 3. Amino acid composition (final generation) ─────────────────────────────
AAS = list('ACDEFGHIKLMNPQRSTVWY')

def aa_freq(df, last_n_gens=10):
    final = df[df.gen >= df.gen.max() - last_n_gens]
    counts = {aa: 0 for aa in AAS}
    total = 0
    for seq in final.sequence.dropna():
        for aa in str(seq):
            if aa in counts:
                counts[aa] += 1
                total += 1
    return {aa: counts[aa]/total*100 for aa in AAS} if total else {aa: 0 for aa in AAS}

freq_a = aa_freq(A)
freq_b = aa_freq(B)

x = np.arange(len(AAS))
w = 0.38

fig, ax = plt.subplots(figsize=(13, 5))
bars_a = ax.bar(x - w/2, [freq_a[aa] for aa in AAS], w, color=CA, alpha=0.85, label=LABEL_A)
bars_b = ax.bar(x + w/2, [freq_b[aa] for aa in AAS], w, color=CB, alpha=0.85, label=LABEL_B)
ax.set_xticks(x)
ax.set_xticklabels(AAS, fontsize=11)
ax.set_ylabel('Frequency in final population (%)', fontsize=11)
ax.set_title('Amino Acid Composition — Final 10 Generations', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, axis='y', alpha=0.3)
ax.set_ylim(0, max(max(freq_a.values()), max(freq_b.values())) * 1.2)

# highlight top AA per run
top_a = sorted(AAS, key=lambda aa: freq_a[aa], reverse=True)[:3]
top_b = sorted(AAS, key=lambda aa: freq_b[aa], reverse=True)[:3]
for aa in top_a:
    i = AAS.index(aa)
    ax.annotate(f'{freq_a[aa]:.1f}%', (i - w/2, freq_a[aa]+0.1), ha='center', fontsize=7.5, color=CA, fontweight='bold')
for aa in top_b:
    i = AAS.index(aa)
    ax.annotate(f'{freq_b[aa]:.1f}%', (i + w/2, freq_b[aa]+0.1), ha='center', fontsize=7.5, color=CB, fontweight='bold')

plt.tight_layout()
plt.savefig('comparison_plots/3_aa_composition.png', dpi=DPI, bbox_inches='tight')
plt.close()


# ── 4. Score distribution (final generation) ─────────────────────────────────
final_a = A[A.gen == A.gen.max()]['score'].dropna()
final_b = B[B.gen == B.gen.max()]['score'].dropna()

fig, ax = plt.subplots(figsize=(10, 5))
bins = np.linspace(0, 1, 40)
ax.hist(final_a, bins=bins, color=CA, alpha=0.7, label=f'{LABEL_A}  (mean={final_a.mean():.3f})', density=True)
ax.hist(final_b, bins=bins, color=CB, alpha=0.7, label=f'{LABEL_B}  (mean={final_b.mean():.3f})', density=True)

ax.axvline(final_a.mean(), color=CA, lw=2, ls='--', alpha=0.9)
ax.axvline(final_b.mean(), color=CB, lw=2, ls='--', alpha=0.9)

ax.set_xlabel('Score', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('Score Distribution — Final Generation (gen 499)', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('comparison_plots/4_score_distribution.png', dpi=DPI, bbox_inches='tight')
plt.close()

# ── 5. Sequence length over generations ─────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4.5))
for df, color, label in [(A, CA, LABEL_A), (B, CB, LABEL_B)]:
    mx, mn, q25, q75 = gen_stat(df, 'seq_len')
    gens = mx.index
    ax.fill_between(gens, smooth(q25), smooth(q75), color=color, alpha=0.2)
    ax.plot(gens, smooth(mx), color=color, lw=2.5, label=f'{label} — best length/gen')
    ax.plot(gens, smooth(mn), color=color, lw=1.2, ls='--', alpha=0.7)
ax.set_xlabel('Generation', fontsize=11)
ax.set_ylabel('Sequence length (aa)', fontsize=11)
ax.set_title('Sequence Length Over Generations', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.axhline(30, color='grey', lw=1, ls=':', label='30 aa (typical AMP)')
plt.tight_layout()
plt.savefig('comparison_plots/5_sequence_length.png', dpi=DPI, bbox_inches='tight')
plt.close()

print("Plots saved to comparison_plots/")
for f in sorted(os.listdir('comparison_plots')):
    print(f"  {f}")
