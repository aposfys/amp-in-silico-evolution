"""Generic PFES loop diagram — concept-level, no tool names in boxes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig = plt.figure(figsize=(14, 8))
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 14); ax.set_ylim(0, 8); ax.axis('off')
fig.patch.set_facecolor('white')

# ── colours ────────────────────────────────────────────────────────────────
C = ["#1565C0","#2E7D32","#E65100","#6A1B9A"]  # blue green orange purple
L = ["#BBDEFB","#C8E6C9","#FFE0B2","#E1BEE7"]
DARK = "#263238"

# ── subtle background oval ──────────────────────────────────────────────────
ax.add_patch(mpatches.Ellipse((7,4), 13.2, 7.4,
    facecolor="#F8F9FA", edgecolor="#ECEFF1", lw=1.5, zorder=0))

# ── centre label ────────────────────────────────────────────────────────────
ax.text(7, 4, "In silico\nevolution", ha='center', va='center',
        fontsize=16, fontweight='bold', color='#455A64',
        bbox=dict(boxstyle='round,pad=0.6', facecolor='white',
                  edgecolor='#90A4AE', lw=1.5, zorder=5))

# ── nodes (top, right, bottom, left) ────────────────────────────────────────
nodes = [
    (0, "MUTATE",            "One random mutation\nper sequence",                    7,   6.9),
    (1, "PREDICT\nSTRUCTURE","3D structure per sequence\nConfidence scores (pLDDT, pTM)",12,  4.0),
    (2, "SCORE",             "Multiplicative fitness:\nall terms must pass",         7,   1.1),
    (3, "SELECT",            "Boltzmann sampling\nFavour high-scorers, keep diversity",2,   4.0),
]
BW, BH = 3.3, 1.15

for idx, label, sub, cx, cy in nodes:
    # shadow
    ax.add_patch(FancyBboxPatch((cx-BW/2+0.06, cy-BH/2-0.06), BW, BH,
        boxstyle="round,pad=0.14", facecolor='#B0BEC5', lw=0, zorder=2, alpha=0.4))
    # box
    ax.add_patch(FancyBboxPatch((cx-BW/2, cy-BH/2), BW, BH,
        boxstyle="round,pad=0.14", facecolor=L[idx], edgecolor=C[idx], lw=2.8, zorder=3))
    ax.text(cx, cy+0.08, label, ha='center', va='center',
            fontsize=16, fontweight='bold', color=C[idx], zorder=4)
    # sub-label
    ax.text(cx, cy-BH/2-0.32, sub, ha='center', va='top',
            fontsize=9.5, color=C[idx], style='italic', zorder=4,
            bbox=dict(boxstyle='round,pad=0.22', facecolor='white',
                      edgecolor=C[idx], lw=0.7, alpha=0.9))

# ── arrows (clockwise) ───────────────────────────────────────────────────────
akw = dict(arrowstyle='->', color='#455A64', lw=2.4,
           connectionstyle='arc3,rad=0.18', shrinkA=30, shrinkB=30)
pairs = [
    ((7+BW/2, 6.9), (12, 4+BH/2)),
    ((12, 4-BH/2),  (7+BW/2, 1.1)),
    ((7-BW/2, 1.1), (2, 4-BH/2)),
    ((2, 4+BH/2),   (7-BW/2, 6.9)),
]
for (x1,y1),(x2,y2) in pairs:
    ax.annotate("", xy=(x2,y2), xytext=(x1,y1),
                arrowprops=dict(**akw), zorder=3)

# ── flow state labels on arrows ──────────────────────────────────────────────
flow = [
    (10.1, 5.88, "100 new\ncandidates"),
    (10.3, 2.16, "100 structures\n+ metrics"),
    (3.9,  2.16, "200 sequences\nscored"),
    (3.7,  5.88, "100 survivors\n→ next gen"),
]
for fx,fy,txt in flow:
    ax.text(fx, fy, txt, ha='center', va='center', fontsize=10,
            color=DARK, bbox=dict(boxstyle='round,pad=0.3',
            facecolor='#ECEFF1', edgecolor='#B0BEC5', lw=0.9))

# ── top banner ───────────────────────────────────────────────────────────────
ax.text(7, 7.65,
        "300 generations  ·  population 100 sequences  ·  30,000 structure predictions",
        ha='center', va='center', fontsize=12, color=DARK,
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                  edgecolor='#CFD8DC', lw=1.2))

# ── bottom note ──────────────────────────────────────────────────────────────
ax.text(7, 0.22,
        "Structure predictor: ESM3  ·  AMP classifier: MACREL  ·  Selection: Boltzmann (β = 20)",
        ha='center', va='center', fontsize=10, color='#78909C', style='italic')

plt.savefig("loop_diagram.png", dpi=220, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved → loop_diagram.png")
