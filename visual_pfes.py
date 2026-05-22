import argparse
import os, re
import pandas as pd
import numpy as np
import shutil
import gzip
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors


import warnings


parser = argparse.ArgumentParser(description="Analyse PFES")
parser.add_argument('-l', '--log', type=str, help='log file name', default='progress.log') 
parser.add_argument('-s', '--pdbdir', type=str, help='directory with pdb files', default='structures')
parser.add_argument('-t', '--traj', type=str, help='make backbone trajectory', default='pfestraj.pdb')
parser.add_argument('-o', '--outdir', type=str, help='output directory name', default='visual_pfes_results')
parser.add_argument('-b', '--start', type=int, help='first point to read from trajectory', default=0)
parser.add_argument('-e', '--end', type=int, help='last point to read from trajectory', default=99999999)
parser.add_argument('--notraj', action='store_false', )
parser.add_argument('--noplots', action='store_false', )


args = parser.parse_args()


#class VisualPFES():


def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(data, key=alphanum_key)

def extract_lineage(log):
    traj_len = len(log)
    #pop_size = len(log[log.gndx == 'gndx0'])

    print(f'Processing a trajectory with {traj_len} mutations')
    lineage = log.drop_duplicates('gndx').tail(1)
    df = lineage
    ndx = df.id.to_string(index=False).strip()
        
    def return_ancestor(log, node):
        parent = log[log.id == node]
        parent = parent.drop_duplicates('sequence')
        return parent

    
    pbar = tqdm(desc='Extracting lineage')
    while not df.empty:
        ndx = return_ancestor(log, ndx)
        df = ndx
        lineage = pd.concat([lineage, df], axis=0)
        ndx = ndx.prev_id.to_string(index=False).strip()
        pbar.update(1)
    pbar.close()
    lineage = lineage.sort_index()
    r = lineage.tail(1).iloc[-1]
    W = 60
    print(f"\n  {'─'*W}")
    print(f"  Best sequence in lineage  (gen {r.gndx})")
    print(f"  {'─'*W}")
    print(f"  Score:    {float(r.score):.4f}   │  Length: {int(r.seq_len)}")
    print(f"  pLDDT:    {float(r.mean_plddt):.4f}   │  pTM:    {float(r.ptm):.4f}")
    if float(r.num_inter_conts) > 1:
        print(f"  iPLDDT:   {float(r.iplddt):.4f}   │  InterConts: {int(r.num_inter_conts)}")
    _amp_col = 'amp_prob' if 'amp_prob' in lineage.columns else ('s_amp' if 's_amp' in lineage.columns else None)
    if _amp_col:
        hemo_str = f"   │  hemo: {float(r.hemo_prob):.4f}" if 'hemo_prob' in lineage.columns else ""
        print(f"  AMP prob: {float(r[_amp_col]):.4f}{hemo_str}")
    print(f"  Contacts: {int(r.num_conts)}")
    seq = str(r.sequence)
    ss  = str(r.ss)
    print(f"  {'─'*W}")
    for i in range(0, max(len(seq), len(ss)), 60):
        print(f"  seq  {seq[i:i+60]}")
        print(f"  ss   {ss[i:i+60]}")
    print(f"  {'─'*W}\n")
    return lineage

#======================= make separate plots =======================#
def make_plots(log, bestlog, lineage):

    ms=0.1
    lw=1.4
    dpi=500

    os.makedirs(plotdir, exist_ok=True)
    for colname in log.select_dtypes(include=[np.number]).columns: 
        if not colname in ['seq', 'sequence', 'ss', 'genindex','dssp', 'mutation', 'index', 'id', 'prev_id', 'gndx', 'sel_mode']:
                fig, ax1 = plt.subplots(figsize=(9, 3))
                ax1.plot(log[colname].astype(float),'.', markersize=ms,    color='silver', label='all mutations')
                ax1.plot(bestlog[colname],'-', linewidth=lw, label='best of the generation')
                ax1.plot(lineage[colname],'-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={len(lineage[colname])})')
                ax1.legend(loc ="lower right")
                ax1.grid(True, which="both",linestyle='--', linewidth=0.3)
                ax1.set(xlabel="Total number of mutations", ylabel=colname.capitalize())
                #ax2 = ax1.twiny()
                #ax2.plot(lineage[colname].tolist(),'-', linewidth=lw, color='mediumslateblue')
                #ax2.set(xlabel="Lineage length")
                fig.tight_layout()
                fig.savefig(os.path.join(plotdir, colname + '.png'), dpi=dpi)
                fig.clf()

#======================= Summary plot =======================#
def make_summary_plot(log, bestlog, lineage):

    ms  = 0.1
    lw  = 1.0
    dpi = 500

    n_muts = len(log)
    n_gens = log.gndx.nunique()
    L      = len(lineage)

    # Rolling-mean window: ~5 % of total mutations, at least 10 points
    win = max(10, n_muts // 20)

    fig, axs = plt.subplots(3, 2, figsize=(10, 8))
    fig.suptitle(
        f"PFES run  │  {n_muts:,} sequences evaluated  │  "
        f"{n_gens} generations  │  lineage length {L}",
        fontsize=9)

    def _panel(ax, col, ylabel, xlabel=None, hide_x=True):
        v_all  = log[col].astype(float)
        v_roll = v_all.rolling(win, min_periods=1).mean()
        v_best = bestlog[col].astype(float)
        v_lin  = lineage[col].astype(float)
        ax.plot(v_all,  '.', markersize=ms, color='silver',         label='all evaluated')
        ax.plot(v_roll, '-', linewidth=0.9, color='steelblue', alpha=0.8,
                label=f'population rolling mean (w={win})')
        ax.plot(v_best, '-', linewidth=lw,  color='darkorange',     label='best of generation')
        ax.plot(v_lin,  '-', linewidth=lw,  color='mediumslateblue',label=f'lineage (L={L})')
        ax.set(ylabel=ylabel, xlabel=xlabel)
        ax.grid(True, which='both', linestyle='--', linewidth=0.4)
        if hide_x:
            ax.set_xticklabels([])

    _xlabel = 'Evaluated sequences (cumulative)'

    _panel(axs[0, 0], 'mean_plddt',
           'mean pLDDT\n(ESM3 per-residue confidence, 0–1)')
    _panel(axs[1, 0], 'ptm',
           'pTM\n(ESM3 predicted TM-score proxy, 0–1)')
    _panel(axs[2, 0], 'score',
           'Fitness score\n(product of all scoring terms)',
           xlabel=_xlabel, hide_x=False)
    _panel(axs[0, 1], 'seq_len',
           'Sequence length\n(amino acids)')

    # Adaptive middle-right panel
    if 'num_inter_conts' in bestlog.columns and bestlog.num_inter_conts.max() != 1:
        _panel(axs[1, 1], 'iplddt',
               'iPLDDT\n(interface pLDDT, 0–1)')
    elif 'hemo_prob' in log.columns:
        _panel(axs[1, 1], 'hemo_prob',
               'Hemolytic probability\n(MACREL, lower = safer, 0–1)')
    else:
        _panel(axs[1, 1], 'max_beta_penalty',
               'β-sheet penalty\n(sigmoid suppression, 0–1)')

    # Bottom-right panel
    _amp_col = 'amp_prob' if 'amp_prob' in log.columns else ('s_amp' if 's_amp' in log.columns else None)
    if _amp_col:
        ylabel = ('AMP probability\n(MACREL ML classifier, 0–1)'
                  if 'hemo_prob' in log.columns
                  else 'AMP score\n(biophysical s_amp, 0–1)')
        _panel(axs[2, 1], _amp_col, ylabel, xlabel=_xlabel, hide_x=False)
    else:
        _panel(axs[2, 1], 'num_conts',
               'Intra-chain contacts\n(Cα–Cα ≤ 8 Å)',
               xlabel=_xlabel, hide_x=False)

    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.02),
               ncol=4, fontsize=8)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.1)
    fig.savefig(os.path.join(outdir, 'Summary.png'), dpi=dpi, bbox_inches='tight')


#======================= Evolution rate plot =======================#
def make_evolution_plot(log, bestlog, lineage):
    """
    Two-panel figure:
      Top   — per-generation population mean ± 1 std (shaded), best score,
               and lineage score per generation.
      Bottom — evolution rate: Δ(best score) per generation as bars,
               with a rolling mean overlaid.
    """
    dpi = 500

    # Extract numeric generation index from gndx string (e.g. 'gndx12' → 12)
    log2 = log.copy()
    log2['_gen'] = log2['gndx'].str.extract(r'(\d+)', expand=False).astype(int)
    gen_stats = log2.groupby('_gen')['score'].agg(
        mean='mean', std='std', best='max'
    ).reset_index()
    gen_stats['std'] = gen_stats['std'].fillna(0)
    g = gen_stats['_gen'].values

    # Lineage score per generation (max score reached in each generation along the lineage)
    lin2 = lineage.copy()
    lin2['_gen'] = lin2['gndx'].str.extract(r'(\d+)', expand=False).astype(int)
    lin_gen = lin2.groupby('_gen')['score'].max().reset_index()

    # Rolling window for the rate panel
    rate_win = max(5, len(g) // 20)

    fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(
        f"PFES evolution rate  │  {log2['_gen'].nunique()} generations  │  "
        f"pop mean fitness and improvement per generation",
        fontsize=9)

    # --- Panel 1: population fitness per generation ---
    ax = axs[0]
    ax.fill_between(g,
                    gen_stats['mean'] - gen_stats['std'],
                    gen_stats['mean'] + gen_stats['std'],
                    alpha=0.20, color='steelblue',
                    label='population mean ± 1 std')
    ax.plot(g, gen_stats['mean'], '-', color='steelblue',       linewidth=1.2,
            label='population mean score')
    ax.plot(g, gen_stats['best'], '-', color='darkorange',      linewidth=1.2,
            label='best score of generation')
    ax.plot(lin_gen['_gen'], lin_gen['score'],
            '-', color='mediumslateblue', linewidth=1.0,
            label=f'lineage score (L={len(lineage)})')
    ax.set(ylabel='Fitness score\n(per generation)')
    ax.grid(True, linestyle='--', linewidth=0.4)
    ax.legend(fontsize=8, loc='upper left')

    # --- Panel 2: evolution rate ---
    ax2 = axs[1]
    delta = pd.Series(gen_stats['best']).diff().fillna(0).values
    bar_colors = ['steelblue' if v >= 0 else 'tomato' for v in delta]
    ax2.bar(g, delta, color=bar_colors, alpha=0.55, width=0.8,
            label='Δ best score / generation')
    ax2.axhline(0, color='black', linewidth=0.5)
    roll_rate = pd.Series(delta).rolling(rate_win, min_periods=1).mean()
    ax2.plot(g, roll_rate.values, '-', color='darkorange', linewidth=1.5,
             label=f'rolling mean ({rate_win} gen)')
    ax2.set(xlabel='Generation',
            ylabel='Δ best score / generation\n(positive = improvement)')
    ax2.grid(True, linestyle='--', linewidth=0.4)
    ax2.legend(fontsize=8, loc='upper left')

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'Evolution.png'), dpi=dpi, bbox_inches='tight')


#======================= Score component breakdown =======================#
def make_score_components_plot(log, lineage):
    """
    Multi-line chart of every individual scoring term along the best lineage,
    with the combined score overlaid in black.  Immediately shows which term
    is the rate-limiting factor at each evolutionary stage.
    """
    lin = lineage.reset_index(drop=True)
    x   = np.arange(len(lin))          # lineage step index

    # Build ordered dict: display_name -> Series of values (all in [0,1])
    terms = {}
    if 'mean_plddt' in lin.columns:
        terms['pLDDT  (ESM3 fold confidence)'] = lin['mean_plddt'].astype(float)
    if 'ptm' in lin.columns:
        terms['pTM  (ESM3 TM-score proxy)']    = lin['ptm'].astype(float)
    if 'prot_len_penalty' in lin.columns:
        terms['length penalty']                 = lin['prot_len_penalty'].astype(float)
    if 'max_alpha_penalty' in lin.columns:
        terms['helix penalty']                  = lin['max_alpha_penalty'].astype(float)
    if 'max_beta_penalty' in lin.columns:
        terms['β-sheet penalty']                = lin['max_beta_penalty'].astype(float)
    _amp_col = 'amp_prob' if 'amp_prob' in lin.columns else ('s_amp' if 's_amp' in lin.columns else None)
    if _amp_col:
        lbl = 'AMP probability  (MACREL)' if 'hemo_prob' in lin.columns else 'AMP score  (s_amp)'
        terms[lbl]                              = lin[_amp_col].astype(float)
    if 'hemo_prob' in lin.columns:
        terms['(1 − hemolytic prob)']           = (1 - lin['hemo_prob'].astype(float))
    if 'num_conts' in lin.columns and 'seq_len' in lin.columns:
        ct = (lin['num_conts'].astype(float) + lin['seq_len'].astype(float)) / lin['seq_len'].astype(float)
        terms['contact term  ((N_conts + L) / L)'] = ct.clip(upper=1.0)

    if not terms:
        return

    dpi = 500
    palette = plt.cm.tab10.colors
    fig, ax = plt.subplots(figsize=(11, 4))

    for idx, (name, vals) in enumerate(terms.items()):
        ax.plot(x, vals.values, '-', linewidth=1.4,
                color=palette[idx % len(palette)], label=name)

    ax.plot(x, lin['score'].astype(float).values,
            'k-', linewidth=2.0, label='combined score  (product)', zorder=10)

    ax.set(xlabel='Lineage step  (best ancestral sequence at each generation)',
           ylabel='Term value  [0 – 1]',
           title='Score component breakdown along the lineage\n'
                 'combined score = product of all displayed terms')
    ax.set_ylim(0, 1.08)
    ax.grid(True, linestyle='--', linewidth=0.4)
    ax.legend(fontsize=8, loc='upper left',
              bbox_to_anchor=(1.01, 1), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'Score_components.png'), dpi=dpi, bbox_inches='tight')


#======================= Fitness landscape =======================#
def make_fitness_landscape_plot(log, lineage):
    """
    2D scatter of structural quality (pLDDT) vs AMP fitness (AMP prob / contacts),
    every sequence ever evaluated coloured by generation number.
    The lineage path is overlaid in red with start (green) / end (red star) markers.
    Shows the evolutionary trajectory through fitness space.
    """
    _amp_col = 'amp_prob' if 'amp_prob' in log.columns else ('s_amp' if 's_amp' in log.columns else None)
    has_amp = _amp_col is not None
    if has_amp:
        x_col  = _amp_col
        xlabel = ('AMP probability  (MACREL ML classifier, 0–1)'
                  if 'hemo_prob' in log.columns else 'AMP score  (biophysical s_amp, 0–1)')
    else:
        x_col  = 'num_conts'
        xlabel = 'Intra-chain contacts  (Cα–Cα ≤ 8 Å)'

    log2 = log.copy()
    log2['_gen'] = log2['gndx'].str.extract(r'(\d+)', expand=False).astype(int)

    dpi = 500
    fig, ax = plt.subplots(figsize=(7, 5))

    sc = ax.scatter(
        log2[x_col].astype(float), log2['mean_plddt'].astype(float),
        c=log2['_gen'], cmap='viridis', s=2, alpha=0.35, linewidths=0,
        label='all evaluated sequences')
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.04)
    cb.set_label('Generation', fontsize=8)

    # Lineage path
    lin_x = lineage[x_col].astype(float).values
    lin_y = lineage['mean_plddt'].astype(float).values
    ax.plot(lin_x, lin_y, 'w-', linewidth=2.5, zorder=4)
    ax.plot(lin_x, lin_y, '-', color='tomato', linewidth=1.5, zorder=5,
            label=f'lineage path  (L={len(lineage)})')
    ax.scatter(lin_x[0],  lin_y[0],  color='limegreen', s=80,  zorder=7,
               marker='o', edgecolors='white', linewidths=0.8, label='lineage start')
    ax.scatter(lin_x[-1], lin_y[-1], color='tomato',    s=120, zorder=7,
               marker='*', edgecolors='white', linewidths=0.8, label='lineage end')

    ax.set(xlabel=xlabel,
           ylabel='mean pLDDT  (ESM3 fold confidence, 0–1)',
           title='Evolutionary trajectory in fitness space\n'
                 '(colour = generation number; each point = one evaluated sequence)')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, linestyle='--', linewidth=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'Fitness_landscape.png'), dpi=dpi, bbox_inches='tight')


#======================= Amino acid composition =======================#
def make_aa_composition_plot(lineage):
    """
    Heatmap of amino acid fractions along the best lineage.
    Shows how composition drifts over evolutionary time — for AMPs,
    expect R/K to rise (positive charge) and hydrophobic residues to stabilise.
    """
    AAs  = list('ACDEFGHIKLMNPQRSTVWY')
    seqs = lineage['sequence'].tolist()
    L    = len(seqs)
    mat  = np.zeros((20, L))
    for j, seq in enumerate(seqs):
        seq = str(seq)
        if len(seq) == 0:
            continue
        for i, aa in enumerate(AAs):
            mat[i, j] = seq.count(aa) / len(seq)

    dpi = 500
    fig, ax = plt.subplots(figsize=(11, 4))
    im = ax.imshow(mat, aspect='auto', origin='lower', cmap='YlOrRd',
                   interpolation='nearest', vmin=0, vmax=0.35)
    ax.set_yticks(range(20))
    ax.set_yticklabels(AAs, fontsize=8)
    ax.set(xlabel='Lineage step  (best ancestral sequence at each generation)',
           ylabel='Amino acid',
           title='Amino acid composition along the lineage\n'
                 '(fraction of each residue type in the sequence at each step)')
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cb.set_label('Fraction', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'AA_composition.png'), dpi=dpi, bbox_inches='tight')


#======================= Score distribution per generation =======================#
def make_score_distribution_plot(log):
    """
    Box plots of fitness score distributions per generation.
    Shows how selection pressure narrows the population over time and
    whether the whole population is shifting upward or just outliers.
    Sampled to ≤40 generations for readability.
    """
    log2 = log.copy()
    log2['_gen'] = log2['gndx'].str.extract(r'(\d+)', expand=False).astype(int)

    gen_ids = sorted(log2['_gen'].unique())
    max_show = 40
    if len(gen_ids) > max_show:
        step = max(1, len(gen_ids) // max_show)
        gen_ids = gen_ids[::step]

    data     = [log2[log2['_gen'] == g]['score'].astype(float).values for g in gen_ids]
    box_w    = max(0.6, gen_ids[-1] / len(gen_ids) * 0.75) if len(gen_ids) > 1 else 1

    dpi = 500
    fig, ax = plt.subplots(figsize=(12, 4))
    bp = ax.boxplot(
        data, positions=gen_ids, widths=box_w,
        patch_artist=True, showfliers=False,
        boxprops=dict(facecolor='steelblue', alpha=0.45, linewidth=0.6),
        medianprops=dict(color='darkorange', linewidth=1.8),
        whiskerprops=dict(linewidth=0.7),
        capprops=dict(linewidth=0.7))
    ax.set(xlabel='Generation',
           ylabel='Fitness score  (all evaluated sequences in generation)',
           title='Score distribution per generation\n'
                 '(box = IQR, orange bar = median, whiskers = 5–95th percentile)')
    ax.grid(True, linestyle='--', linewidth=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'Score_distribution.png'), dpi=dpi, bbox_inches='tight')


#======================= seconday structure plot =======================#
def make_ss_plot(lineage):

    dpi=500

    max_seq_len = int(max(lineage.seq_len))
    lineage_len = len(lineage)

    sse = np.empty((lineage_len, max_seq_len), dtype='U1')
    i=0
    for ss in lineage.ss:
        if not isinstance(ss, str):
            ss = ""
        sse[i] = list(ss + "X"*(max_seq_len-len(ss)))
        i+=1

    def sse_to_num(sse):
        num = np.empty(sse.shape, dtype=int)
        num[sse == 'F'] = 0
        num[sse == 'f'] = 0
        num[sse == 'g'] = 0
        num[sse == 's'] = 0
        num[sse == 'P'] = 0 
        num[sse == 'C'] = 0 
        num[sse == 'E'] = 1 
        num[sse == 'B'] = 2 
        num[sse == 'S'] = 3 
        num[sse == 'T'] = 4 
        num[sse == 'H'] = 5 
        num[sse == 'G'] = 6 
        num[sse == 'I'] = 7 
        num[sse == 'X'] = 8 
        return num


    sse_digit = sse_to_num(sse)

    color_assign = {
        r"coil": "darkgrey",
        r"$\beta$-sheet": "yellow",
        r"$\beta$-bridge": "y",
        r"bend": "orange",
        r"turn": "brown",
        r"$\alpha$-helix": "purple",
        r"$3_{10}$-helix": "blue",
        r"$\pi$-helix": "mediumpurple",
        r"": "white"
        }

    cmap = colors.ListedColormap(color_assign.values())
    if len(lineage) > 2000:
        ticks = np.arange(0, len(lineage)+1, 1000)
    else: 
        ticks = np.arange(0, len(lineage)+1, 100)

    plt.figure(figsize=(9, 3), dpi=dpi)
    plt.imshow(sse_digit.T, origin='lower', cmap=cmap,  interpolation='nearest', aspect='auto')
    plt.xticks(ticks, ticks.astype(int))
    plt.xlabel("Lineage length")
    plt.ylabel("Residues")

    custom_lines = [
        Line2D([0], [0], color=cmap(i), lw=4) for i in range(len(color_assign)-1)]

    plt.legend(
        custom_lines, color_assign.keys(), loc="upper center",
        bbox_to_anchor=(0.5, 1.15), ncol=len(color_assign), fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'Secondary_structures.png'), dpi=dpi) 

def backbone_traj(trajlog, pdbdir):
    """
    make trajectory from C-alpha atoms
    """
    trajpdb = os.path.join(outdir, 'trajpdb/')
    trajlog = trajlog.drop_duplicates(subset = 'sequence')
    pfeslen = len(trajlog)

    os.makedirs(trajpdb, exist_ok=True)
    if os.path.isdir(trajpdb) and len(os.listdir(trajpdb)) != pfeslen:
        shutil.rmtree(trajpdb)
        os.makedirs(trajpdb, exist_ok=True)
        print(f"  Copying {pfeslen} unique best-fold structures...")
        for gndx, pdbid in tqdm(zip(trajlog.gndx, trajlog.id), total=len(trajlog)):
            try:
                shutil.copy(pdbdir +'/' + pdbid +'.pdb.gz', trajpdb +'/'+ gndx + '.pdb.gz')
            except FileNotFoundError:
                print(pdbid +'.pdb.gz is missing' )
                pass
    else:
        print(f"  {pfeslen} best-fold structures already present, skipping copy")

    print("  Extracting backbone coordinates...")
    i=0
    PDB_A, PDB_B, lastBB_A, lastBB_B = [], [], [], []
    for pdb in tqdm(sorted_alphanumeric(os.listdir(trajpdb))):
        with gzip.open(os.path.join(trajpdb, pdb), 'rb') as file:
            pdb_txt = file.read().decode()
        bb_chain_A, bb_chain_B = [], []
        for line in pdb_txt.splitlines():
            col = line.split()
            if (col[0] == 'ATOM' and (
                col[2] == 'N' or 
                col[2] == 'CA' or 
                col[2] == 'C' or 
                col[2] == 'O') and 
                col[4] == 'A'):
                bb_chain_A.append(line + '\n') 
            if (col[0] == 'ATOM' and (
                col[2] == 'N' or 
                col[2] == 'CA' or 
                col[2] == 'C' or 
                col[2] == 'O') and 
                col[4] == 'B'):
                bb_chain_B.append(line + '\n')

        lastresidueA=''.join([str(elem) for elem in bb_chain_A[-4:]]) # keep the last four lines to repeat them and make the number of atoms in all models equal
        lastresidueB=''.join([str(elem) for elem in bb_chain_B[-4:]]) # keep the last four lines to repeat them and make number of atoms in all models equal

        PDB_A.append(bb_chain_A) #save chain A
        PDB_B.append(bb_chain_B) #save chain B
        lastBB_A.append(lastresidueA) #save bb of the last residue of chain A 
        lastBB_B.append(lastresidueB) #save bb of the last residue of chain B

    topmax_A = max([len(i) for i in PDB_A])
    topmax_B = max([len(i) for i in PDB_B])
    for pdbA, pdbB in zip(PDB_A, PDB_B): 
        if len(pdbA) == topmax_A:
            toppdb = ''.join(pdbA + pdbB)
            break
    
    print("  Preparing backbone trajectory...")
    with open(outdir+'/.tmp.pdb', 'w') as f:
        i=1
        f.write(f'MODEL        {i}\n' + toppdb + 'TER\nENDMDL\n')
        for chA, chB, lastresBB_A, lastresBB_B in tqdm(zip(PDB_A, PDB_B, lastBB_A, lastBB_B), total=len(PDB_A)): 
            i+=1
            dumlinesA = lastresBB_A  * int((topmax_A - len(chA)) / 4)
            dumlinesB = lastresBB_B  * int((topmax_B - len(chB)) / 4)

            f.write(f'MODEL        {i}\n' + ''.join(chA) + dumlinesA + 'TER\n' + ''.join(chB) + dumlinesB + 'ENDMDL\n')
    try: 
        import MDAnalysis as mda
        from MDAnalysis.analysis import align

        print("  Writing aligned backbone trajectory...")
        traj = mda.Universe(outdir+'/.tmp.pdb')
        top = traj.select_atoms('protein')

        warnings.filterwarnings("ignore")
        if 'num_inter_conts' in trajlog.columns and trajlog.num_inter_conts.max() != 1:
            chianID = 'chainID B'
        else:
            chianID = 'chainID A'

        align.AlignTraj(traj,  # trajectory to align
                        top,  # reference
                        select=chianID,  # selection of atoms to align
                        filename=trajpath,  # file to write the trajectory to
                        ).run()

        os.remove(outdir+'/.tmp.pdb')
    except:
        print("  Warning: MDAnalysis not available — writing unaligned trajectory")
        shutil.copy(outdir+'/.tmp.pdb', trajpath)



outdir = args.outdir 
os.makedirs(outdir, exist_ok=True)

pdbdir = args.pdbdir
plotdir = os.path.join(outdir, 'plots/')
trajpath = os.path.join(outdir, args.traj)

try:
    log = pd.read_csv(args.log, sep='\t', comment='#')
except FileNotFoundError:
    print(f'Error: Could not find log file at {args.log}')
    exit(1)

log = log.iloc[args.start:args.end]

bestlog = log.groupby('gndx').head(1)
bestlog.to_csv(os.path.join(outdir, 'bestlog.tsv'), sep='\t', index=False, header=True)


W = 60
print(f"\n  {'═'*W}")
print(f"  Extracting lineage...")
print(f"  {'─'*W}")
lineage = extract_lineage(log)
lineage.to_csv(os.path.join(outdir, 'lineage.tsv'), sep='\t', index=False, header=True)
print(f"  Lineage length: {len(lineage)}  │  saved → {os.path.join(outdir, 'lineage.tsv')}")

if args.noplots:
    print(f"\n  {'─'*W}")
    print(f"  Generating plots → {plotdir}")
    print(f"  {'─'*W}")
    make_plots(log, bestlog, lineage)
    print(f"  ✓ individual column plots")
    make_summary_plot(log, bestlog, lineage)
    print(f"  ✓ summary plot")
    make_evolution_plot(log, bestlog, lineage)
    print(f"  ✓ evolution rate plot")
    make_score_components_plot(log, lineage)
    print(f"  ✓ score components plot")
    make_fitness_landscape_plot(log, lineage)
    print(f"  ✓ fitness landscape plot")
    make_aa_composition_plot(lineage)
    print(f"  ✓ amino acid composition plot")
    make_score_distribution_plot(log)
    print(f"  ✓ score distribution plot")
    make_ss_plot(lineage)
    print(f"  ✓ secondary structure plot")

if args.notraj:
    print(f"\n  {'─'*W}")
    print(f"  Building backbone trajectory → {trajpath}")
    print(f"  {'─'*W}")
    backbone_traj(lineage, pdbdir)
    print(f"  ✓ trajectory written")

print(f"\n  {'═'*W}")
print(f"  Done.  Results in: {outdir}")
print(f"  {'═'*W}\n")
