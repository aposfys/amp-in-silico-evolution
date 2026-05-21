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
    if 's_amp' in lineage.columns:
        hemo_str = f"   │  hemo: {float(r.hemo_prob):.4f}" if 'hemo_prob' in lineage.columns else ""
        print(f"  AMP prob: {float(r.s_amp):.4f}{hemo_str}")
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
    if 's_amp' in log.columns:
        ylabel = ('AMP probability\n(MACREL ML classifier, 0–1)'
                  if 'hemo_prob' in log.columns
                  else 'AMP score\n(biophysical s_amp, 0–1)')
        _panel(axs[2, 1], 's_amp', ylabel, xlabel=_xlabel, hide_x=False)
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


#======================= seconday structure plot =======================#
def make_ss_plot(lineage):

    dpi=500

    max_seq_len = int(max(lineage.seq_len))
    lineage_len = len(lineage)

    sse = np.empty((lineage_len, max_seq_len), dtype='U1')
    i=0
    for ss in lineage.ss:
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
