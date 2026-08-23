import argparse
import os
import sys
import shutil
import pandas as pd
import numpy as np
import typing as T
import threading
import time
from datetime import datetime

os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"
os.environ["TOKENIZERS_PARALLELISM"] = "false"
import torch
import esm

# Monkeypatch ESM3 fp32_autocast_context to support MPS
from contextlib import nullcontext
try:
    from esm.utils import misc
    original_fp32_autocast_context = misc.fp32_autocast_context
    def mps_compatible_fp32_autocast_context(device_type: str):
        if device_type == "mps":
            return nullcontext()
        return original_fp32_autocast_context(device_type)
    misc.fp32_autocast_context = mps_compatible_fp32_autocast_context
except ImportError:
    pass

from evolution import Evolver
from score import get_nconts, cbiplddt, calculate_samp, calculate_hemo_proxy, macrel_score_batch_src
from psique import pypsique

from esm.sdk.api import ESMProtein, GenerationConfig
from esm.models.esm3 import ESM3


def backup_output(outpath):
    print(f'\nSaving output files to {outpath}')
    if os.path.isdir(outpath):
        parent = os.path.dirname(outpath) or '.'
        base = os.path.basename(outpath)
        backup_list = []
        last_backup = 0
        for dir_name in os.listdir(parent):
            if dir_name.startswith(base + '.'):
                suffix = dir_name[len(base) + 1:]
                if suffix.isdigit():
                    backup_list.append(int(suffix))
                    last_backup = max(backup_list)
        print(f'\n{outpath} already exists, renaming it to {outpath}.{last_backup + 1}')
        os.replace(outpath, outpath + '.' + str(last_backup + 1))


def create_batched_sequence_dataset(sequences: T.List[T.Tuple[str, str]], max_tokens_per_batch: int = 1524
) -> T.Generator[T.Tuple[T.List[str], T.List[str]], None, None]:
    batch_headers, batch_sequences, num_tokens, num_sequences= [], [], 0, 0
    for header, seq in sequences:
        if (len(seq) + num_tokens > max_tokens_per_batch) and num_tokens > 0:
            yield batch_headers, batch_sequences
            batch_headers, batch_sequences, num_tokens, num_sequences= [], [], 0, 0
        batch_headers.append(header)
        batch_sequences.append(seq)
        num_tokens += len(seq)
        num_sequences += 1
        if num_sequences > args.pop_size / 2: #TODO test this with args.pop_size / 4 and lartge pop size
           yield batch_headers, batch_sequences
           batch_headers, batch_sequences, num_tokens, num_sequences= [], [], 0, 0
    if batch_headers:
        yield batch_headers, batch_sequences

def pdbtxt2bbcoord(pdb_txt, chain='A'):
    # can extract this directly from esm output
    # positions contains coordinates, and aatype contains the sequence
    coords3 = np.array([line[30:54].split()  for line in pdb_txt.splitlines() if line[:4] == "ATOM" and 
                        line[20:22].strip() == chain and 
                        ((line[11:16].strip() == "N") | 
                         (line[11:16].strip()== "CA") | 
                         (line[11:16].strip() == "C"))], dtype='float32')
    coords33 = coords3.reshape(int(coords3.shape[0]/3),3,3)
    return(coords33)

def esm2data(esm_out):
    # ESM3 outputs a list of ESMProtein objects
    pdbs = []
    ptms = []
    mean_plddts = []

    for prot in esm_out:
        # PDB text
        pdbs.append(prot.to_pdb_string())
        
        # pLDDT
        if prot.plddt is not None:
            if isinstance(prot.plddt, torch.Tensor):
                plddt_val = prot.plddt.mean().item()
            else:
                plddt_val = float(np.mean(prot.plddt))
            # Scale pLDDT to 0-1 if it is 0-100
            if plddt_val > 1.0:
                plddt_val /= 100.0
            mean_plddts.append(plddt_val)
        else:
            mean_plddts.append(0.0)

        # pTM
        if prot.ptm is not None:
            ptms.append(prot.ptm.item() if isinstance(prot.ptm, torch.Tensor) else float(prot.ptm))
        else:
            ptms.append(mean_plddts[-1])

    return(pdbs, ptms, mean_plddts) #return score instead

    #calculate the number of contacts
    # bins = np.append(0,np.linspace(2.3125,21.6875,63))
    # #you do not need softmax to keep the actual values 
    # sm_contacts = softmax(output["distogram_logits"],-1)
    # sm_contacts = sm_contacts[...,bins<8].sum(-1)
    # mask = output["atom37_atom_exists"][0,:,1] == 1
    # contact_map = sm_contacts[0][mask,:][:,mask]
    # num_conts = []
    """
    Return the number of contacts and individual plddts (write it in the log). 
    In the case of dimers, the number of interchain interactions with indexes is also returned. 
    Use indexes to calculate iPLDDT

    """

def sigmoid(x,L0=0,c=0.1):
    return 1 / (1+2.71828182**(c * (L0-x)))


_W = 74  # output column width

def _print_startup(args, evolver, date_now, time_now):
    print(f"\n{'═'*_W}")
    print(f"  PFES v0.1  │  {date_now}  {time_now}  │  {os.getcwd()}")
    print(f"{'─'*_W}")
    print(f"  mode:       {args.evolution_mode}  │  selection: {args.selection_mode}  (β={args.beta})")
    print(f"  pop size:   {args.pop_size}  │  generations: {args.num_generations}  │  evoldict: {args.evoldict}")
    init_display = args.initial_seq if len(args.initial_seq) <= 40 else args.initial_seq[:40] + '…'
    print(f"  init seq:   {init_display}  (rand_len={args.random_seq_len})")
    print(f"  penalties:  len≥{args.prot_len_penalty}  helix≥{args.helix_len_penalty}  beta≥{args.beta_len_penalty}")
    _hemo = "×(1-hemo)" if args.hemo_in_score else "hemo=attribute-only (not in score)"
    print(f"  scoring:    pLDDT×pTM×len_pen×helix_pen×beta_pen×AMP{'×(1-hemo)' if args.hemo_in_score else ''}×contacts  [MACREL AMP + PFES; {_hemo}]")
    print(f"  output:     {args.outpath}/{args.log}")
    print(f"{'═'*_W}\n")

def _print_gen_summary(gen_i, num_gen, new_gen, elapsed, best_so_far=0.0):
    top = new_gen.sort_values('score', ascending=False)
    n_total = len(top)
    this_best = float(top.iloc[0].score) if not top.empty else 0.0
    pop_mean  = float(top.score.mean()) if not top.empty else 0.0
    delta = this_best - best_so_far
    sign = '↑' if delta >= 0 else '↓'
    has_samp = 'amp_prob' in top.columns or 's_amp' in top.columns
    has_hemo = 'hemo_prob' in top.columns
    print(f"\n  ── Gen {gen_i + 1}/{num_gen}  ({elapsed:.1f}s)  "
          f"best={this_best:.4f} {sign}{abs(delta):.4f}  pop_mean={pop_mean:.4f}  "
          f"{'─' * max(2, _W - 55)}")
    hdr = f"  {'score':>6}  {'pLDDT':>6}  {'pTM':>5}"
    hdr += f"  {'AMP':>5}" if has_samp else ""
    hdr += f"  {'hemo':>5}" if has_hemo else ""
    hdr += f"  {'len':>4}  {'mutation':<14}  sequence"
    print(hdr)
    sep = f"  {'─'*6}  {'─'*6}  {'─'*5}"
    sep += f"  {'─'*5}" if has_samp else ""
    sep += f"  {'─'*5}" if has_hemo else ""
    sep += f"  {'─'*4}  {'─'*14}  {'─'*20}"
    print(sep)
    for _, row in top.head(10).iterrows():
        seq = str(row.sequence)
        seq_disp = seq[:20] + '…' if len(seq) > 20 else seq
        mut = str(row.mutation)[:14]
        line = f"  {float(row.score):>6.4f}  {float(row.mean_plddt):>6.3f}  {float(row.ptm):>5.3f}"
        amp_col = 'amp_prob' if 'amp_prob' in row.index else 's_amp'
        line += f"  {float(row[amp_col]):>5.3f}" if has_samp else ""
        line += f"  {float(row.hemo_prob):>5.3f}" if has_hemo else ""
        line += f"  {int(row.seq_len):>4}  {mut:<14}  {seq_disp}"
        print(line)
    if n_total > 10:
        print(f"  … {n_total - 10} more sequences not shown")
    print()
    return this_best


def _print_run_summary(log_path, total_time):
    try:
        df = pd.read_csv(log_path, sep='\t', comment='#').dropna(subset=['score', 'mean_plddt'])
        best = df.loc[df.score.idxmax()]
        n_gens = df.gndx.nunique()
        print(f"\n  {'═'*_W}")
        print(f"  Run complete  │  {n_gens} generations  │  {len(df):,} sequences evaluated  │  {total_time:.1f}s")
        print(f"  Best score:   {best.score:.4f}  │  pLDDT={best.mean_plddt:.3f}  pTM={best.ptm:.3f}", end="")
        amp_col = 'amp_prob' if 'amp_prob' in best.index else ('s_amp' if 's_amp' in best.index else None)
        if amp_col:
            print(f"  AMP={best[amp_col]:.3f}  hemo={best.hemo_prob:.3f}", end="")
        print(f"\n  Best sequence: {best.sequence}  (gen {best.gndx})")
        print(f"  {'═'*_W}\n")
    except Exception:
        pass


#==============================================================================================#
#================================== EXTRACT AND SCORE =========================================#
#==============================================================================================#

def extract_results(gen_i, headers, sequences, pdbs, ptms, mean_plddts, macrel_scores) -> None:
    global new_gen #this will be modified in the fold_evolver()

    for meta_id, seq, pdb_txt, ptm, mean_plddt, in zip(headers, sequences, pdbs, ptms, mean_plddts): #which plddt is better? this is plddt for both A and B chains in case of inter_chain

        all_seqs = seq.split(':')
        seq = all_seqs[0]
        seq_len = len(seq)

        # meta_id format: "{id}_{prev_id}_{mutation}"
        # prev_id may contain underscores (e.g. "init_seq"), mutation never does
        id, rest = meta_id.split('_', 1)
        prev_id, mutation = rest.rsplit('_', 1)

        with open(pdb_path + id + '.pdb', 'wb') as f:
            f.write(pdb_txt.encode())

        #=======================================================================#
        #================================SCORING================================#
        # Eq. 5 of Sahakyan et al.: Cbeta contacts within 6 A, |i-j| > 5,
        # both residues above the pLDDT floor. 6.0/0.5 rather than the
        # upstream 6.0/50 because ESM3 writes pLDDT on a 0-1 scale;
        # get_nconts normalises either scale, so both forms are safe.
        num_conts, _mean_plddt_ = get_nconts(pdb_txt, 'A', 6.0, 0.5)

        if args.evolution_mode == "single_chain": #if there are two or more chains, then calculate the number of interacting contacts
            num_inter_conts, iplddt = 1, 1
        else:
            # Eq. 5 / Methods: interface Cbeta contacts within 6 A, pLDDT floor
            # 0.5 (the paper's >50 on ESM3's 0-1 scale) -- matched to the
            # intra-chain call above rather than the old 8.0 / 0.4.
            num_inter_conts, iplddt = cbiplddt(pdb_txt, 'A', 'B', 6.0, 0.5)

        ss, max_helix, max_beta = pypsique(pdb_txt, 'A')
        #Rg, aspher = get_aspher(pdb_txt)
        #dG = dGscore(pdbtxt2bbcoord(pdb_txt), seq)
        prot_len_penalty =  1 - sigmoid(seq_len, args.prot_len_penalty, 0.12)
        max_alpha_penalty = 1 - sigmoid(max_helix, args.helix_len_penalty, 0.5)
        max_beta_penalty = 1 - sigmoid(max_beta, args.beta_len_penalty, 0.5)

        # MACREL AMP probability. HemoPI2 hemolysis is ALWAYS computed and logged as an
        # attribute (the hemo_prob column). By DEFAULT it does NOT drive selection
        # (hemo_factor = 1); pass --hemo-in-score to penalise the score via (1 - hemo_prob).
        # amp_src records WHICH scorer answered: 'macrel' or 'proxy'. MACREL is
        # defined for 10-100 residues and the biophysical surrogate is
        # substituted silently outside it, so on an arm with no length term the
        # amp_prob column can change identity mid-run. Logging the source keeps
        # that in the data rather than only on stderr.
        amp_prob, hemo_prob, amp_src = macrel_scores.get(
            seq, (calculate_samp(seq), calculate_hemo_proxy(seq), 'proxy'))
        hemo_factor = (1 - hemo_prob) if args.hemo_in_score else 1.0

        if args.evolution_mode == "single_chain":
            # iplddt and inter-chain contact term are always 1.0 for single chains — omit them.
            score = np.prod([mean_plddt,
                             ptm,
                             prot_len_penalty,
                             max_beta_penalty,
                             max_alpha_penalty,
                             amp_prob,
                             hemo_factor,
                             (num_conts + seq_len) / seq_len])
        else:
            score = np.prod([mean_plddt,
                             ptm,
                             iplddt,
                             prot_len_penalty,
                             max_beta_penalty,
                             max_alpha_penalty,
                             amp_prob,
                             hemo_factor,
                             (num_conts + seq_len) / seq_len,
                             (num_inter_conts + seq_len) / (seq_len + 1)])
        #================================SCORING================================#
        #=======================================================================#

        iterlog = pd.DataFrame({'gndx': gen_i,
                                'id': id,
                                'seq_len': seq_len,
                                'prot_len_penalty': round(prot_len_penalty, 2),
                                'max_alpha_penalty': round(max_alpha_penalty, 2),
                                'max_beta_penalty': round(max_beta_penalty, 2),
                                'ptm': round(ptm, 2),
                                'mean_plddt': round(mean_plddt, 2),
                                'num_conts': num_conts,
                                'iplddt': iplddt,
                                'num_inter_conts': num_inter_conts,
                                'sel_mode': args.selection_mode,
                                #'dG': round(dG, 3),
                                #'ptm_full': ptm_full,
                                #'cd' contact_density
                                'amp_prob': round(amp_prob, 3),
                                'amp_src': amp_src,
                                'hemo_prob': round(hemo_prob, 3),
                                'score': round(score, 3),
                                'sequence': seq,
                                'mutation': mutation,
                                'prev_id': prev_id,
                                'ss': ss}, index=[0])
        
        if new_gen.empty:
            new_gen = iterlog
        else:
            new_gen = pd.concat([new_gen, iterlog], axis=0, ignore_index=True) 
        os.system(f"gzip '{pdb_path}{id}.pdb' &")


def multimer_evolver(model, args):  
    print("evolution of interacting dimers")

global new_gen #this will be modified in the extract_results() 

#============================================================================# 
#================================FOLD_EVOLVER================================# 
 

def fold_evolver(args, model, evolver, logheader, init_gen, device) -> None:

    os.makedirs(pdb_path, exist_ok=True)
    with open(os.path.join(args.outpath, args.log), 'w') as f:
        f.write(logheader)

    condition = True
    
    #creare an initial pool of sequences with pop_size
    columns=['gndx',
             'id',
             'seq_len',
             'prot_len_penalty',
             'max_alpha_penalty',
             'max_beta_penalty',
             'ptm',
             'mean_plddt',
             'num_conts',
             'iplddt',
             'num_inter_conts',
             'sel_mode',
             #'dG',
             'amp_prob',
             'amp_src',
             'hemo_prob',
             'score',
             'sequence',
             'mutation',
             'prev_id',
             'ss']


    ancestral_memory = pd.DataFrame(columns=columns)
    ancestral_memory.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=True, sep='\t') #write header of the progress log

    #mutate seqs from init_gen and select the best N seqs for the next generation
    best_so_far = 0.0
    run_start = time.time()
    for gen_i in range(args.num_generations):
        n = 0
        global new_gen #this will be modified in the extract_results()
        new_gen = pd.DataFrame(columns=columns)
        gen_start = time.time()
        generated_sequences = []

        for prev_id, sequence in zip(init_gen.id, init_gen.sequence):
            seq, mutation_data= evolver.mutate(sequence)

            #check if the mutated sequence was already predicted
            seqmask = ancestral_memory.sequence == seq

            #if --norepeat and seq is in the ancestral_memory mutate it again
            if args.norepeat and seqmask.any():
                while seqmask.any():
                    seq, mutation_data = evolver.mutate(seq)
                    seqmask = ancestral_memory.sequence == seq

            id = "g{0}seq{1}_{2}_{3}".format(gen_i, n, prev_id, mutation_data); n+=1

            if seqmask.any(): #if sequence already exits do not predict a structure again
                repeat = ancestral_memory[seqmask].drop_duplicates(subset=['sequence'], keep='last')
                if new_gen.empty:
                    new_gen = repeat
                else:
                    new_gen = pd.concat([new_gen, repeat])
            else:
                generated_sequences.append((id, seq))

        batched_sequences = list(create_batched_sequence_dataset(generated_sequences, args.max_tokens_per_batch))
        if not batched_sequences:
            continue

        print(f"  Gen {gen_i + 1:>4}/{args.num_generations}  │  folding {len(generated_sequences)} sequences...", end='  ', flush=True)

        use_threads = device.type != 'cpu'
        trd = None
        pdbs, ptms, mean_plddts = [], [], []

        for headers, sequences in batched_sequences:
            # 1. Start scoring thread for the previous batch (GPU/MPS only)
            if use_threads and trd is not None:
                trd.start()

            # 2. Fold current batch
            if torch.backends.mps.is_available():
                torch.mps.synchronize()
            with torch.no_grad():
                esm_proteins = [ESMProtein(sequence=s) for s in sequences]
                configs = [GenerationConfig(track="structure", num_steps=args.num_recycles) for _ in sequences]
                pdbs, ptms, mean_plddts = esm2data(model.batch_generate(esm_proteins, configs))

            # 3. Wait for previous scoring thread (GPU/MPS only)
            if use_threads and trd is not None:
                trd.join()

            # 4. Run MACREL on current batch
            macrel_scores = macrel_score_batch_src(sequences)

            if use_threads:
                # 5a. GPU/MPS: overlap next fold with this scoring
                trd = threading.Thread(target=extract_results, args=(gen_i, headers, sequences, pdbs, ptms, mean_plddts, macrel_scores))
            else:
                # 5b. CPU: run scoring immediately to avoid competing for cores
                extract_results(gen_i, headers, sequences, pdbs, ptms, mean_plddts, macrel_scores)

        # Flush the last thread (GPU/MPS only)
        if use_threads and trd is not None:
            trd.start()
            trd.join()

        this_best = _print_gen_summary(gen_i, args.num_generations, new_gen, time.time() - gen_start, best_so_far)
        best_so_far = max(best_so_far, this_best)

        if ancestral_memory.empty:
            ancestral_memory = init_gen
        else:
            ancestral_memory = pd.concat([ancestral_memory, init_gen])

        #select the next generation
        init_gen = evolver.select(new_gen, init_gen, args.pop_size, args.selection_mode, args.norepeat, args.beta)
        init_gen.gndx = f'gndx{gen_i}' #assign a new gen index
        init_gen.dropna(subset=['mean_plddt']).to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=False, sep='\t')

        #Change the selection with a condition (plddt, ptm)
        if args.strong_selection_by_condition:
            if (init_gen['mean_plddt'] > 0.6).any() & (init_gen['ptm'] > 0.5).any() & condition:
                args.selection_mode = 'strong'
                condition = False
                print(f"  → Selection switched to STRONG (pLDDT > 0.6, pTM > 0.5 reached at gen {gen_i + 1})")
                with open(os.path.join(args.outpath, args.log), mode='a') as f:
                    f.write("#changing the selection mode to strong")

        #Change the selection mode after n generations
        if args.strong_selection_after_n_gen > 0:
            if (gen_i > args.strong_selection_after_n_gen) & condition:
                args.selection_mode = 'strong'
                evolver = Evolver('flatoptim')
                condition = False
                print(f"  → Selection switched to STRONG after {gen_i + 1} generations")
                with open(os.path.join(args.outpath, args.log), mode='a') as f:
                    f.write("#changing the selection mode to strong")

        #stop simulation by a condition
        if args.stop_by_condition:
            if (init_gen['mean_plddt'] > 0.85).any() & (init_gen['ptm'] > 0.75).any():
                print(f"\n  ✓ Stopping condition reached at gen {gen_i + 1}  (pLDDT > 0.85, pTM > 0.75)\n")
                break

    _print_run_summary(os.path.join(args.outpath, args.log), time.time() - run_start)


#================================FOLD_EVOLVER================================#
#============================================================================# 





#==================================================================================#
#================================INTER_FOLD_EVOLVER================================# 

def inter_fold_evolver(args, model, evolver, logheader, init_gen, device) -> None:
    if not args.initial_seq2:
        print("  Error: --initial_seq2 / -iseq2 is required for inter_chain mode.")
        sys.exit(1)

    #evolution of an interacting chain
    NZ_CP011286=":LNIIKLFHGHKYCLIFYVLP" #intergenic region from Yersinia
    PDB_1RFP=":QCRRLCYKQRCVTYCRGR" # 1RFP contains S-S bond
    PDB_6SVE=":WEKRMSRNSGRVYYFNHITNASQF" #WW domain
    PDB_5YIW=":GAMDMSWTDERVSTLKKLWLDGLSASQIAKQLGGVTRNAVIGKVHRLGL" #HTH
    PDB_4REX=":DVPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPRKAMLSQ" #4REX (170 to 207) 
    PDB_6M6W=":MNDIIINKIATIKRCIKRIQQVYGDGSQFKQDFTLQDSVILNLQRCCEACIDIANHINRQQQLGIPQSSRDSFTLLAQNNLITQPLSDNLKKMVGLRNIAVHDAQELNLDIVVHVVQHHLEDFEQFIDVIKAE" #HEPN toxin
    PDB_4OO8=":GQKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDHIVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLITQRKFDNLTKAERGGL" #CAS9 HNH
    PDB_5VGB=":GAASEIEKRQEENRKDREKAAAKFREYFPNFVGEPKSKDILKLRLYEQQHGKCLYSGKEINLGRLNEKGYVEIDHALPFSRTWDDSFNNKVLVLGSENQNKGNQTPYEYFNGKDNSREWQEFKARVETSRFPRSKKQRILLQ" #CAS9 HNH
    PDB_5O56=":SKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDHIVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLITQRKFDNLTKAERG"
    seq2 = ':' + args.initial_seq2

    os.makedirs(pdb_path, exist_ok=True)
    with open(os.path.join(args.outpath, args.log), 'w') as f:
        f.write(logheader)


    #creare an initial pool of sequences with pop_size
    columns = ['gndx',
               'id',
               'seq_len',
               'prot_len_penalty',
               'max_alpha_penalty',
               'max_beta_penalty',
               'ptm',
               'mean_plddt',
               'num_conts',
               'iplddt',
               'num_inter_conts',
               'sel_mode',
               #'dG',
               'amp_prob',
               'amp_src',
               'hemo_prob',
               'score',
               'sequence',
               'mutation',
               'prev_id',
               'ss']
      
    ancestral_memory = pd.DataFrame(columns=columns)
    ancestral_memory.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=True, sep='\t') #write header of the progress log

    #mutate seqs from init_gen and select the best n seqs for the next generation
    best_so_far = 0.0
    run_start = time.time()
    for gen_i in range(args.num_generations):
        n = 0
        global new_gen #this will be modified in the extract_results()
        new_gen = pd.DataFrame(columns=columns)
        gen_start = time.time()
        generated_sequences = []

        for prev_id, sequence in zip(init_gen.id, init_gen.sequence):
            seq, mutation_data= evolver.mutate(sequence)

            #check if the mutated sequence was already predicted
            seqmask = ancestral_memory.sequence == seq

            #if --norepeat and seq is in the ancestral_memory mutate it again
            if args.norepeat and seqmask.any():
                while seqmask.any():
                    seq, mutation_data = evolver.mutate(seq)
                    seqmask = ancestral_memory.sequence == seq

            id = "g{0}seq{1}_{2}_{3}".format(gen_i, n, prev_id, mutation_data); n+=1

            if seqmask.any(): #if sequence already exits do not predict a structure again
                repeat = ancestral_memory[seqmask].drop_duplicates(subset=['sequence'], keep='last')
                if new_gen.empty:
                    new_gen = repeat
                else:
                    new_gen = pd.concat([new_gen, repeat])
            else:
                generated_sequences.append((id, seq + seq2))

        batched_sequences = list(create_batched_sequence_dataset(generated_sequences, args.max_tokens_per_batch))
        if not batched_sequences:
            continue

        print(f"  Gen {gen_i + 1:>4}/{args.num_generations}  │  folding {len(generated_sequences)} sequences...", end='  ', flush=True)

        use_threads = device.type != 'cpu'
        trd = None
        pdbs, ptms, mean_plddts = [], [], []

        for headers, sequences in batched_sequences:
            # 1. Start scoring thread for the previous batch (GPU/MPS only)
            if use_threads and trd is not None:
                trd.start()

            # 2. Fold current batch
            if torch.backends.mps.is_available():
                torch.mps.synchronize()
            with torch.no_grad():
                linker = "GP" + "G"*30 + "PG"
                esm_proteins = [ESMProtein(sequence=s.replace(':', linker)) for s in sequences]
                configs = [GenerationConfig(track="structure", num_steps=args.num_recycles) for _ in sequences]
                pdbs, ptms, mean_plddts = esm2data(model.batch_generate(esm_proteins, configs))

            # 3. Wait for previous scoring thread (GPU/MPS only)
            if use_threads and trd is not None:
                trd.join()

            # 4. Run MACREL on chain A only
            chain_a_seqs = [s.split(':')[0] for s in sequences]
            macrel_scores = macrel_score_batch_src(chain_a_seqs)

            if use_threads:
                # 5a. GPU/MPS: overlap next fold with this scoring
                trd = threading.Thread(target=extract_results, args=(gen_i, headers, sequences, pdbs, ptms, mean_plddts, macrel_scores))
            else:
                # 5b. CPU: run scoring immediately to avoid competing for cores
                extract_results(gen_i, headers, sequences, pdbs, ptms, mean_plddts, macrel_scores)

        # Flush the last thread (GPU/MPS only)
        if use_threads and trd is not None:
            trd.start()
            trd.join()

        this_best = _print_gen_summary(gen_i, args.num_generations, new_gen, time.time() - gen_start, best_so_far)
        best_so_far = max(best_so_far, this_best)

        if ancestral_memory.empty:
            ancestral_memory = init_gen
        else:
            ancestral_memory = pd.concat([ancestral_memory, init_gen])

        #select the next generation
        init_gen = evolver.select(new_gen, init_gen, args.pop_size, args.selection_mode, args.norepeat, args.beta)
        init_gen.gndx = f'gndx{gen_i}' #assign a new gen index
        init_gen.dropna(subset=['mean_plddt']).to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=False, sep='\t')

    _print_run_summary(os.path.join(args.outpath, args.log), time.time() - run_start)


#================================INTER_FOLD_EVOLVER================================#
#==================================================================================#


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Sample sequences based on a given structure.'
    )
    parser.add_argument(
            '-em', '--evolution_mode', type=str,
            help='evolution mode: single_chain, inter_chain, multimer',
            default='single_chain',
    )
    parser.add_argument(
            '-sm', '--selection_mode', type=str,
            help='selection mode\n options: strong, weak ',
            default="weak"
    )
    parser.add_argument(
            '-b', '--beta', type=float,
            help='selection strength',
            default=1,
    )
    parser.add_argument(
            '-iseq', '--initial_seq', type=str,
            help='sequence to initiate with. "random"/"randoms" = random seqs (length via --random_seq_len); '
                 'a literal AA sequence = start from it; '
                 '"db[:spec]" = seed from the known-AMP database (spec: random|diverse|low|high|<name>, e.g. db:magainin_2; "dbaasp" alias accepted)',
            default='random'
    )
    parser.add_argument(
            '-iseq2', '--initial_seq2', type=str,
            help='second sequence'
    )
    parser.add_argument(
            '-l', '--log', type=str,
            help='log output',
            default='progress.log',
    )   
    parser.add_argument(
            '-o' ,'--outpath', type=str,
            help='output filepath for saving sampled sequences',
            default='output',
    )
    parser.add_argument(
            '-ng', '--num_generations', type=int,
            help='number of generations',
            default=100,
    )
    parser.add_argument(
            '-ps', '--pop_size', type=int,
            help='population size',
            default=10,
    )
    parser.add_argument(
            '-ed', '--evoldict', type=str,
            help='mutation rate dictionary: flatrates, codonrates, flatoptim, uniprotrates',
            default='flatrates',
    )
    parser.add_argument(
            '-pl0', '--prot_len_penalty', type=int,
            help='length penalty threshold (default 30 for AMPs)',
            default=30,
    )
    parser.add_argument(
            '-hl0', '--helix_len_penalty', type=int,
            help='helix length threshold above which the alpha-helix penalty kicks in',
            default=20,
    )
    parser.add_argument(
            '-bl0', '--beta_len_penalty', type=int,
            help='beta strand length threshold above which the beta-sheet penalty kicks in',
            default=12,
    )
    parser.add_argument(
            '--random_seq_len', type=int,
            help='a sequence to initiate with',
            default=24,
    )
    parser.add_argument(                      
            '--norepeat', action='store_true', 
            help='do not generate and/or select the same sequences more than once', 
    )
    parser.add_argument(
            '--nobackup', action='store_true', 
            help='overwrite files if exists',
    )
    parser.add_argument(
            '--stop_by_condition', action='store_true', 
            help='experimental',
    )
    parser.add_argument(
            '--strong_selection_by_condition', action='store_true', 
            help='experimental',
    )
    parser.add_argument(
            '--strong_selection_after_n_gen', type=int,
            help='switch to strong (top-N) selection after this many generations (default 4500 = effectively disabled)',
            default=4500,
    )
    # parser.add_argument(
    #         '--continue', action='store_true', 
    #         help='',
    # )
    parser.add_argument(
            '--num-recycles',
            type=int,
            default=1,
            help="Number of ESM3 denoising steps per structure prediction (1=fast, 4-8=higher quality). Default 1.",
    )
    parser.add_argument(
            '--max-tokens-per-batch',
            type=int,
            default=512, # 2048+ works fine with A100/V100; 512 is safe for CPU
            help="Maximum number of tokens per forward-pass. Lower this if you run out of memory. "
            "Default 512 is conservative and safe for CPU runs."
    )
    parser.add_argument(
            '--start', choices=['random', 'randoms', 'existing', 'mix', 'file', 'seq'],
            default=None,
            help='HIGH-LEVEL choice of the starting population (overrides -iseq): '
                 'random = one random sequence copied across the population; '
                 'randoms = a different random sequence per slot; '
                 'existing = known AMPs from the database (diverse); '
                 'mix = existing AMPs + random sequences at the existing mean length (see --mix-frac); '
                 'file = a fixed population from --start-file; '
                 'seq = a literal sequence given via -iseq.',
    )
    parser.add_argument(
            '--mix-frac', type=float, default=0.5,
            help='fraction of the population drawn from existing AMPs when --start mix (default 0.5)',
    )
    parser.add_argument(
            '--start-file', type=str, default=None,
            help='FASTA of a fixed initial population (one record per member) when --start file',
    )
    parser.add_argument(
            '--hemo-in-score', action='store_true',
            help='include HemoPI2 hemolysis in the fitness score via (1 - hemo_prob). '
                 'DEFAULT: hemolysis is computed & logged as an attribute (hemo_prob column) '
                 'but does NOT drive selection.',
    )
    parser.add_argument(
            '--threads', type=int, default=None,
            help='cap CPU threads for ESM3 folding (default: all cores). Set this when '
                 'running several PFES processes on one node, e.g. --threads 6 for 8 '
                 'parallel runs on a 48-core machine.',
    )

    args = parser.parse_args()

    # --start is a friendly front-end for -iseq: translate it into the initial_seq
    # string the seeding code understands. -iseq still works on its own.
    if args.start:
        if args.start in ('random', 'randoms'):
            args.initial_seq = args.start
        elif args.start == 'existing':
            args.initial_seq = 'db:diverse'
        elif args.start == 'mix':
            args.initial_seq = f'db:mix:{args.mix_frac}'
        elif args.start == 'file':
            if not args.start_file:
                print("  Error: --start file requires --start-file <path.faa>")
                sys.exit(1)
            args.initial_seq = f'file:{args.start_file}'
        elif args.start == 'seq':
            pass  # use whatever -iseq provided

    evolver = Evolver(args.evoldict)

    now = datetime.now() # current date and time
    date_now = now.strftime("%d-%b-%Y")
    time_now = now.strftime("%H:%M:%S")
    


    logheader = f'''#======================== PFESv0.1 ========================#
#====================== {date_now} =======================#
#======================== {time_now} ========================#
#WD: {os.getcwd()}
#$pfes.py {' '.join(sys.argv[1:])}
#
#====================  pfes input params ==================#
#
#--evolution_mode, -em \t\t = {args.evolution_mode}
#--selection_mode, -sm\t\t = {args.selection_mode}
#--initial_seq, -iseq\t\t = {args.initial_seq}
#--pop_size, -ps\t\t = {args.pop_size}
#--evoldict, -ed\t\t = {args.evoldict}
#--log, -l\t\t\t = {args.log}
#--outpath, -o\t\t\t = {args.outpath}
#--random_seq_len\t\t = {args.random_seq_len}
#--beta, -b\t\t\t = {args.beta}
#--helix_len_penalty, -hl0\t = {args.helix_len_penalty}
#--prot_len_penalty, -pl0\t = {args.prot_len_penalty}
#--num_generations, -ng\t\t = {args.num_generations}
#--strong_selection_after_n_gen\t\t = {args.strong_selection_after_n_gen}
#--norepeat\t\t\t = {args.norepeat}
#--nobackup\t\t\t = {args.nobackup}
#--num-recycles\t\t\t = {args.num_recycles}
#--max-tokens-per-batch\t\t = {args.max_tokens_per_batch}
# evolution dictionary = {evolver.evoldict}
# evolution dictionary normalized = {evolver.evoldict_normal}
#==========================================================#
'''
    
    _print_startup(args, evolver, date_now, time_now)

    #backup if output directory exists
    if args.nobackup:
        if os.path.isdir(args.outpath):
            print(f'\nWARNING! Directory {args.outpath} exists, it will be replaced!')
            shutil.rmtree(args.outpath)
        os.makedirs(args.outpath)
    else:
        backup_output(args.outpath)

    pdb_path = args.outpath + '/structures/' 

    #create the initial generation
    if args.initial_seq.startswith('file:'):
        # Load a FIXED initial population from a FASTA (one record per member).
        # Build it once with `python amp_db.py --make-init ... -o init_pop.faa`
        # and pass the SAME file to every branch for an identical start.
        import amp_db
        init_path = args.initial_seq.split(':', 1)[1]
        if not os.path.isfile(init_path):
            print(f"  Error: init population file not found: {init_path}")
            print(f"  Build it with:  python amp_db.py --make-init --pop {args.pop_size} -o {init_path}")
            sys.exit(1)
        recs = amp_db.parse_fasta(init_path)
        if not recs:
            print(f"  Error: no sequences in {init_path}")
            sys.exit(1)
        seqs = [r['seq'] for r in recs]
        names = [r['name'] for r in recs]
        if len(seqs) < args.pop_size:        # cycle if the file is smaller than pop
            reps = args.pop_size // len(seqs) + 1
            seqs = (seqs * reps)[:args.pop_size]
            names = (names * reps)[:args.pop_size]
        else:
            seqs, names = seqs[:args.pop_size], names[:args.pop_size]
        print(f"  loaded fixed init population: {len(seqs)} members from "
              f"{args.initial_seq.split(':', 1)[1]}")
        init_gen = pd.DataFrame({'id': [f'init_{n}' for n in names],
                                 'sequence': seqs,
                                 'score': [0.001] * len(seqs)})
    elif args.initial_seq == 'random':
        randomsequence = evolver.randomseq(args.random_seq_len)
        init_gen = pd.DataFrame({'id': ['init_seq'] * args.pop_size,
                                 'sequence': [randomsequence] * args.pop_size,
                                 'score': [0.001] * args.pop_size})
    elif args.initial_seq == 'randoms':
        init_gen = pd.DataFrame({'id': [f'init_seq{i}' for i in range(args.pop_size)], 
                                 'sequence': [evolver.randomseq(args.random_seq_len) for i in range(args.pop_size)],
                                 'score': [0.001] * args.pop_size})
    #elif args.initial_seq == 'c':
    #    init_gen = pd.read_csv('test.chk', sep='\t')
    elif args.initial_seq in ('db', 'dbaasp') or args.initial_seq.startswith(('db:', 'dbaasp:')):
        import amp_db
        spec = args.initial_seq.split(':', 1)[1] if ':' in args.initial_seq else 'random'
        if spec.startswith('mix'):
            # Mixed start: existing AMPs + random sequences whose length equals the
            # existing seeds' MEAN length. spec 'mix' = 50/50; 'mix:<frac>' sets the
            # existing fraction, e.g. -iseq db:mix:0.5 (50 existing / 50 random at -ps 100).
            try:
                frac = float(spec.split(':')[1]) if ':' in spec else 0.5
            except (IndexError, ValueError):
                frac = 0.5
            seeds = amp_db.mixed_population(args.pop_size, frac,
                                            evolver.randomseq, args.random_seq_len)
            n_rand = sum(1 for n, _ in seeds if n.startswith('rand'))
            rand_len = next((len(s) for n, s in seeds if n.startswith('rand')), args.random_seq_len)
            print(f"  seeded MIX: {len(seeds) - n_rand} existing AMPs + {n_rand} random "
                  f"@ {rand_len} aa (= mean length of the existing seeds)")
        else:
            seeds = amp_db.seeds_for_population(spec, args.pop_size)
            uniq = sorted({n for n, _ in seeds})
            print(f"  seeded from AMP DB ({spec}): {', '.join(uniq[:5])}"
                  f"{' …' if len(uniq) > 5 else ''}")
        init_gen = pd.DataFrame({'id': [f'init_{n}' for n, _ in seeds],
                                 'sequence': [s for _, s in seeds],
                                 'score': [0.001] * len(seeds)})
    else:
        init_gen = pd.DataFrame({'id': ['init_seq'] * args.pop_size,
                                 'sequence': [args.initial_seq] * args.pop_size,
                                 'score': [0.001] * args.pop_size})
    

    # HuggingFace token for gated ESM3 model access. Use the HF_TOKEN env var if set,
    # otherwise fall back to a cached `huggingface-cli login` token so you only enter
    # it once and never need to re-export it.
    if "HF_TOKEN" not in os.environ:
        _tok = None
        try:
            from huggingface_hub import get_token as _hf_get_token
            _tok = _hf_get_token()
        except Exception:
            try:
                from huggingface_hub import HfFolder
                _tok = HfFolder.get_token()
            except Exception:
                _tok = None
        if _tok:
            os.environ["HF_TOKEN"] = _tok
        else:
            print(f"  Error: no HuggingFace token found.")
            print(f"  Do ONE of (once):")
            print(f"    huggingface-cli login          # stores it permanently in ~/.cache/huggingface")
            print(f"    export HF_TOKEN=your_token      # add this line to ~/.bashrc to persist")
            sys.exit(1)

    print(f"  Loading ESM3 (esm3-sm-open-v1)...", end='  ', flush=True)
    try:
        model = ESM3.from_pretrained("esm3-sm-open-v1")
    except Exception as e:
        print(f"\n  Error loading ESM3: {e}")
        print(f"  Ensure you have accepted the model license at huggingface.co and HF_TOKEN is valid.")
        sys.exit(1)

    if torch.backends.mps.is_available():
        device = torch.device("mps")
        os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.0"
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
        # --threads caps CPU threads so several runs can share one node without thrashing.
        n_threads = args.threads if args.threads else (os.cpu_count() or 1)
        torch.set_num_threads(n_threads)
        print(f"\n  Note: running on CPU ({n_threads} threads) — expect slow fold times.")
        print(f"        Recommended: -ps 4 -ng 20 --max-tokens-per-batch 256")
    model = model.eval().to(device)
    print(f"ready  [{device}]\n")

    if args.evolution_mode == "single_chain":
        fold_evolver(args, model, evolver, logheader, init_gen, device)
    elif args.evolution_mode == "inter_chain":
        inter_fold_evolver(args, model, evolver, logheader, init_gen, device)
    elif args.evolution_mode == "multimer":
        print("  multimer mode is not yet implemented")
    else:
        print(f"  Unknown evolution mode: '{args.evolution_mode}'  (options: single_chain, inter_chain, multimer)")






