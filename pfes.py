import argparse
import os
import sys
import shutil
import pandas as pd
import numpy as np
import typing as T
import threading
import gzip
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
from score import get_nconts, cbiplddt, calculate_samp, macrel_score_batch
from psique import pypsique

from esm.sdk.api import ESMProtein, GenerationConfig
from esm.models.esm3 import ESM3


def backup_output(outpath):
    print(f'\nSaving output files to {args.outpath}')
    if os.path.isdir(outpath): 
        backup_list = []
        last_backup = int()
        for dir_name in os.listdir():
            if dir_name.startswith(outpath + '.'):
                backup=(dir_name.split('.')[-1])
                if backup.isdigit(): 
                    backup_list.append(backup)
                    last_backup = int(max(backup_list))
        print(f'\n{outpath} already exists, renaming it to {outpath}.{str(last_backup +  1)}')
        os.replace(outpath, outpath + '.' + str(last_backup +  1))


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
    print(f"  output:     {args.outpath}/{args.log}")
    print(f"{'═'*_W}\n")

def _print_gen_summary(gen_i, num_gen, new_gen, elapsed):
    top = new_gen.sort_values('score', ascending=False)
    has_samp = 's_amp' in top.columns
    has_hemo = 'hemo_prob' in top.columns
    print(f"\n  ── Gen {gen_i + 1}/{num_gen}  ({elapsed:.1f}s) {'─' * (_W - 22)}")
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
    for _, row in top.iterrows():
        seq = str(row.sequence)
        seq_disp = seq[:20] + '…' if len(seq) > 20 else seq
        mut = str(row.mutation)[:14]
        line = f"  {float(row.score):>6.3f}  {float(row.mean_plddt):>6.3f}  {float(row.ptm):>5.3f}"
        line += f"  {float(row.s_amp):>5.3f}" if has_samp else ""
        line += f"  {float(row.hemo_prob):>5.3f}" if has_hemo else ""
        line += f"  {int(row.seq_len):>4}  {mut:<14}  {seq_disp}"
        print(line)
    print()


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
        num_conts, _mean_plddt_ = get_nconts(pdb_txt, 'A', 8.0, 0.5) #plddt is better only for chain A and for residues > 50

        if args.evolution_mode == "single_chain": #if there are two or more chains, then calculate the number of interacting contacts
            num_inter_conts, iplddt = 1, 1
        else:
            num_inter_conts, iplddt = cbiplddt(pdb_txt, 'A', 'B', 8.0, 0.4)

        ss, max_helix, max_beta = pypsique(pdb_txt, 'A')
        #Rg, aspher = get_aspher(pdb_txt)
        #dG = dGscore(pdbtxt2bbcoord(pdb_txt), seq)
        prot_len_penalty =  1 - sigmoid(seq_len, args.prot_len_penalty, 0.2)
        max_alpha_penalty = 1 - sigmoid(max_helix, args.helix_len_penalty, 0.5)
        max_beta_penalty = 1 - sigmoid(max_beta, args.beta_len_penalty, 0.6)

        # MACREL AMP probability and hemolytic penalty
        s_amp, hemo_prob = macrel_scores.get(seq, (calculate_samp(seq), 0.0))

        score  = np.prod([mean_plddt,           #[0, 1]
                          ptm,                  #[0, 1]
                          iplddt,               #[0, 1]
                          prot_len_penalty,     #[0, 1]
                          max_beta_penalty,     #[0, 1]
                          max_alpha_penalty,    #[0, 1]
                          s_amp,                # MACREL AMP probability
                          1 - hemo_prob,        # hemolytic penalty
                          #dG, #~[0, inf]
                          (num_conts + seq_len) / seq_len,
                          (num_inter_conts + seq_len) / (seq_len + 1) # change this to sigmod so the number of inter contacts > X would not increase the score
                          ])
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
                                's_amp': round(s_amp, 3),
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
             's_amp',
             'hemo_prob',
             'score',
             'sequence',
             'mutation',
             'prev_id',
             'ss']


    ancestral_memory = pd.DataFrame(columns=columns)
    ancestral_memory.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=True, sep='\t') #write header of the progress log

    #mutate seqs from init_gen and select the best N seqs for the next generation
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
                configs = [GenerationConfig(track="structure", num_steps=1) for _ in sequences]
                pdbs, ptms, mean_plddts = esm2data(model.batch_generate(esm_proteins, configs))

            # 3. Wait for previous scoring thread (GPU/MPS only)
            if use_threads and trd is not None:
                trd.join()

            # 4. Run MACREL on current batch
            macrel_scores = macrel_score_batch(sequences)

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

        _print_gen_summary(gen_i, args.num_generations, new_gen, time.time() - gen_start)

        if ancestral_memory.empty:
            ancestral_memory = init_gen
        else:
            ancestral_memory = pd.concat([ancestral_memory, init_gen])

        #select the next generation
        init_gen = evolver.select(new_gen, init_gen, args.pop_size, args.selection_mode, args.norepeat, args.beta)
        init_gen.gndx = f'gndx{gen_i}' #assign a new gen index
        init_gen.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=False, sep='\t')

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

 
#================================FOLD_EVOLVER================================# 
#============================================================================# 





#==================================================================================#
#================================INTER_FOLD_EVOLVER================================# 

def inter_fold_evolver(args, model, evolver, logheader, init_gen, device) -> None:

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
               's_amp',
               'hemo_prob',
               'score',
               'sequence',
               'mutation',
               'prev_id',
               'ss']
      
    ancestral_memory = pd.DataFrame(columns=columns)
    ancestral_memory.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=True, sep='\t') #write header of the progress log
    
    #mutate seqs from init_gen and select the best n seqs for the next generation
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
                configs = [GenerationConfig(track="structure", num_steps=1) for _ in sequences]
                pdbs, ptms, mean_plddts = esm2data(model.batch_generate(esm_proteins, configs))

            # 3. Wait for previous scoring thread (GPU/MPS only)
            if use_threads and trd is not None:
                trd.join()

            # 4. Run MACREL on chain A only
            chain_a_seqs = [s.split(':')[0] for s in sequences]
            macrel_scores = macrel_score_batch(chain_a_seqs)

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

        _print_gen_summary(gen_i, args.num_generations, new_gen, time.time() - gen_start)

        if ancestral_memory.empty:
            ancestral_memory = init_gen
        else:
            ancestral_memory = pd.concat([ancestral_memory, init_gen])

        #select the next generation
        init_gen = evolver.select(new_gen, init_gen, args.pop_size, args.selection_mode, args.norepeat, args.beta)
        init_gen.gndx = f'gndx{gen_i}' #assign a new gen index
        init_gen.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=False, sep='\t')

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
            help='a sequence to initiate with, if "random" pop_size random sequences will be generated, the length of the random sequences can be assigned with "--random_seq_len"',
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
            help='population size',
            default='flatrates',
    )
    parser.add_argument(
            '-pl0', '--prot_len_penalty', type=int,
            help='length penalty threshold (default 30 for AMPs)',
            default=30,
    )
    parser.add_argument(
            '-hl0', '--helix_len_penalty', type=int,
            help='population size',
            default=20,
    )
    parser.add_argument(
            '-bl0', '--beta_len_penalty', type=int,
            help='population size',
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
            help='',
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
            help="Number of recycles to run. Defaults to number used in training (4).",
    )
    parser.add_argument(
            '--max-tokens-per-batch',
            type=int,
            default=512, # 2048+ works fine with A100/V100; 512 is safe for CPU
            help="Maximum number of tokens per forward-pass. Lower this if you run out of memory. "
            "Default 512 is conservative and safe for CPU runs."
    )

    args = parser.parse_args()
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
    if args.initial_seq == 'random':
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
    else: 
        init_gen = pd.DataFrame({'id': ['init_seq'] * args.pop_size, 
                                 'sequence': [args.initial_seq] * args.pop_size,
                                 'score': [0.001] * args.pop_size})
    

    # Require HuggingFace token for gated ESM3 model access
    if "HF_TOKEN" not in os.environ:
        print(f"  Error: HF_TOKEN environment variable not set.")
        print(f"  Run:  export HF_TOKEN=your_huggingface_token")
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
        n_threads = os.cpu_count() or 1
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






