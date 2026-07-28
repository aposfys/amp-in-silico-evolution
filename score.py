import numpy as np
import sys, os
import math
import pandas as pd


def get_aspher(pdb_txt):

    #from DOI: 10.1016/j.bpj.2018.01.002 (HullRad)

    def distance(ax,ay,az, bx,by,bz):
        #Euclidean distance between two atoms
        return math.sqrt((ax - bx)**2.0 + (ay - by)**2.0 + (az - bz)**2.0)

    def model_from_pdb(pdb_txt):
        
        all_atm_rec = []
        # Get all relevant atoms even if in wrong order
        for line in pdb_txt.splitlines():
            if (line[:4] == 'ATOM'):
                all_atm_rec.append(line)
        # Convert all_atm_rec to multi-item list
        all_atm_array = [['X' for j in range(8)] for i in range(len(all_atm_rec))]
        for row in range(len(all_atm_rec)):
            all_atm_array[row][0] = row					#Atom Index
            all_atm_array[row][1] = (all_atm_rec[row][11:16]).strip()	#Atom Name
            all_atm_array[row][2] = (all_atm_rec[row][17:20]).strip()	#Residue Name
            all_atm_array[row][3] = (all_atm_rec[row][20:22]).strip()	#ChainID
            all_atm_array[row][4] = (all_atm_rec[row][22:26]).strip()	#Residue Number
            all_atm_array[row][5] = (all_atm_rec[row][30:38]).strip()	#x
            all_atm_array[row][6] = (all_atm_rec[row][38:46]).strip()	#y
            all_atm_array[row][7] = (all_atm_rec[row][46:54]).strip()	#z
        return all_atm_array

    all_atm_rec = model_from_pdb(pdb_txt)

    # Radius of Gyration
    # Calc center of mass
    X = 0.0
    Y = 0.0
    Z = 0.0
    tot = 0.0
    for row in range(len(all_atm_rec)):
        X = X + (float(all_atm_rec[row][5]))
        Y = Y + (float(all_atm_rec[row][6]))
        Z = Z + (float(all_atm_rec[row][7]))
        tot += 1
    com_x = (X/tot)
    com_y = (Y/tot)
    com_z = (Z/tot)
    Rg2  = 0.0
    for row in range(len(all_atm_rec)):
        Rg2 += ((distance(com_x, com_y, com_z, float(all_atm_rec[row][5]),\
            float(all_atm_rec[row][6]), float(all_atm_rec[row][7])))**2)
    Rg = math.sqrt(Rg2/tot)

    asphr = 0.0	
    Ixx,Ixy,Ixz,Iyy,Iyz,Izz = 0,0,0,0,0,0
    for row in range(len(all_atm_rec)):
        Ixx += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][5])) - com_x)
        Ixy += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][6])) - com_y)
        Ixz += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][7])) - com_z)
        Iyy += ((float(all_atm_rec[row][6])) - com_y) * ((float(all_atm_rec[row][6])) - com_y)
        Iyz += ((float(all_atm_rec[row][6])) - com_y) * ((float(all_atm_rec[row][7])) - com_z)
        Izz += ((float(all_atm_rec[row][7])) - com_z) * ((float(all_atm_rec[row][7])) - com_z)
    Ixx= Ixx/row
    Iyy= Iyy/row
    Izz= Izz/row
    Ixy= Ixy/row
    Ixz= Ixz/row
    Iyz= Iyz/row
    gyration_tensor = [[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]

    #print(gyration_tensor)
    evals, evecs = np.linalg.eig(gyration_tensor)
    L1 = evals[0]
    L2 = evals[1]
    L3 = evals[2]
    asphr = ((L1 - L2)**2 + (L2 - L3)**2 + (L1 - L3)**2)/(2.0*((L1 + L2 + L3)**2))

    return(Rg, asphr)
 


def get_nconts(pdb_txt, chain="A", distance_cutoff=6.0, plddt_cutoff=0): 
    """
    Calculates number of contaict in a protein.

    """
    
    # Get all C-alpha atoms with specific pLDDT cutoff
    ca_data, plddt = [],[]
    for line in pdb_txt.splitlines():
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            chain_id = line[21]
            atom_name = line[12:16].strip()
            
            if chain_id == chain:
                try:
                    b_factor = float(line[60:66].strip())
                except:
                    b_factor = 0.0
                plddt.append(b_factor)
                
                if atom_name == 'CA' and b_factor > plddt_cutoff:
                    try:
                        res_seq = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        ca_data.append([res_seq, np.array([x, y, z]), b_factor])
                    except ValueError:
                        pass
    
    if len(ca_data) == 0:
        mean_plddt = np.mean(np.array(plddt)) if len(plddt) > 0 else 0.0
        return(1, round(mean_plddt * 0.01, 2))
    else:
        coords = np.array([item[1] for item in ca_data])  # Extract coordinates
        mean_plddt = np.mean(np.array(plddt))
        n_atoms = len(coords)
        #pairs_data = np.zeros((0, 5))

        distances_matrix = np.linalg.norm(coords[:, None] - coords, axis=2)
        row = 0
        for i in range(n_atoms):
            for j in range(i + 4, n_atoms): # do not calc dist between atoms i, ... i+4
                if distances_matrix[i, j] < distance_cutoff:
                    #pairs_data = np.append(pairs_data, [[row, ca_data[i][0], ca_data[j][0], np.mean([ca_data[i][2], ca_data[j][2]]), distances_matrix[i, j]]], axis=0)
                    row += 1
        
        return(row+1, round(mean_plddt * 0.01, 2))

#TODO check how fast this is?
def get_nconts_allatom(pdb_txt, chain="A", distance_cutoff=4.5, plddt_cutoff=0): 
    nconts = "nconts"
    contact_density = "nconts/seq_len"
    return(nconts, contact_density)


def get_inter_nconts(pdb_txt, chainA='A', chainB='B', distance_cutoff=6.0, plddt_cutoff=0): 
    """
    Calculates number of contaict between two protein chains
    returns a tuple (number of contacts, average plddt of residues with plddt > plddt_cutoff)
    """

    # Get all C-beta atoms with specific pLDDT cutoff
    cb_data_A, cb_data_B, = [], []
    for line in pdb_txt.splitlines():
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            chain_id = line[21]
            atom_name = line[12:16].strip()
            try:
                b_factor = float(line[60:66].strip())
            except:
                b_factor = 0.0
                
            if atom_name == 'CA' and b_factor > plddt_cutoff:
                try:
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    if chain_id == chainA:
                        cb_data_A.append([res_seq, np.array([x, y, z]), b_factor])
                    elif chain_id == chainB:
                        cb_data_B.append([res_seq, np.array([x, y, z]), b_factor])
                except ValueError:
                    pass

    if len(cb_data_A) == 0 or len(cb_data_B) == 0:
        return(1, 1)
    else:
        Acoords = np.array([item[1] for item in cb_data_A])
        Bcoords = np.array([item[1] for item in cb_data_B])
        CA_pLDDT_A = np.mean(np.array([item[2] for item in cb_data_A]))
        distances_matrix = np.linalg.norm(Acoords[:, None] - Bcoords, axis=2)
        n_contacts = (distances_matrix <= distance_cutoff).sum()
        return(n_contacts, round(CA_pLDDT_A * 0.01, 2))


def cbiplddt(pdb_txt, chainA='A', chainB='B', distance_cutoff=6.0, plddt_cutoff=0):
    """
    Calculates number of contaict between two protein chains and iPLDDT
    """

    # Get all C-beta atoms with specific pLDDT cutoff
    cbeta_atom = []
    for line in pdb_txt.splitlines():
        if (line.startswith('ATOM  ') or line.startswith('HETATM')) and line[12:16].strip() == "CA":
            cbeta_atom.append(line)
            
    cbeta_array = [['X' for j in range(8)] for i in range(len(cbeta_atom))]
    for row in range(len(cbeta_atom)):
        cbeta_array[row][0] = row					#Index
        cbeta_array[row][1] = (cbeta_atom[row][22:26]).strip()	#Residue Number
        cbeta_array[row][2] = (cbeta_atom[row][30:38]).strip()	#xyz
        cbeta_array[row][3] = (cbeta_atom[row][38:46]).strip()	#xyz
        cbeta_array[row][4] = (cbeta_atom[row][46:54]).strip()	#xyz
        cbeta_array[row][5] = (cbeta_atom[row][60:66]).strip()	#pLDDT 
        cbeta_array[row][6] = (cbeta_atom[row][21:22]).strip()	#ChainID
        cbeta_array[row][7] = (cbeta_atom[row][17:20]).strip()	#Residue Name

    cb_data_A, cb_data_B, = [], []
    for row in range(len(cbeta_array)):
        if cbeta_array[row][6] == chainA and float(cbeta_array[row][5]) > plddt_cutoff:
            cb_data_A.append(cbeta_array[row][0:6])
        if cbeta_array[row][6] == chainB and float(cbeta_array[row][5]) > plddt_cutoff:
            cb_data_B.append(cbeta_array[row][0:6])
            #cb_data_A = np.array(cb_data_A, dtype='float32')
            #cb_data_B = np.array(cb_data_B, dtype='float32')
    if len(cb_data_A) == 0 or len(cb_data_B) == 0: 
        return(1, 0.01)
    else:    
        #Acoords = cb_data_A[:,2:5]        
        #Bcoords = cb_data_B[:,2:5]
        Acoords = np.array([item[2:5] for item in cb_data_A], dtype="float32")
        Bcoords = np.array([item[2:5] for item in cb_data_B], dtype="float32")
        distances_matrix = np.linalg.norm(Acoords[:, None] - Bcoords, axis=2)
        #contact_map = distances_matrix.copy()
        #contact_map[contact_map <= distance_cutoff] = 1
        #contact_map[contact_map > distance_cutoff] = 0
        matrix_mask = distances_matrix <= distance_cutoff
        n_contacts = matrix_mask.sum()
        if n_contacts == 0:
            return(1,0.01)
        else:
            inteface_ndx = np.where(matrix_mask)
            AiPLDDT = np.array([cb_data_A[i][5] for i in np.unique(inteface_ndx[0])],dtype=float)
            BiPLDDT = np.array([cb_data_B[i][5] for i in np.unique(inteface_ndx[1])],dtype=float)
            #iPLDDT = np.concatenate([AiPLDDT, BiPLDDT]).mean()
            iPLDDT = AiPLDDT.mean()
            return(n_contacts, round(iPLDDT * 0.01, 2))


def iplddt_all_atom(pdb_txt, chainA='A', chainB='B', distance_cutoff=6.0,):
    iplddt_all_atom = 'not ready yet'
    return iplddt_all_atom



if __name__ == '__main__':
    input_pdb_path = str(sys.argv[1])

    if os.path.isfile(input_pdb_path):
        with open(input_pdb_path, 'r') as file:
            pdb_txt = file.read()

        print("inner contacts, plddt:" + str(get_nconts(pdb_txt, "A", 8.0, 0)))
        print("inter contacts, iplddt:" + str(cbiplddt(pdb_txt, "A", "B", 8.0, 0)))



# Eisenberg consensus hydrophobicity scale
eisenberg_scale = {
    'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37,
    'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
    'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62,
    'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76
}

def calculate_samp(sequence):
    seq_len = len(sequence)
    if seq_len == 0:
        return 1.0

    # 1. Net Positive Charge Reward (S_charge)
    # count(R,K) + 0.1*count(H) - count(D,E)
    count_R = sequence.count('R')
    count_K = sequence.count('K')
    count_H = sequence.count('H')
    count_D = sequence.count('D')
    count_E = sequence.count('E')
    
    charge = count_R + count_K + 0.1 * count_H - (count_D + count_E)
    s_charge = 1.0 / (1.0 + math.exp(-(charge - 2.0))) # sigmoid(charge - 2)

    # 2. Hydrophobic Ratio Reward (S_hydro)
    hydro_residues = {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'}
    count_hydro = sum(1 for res in sequence if res in hydro_residues)
    ratio = count_hydro / seq_len
    s_hydro = max(1.0 - abs(ratio - 0.50) * 4, 0.01)

    # 3. Amphipathicity Reward (S_amphi)
    # Window of N=11, angle delta=100 degrees
    N = 11
    delta = math.radians(100)
    max_mu_H = 0.0
    
    if seq_len >= N:
        for i in range(seq_len - N + 1):
            window = sequence[i:i+N]
            sum_sin = 0.0
            sum_cos = 0.0
            for j, res in enumerate(window):
                h_i = eisenberg_scale.get(res, 0.0)
                sum_sin += h_i * math.sin(j * delta)
                sum_cos += h_i * math.cos(j * delta)
            
            mu_H = math.sqrt(sum_sin**2 + sum_cos**2) / N
            if mu_H > max_mu_H:
                max_mu_H = mu_H
    else:
        # For sequences shorter than 11, compute over the entire sequence
        sum_sin = 0.0
        sum_cos = 0.0
        for j, res in enumerate(sequence):
            h_i = eisenberg_scale.get(res, 0.0)
            sum_sin += h_i * math.sin(j * delta)
            sum_cos += h_i * math.cos(j * delta)
        max_mu_H = math.sqrt(sum_sin**2 + sum_cos**2) / seq_len if seq_len > 0 else 0.0

    s_amphi = math.tanh(max_mu_H * 5)
    if s_amphi < 0:
        s_amphi = 0.01 # Fallback to prevent negative score inside the power calculation

    s_amp = (s_charge * s_hydro * s_amphi) ** (1/3)
    return s_amp


def calculate_hemo_proxy(sequence):
    """
    Biophysical proxy for hemolytic probability.

    MACREL's hemolytic classifier outputs 0.000 for all cationic amphipathic
    sequences evolved here, providing no selection pressure.

    Formula: sigmoid(hydro_ratio * 10 - min(charge_pH7.4, 8.0) * 0.5 - 2.0)

    The charge is capped at +8 to prevent extreme K/R accumulation from
    producing an artificially near-zero hemo score. Without the cap, sequences
    with charge > 15 (common in MACREL-only runs) saturate the term, allowing
    unlimited hydrophobicity growth with no hemolytic penalty.

    Calibration:
    - mastoparan (hemolytic):  hydro=0.71, charge=3  → 0.97  ✓
    - melittin   (hemolytic):  hydro=0.54, charge=5  → 0.62  (underpredicted)
    - magainin2  (safe AMP):   hydro=0.44, charge=4  → 0.70  (known failure —
        biophysical proxies cannot distinguish magainin-2 from hemolytic peptides
        on hydrophobicity + charge alone; a true ML predictor is needed)
    - best pfes seq 26aa:      hydro=0.50, charge=12 → 0.27  (reasonable)
    - best macrel seq 78aa:    hydro=0.42, charge=19 → 0.27  (capped; was 0.001)

    Returns a value in [0, 1].
    """
    seq_len = len(sequence)
    if seq_len == 0:
        return 0.0

    try:
        import peptides as _pep
        charge = _pep.Peptide(sequence).charge(pH=7.4)
    except Exception:
        charge = (sequence.count('R') + sequence.count('K')
                  + 0.1 * sequence.count('H')
                  - sequence.count('D') - sequence.count('E'))

    hydro_residues = {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'}
    hydro_ratio = sum(1 for r in sequence if r in hydro_residues) / seq_len

    charge_capped = min(charge, 8.0)
    logit = hydro_ratio * 10.0 - charge_capped * 0.5 - 2.0
    hemo = 1.0 / (1.0 + math.exp(-logit))
    return round(min(max(hemo, 0.0), 1.0), 4)


def hemopi2_score_batch(sequences, model=1):
    """
    Predict hemolytic probability with HemoPI2 (Raghava lab, 2024).
    Returns dict {sequence: hemo_probability in [0, 1]}, higher = more hemolytic.

    HemoPI2 replaces the old biophysical calculate_hemo_proxy. It is an ML model
    trained on 1,926 experimentally validated hemolytic peptides (AUROC ~0.92)
    and produces a real, varying signal on the cationic amphipathic peptides PFES
    evolves — unlike MACREL's hemolytic RF, which saturates to ~0.000 out of its
    training distribution and gave no selection gradient.

    Install:  pip install hemopi2
    CLI:      hemopi2_classification -i in.fa -o out.csv -m {1=RF,2=RF+MERCI,3=ESM2,4=ESM+MERCI}
    model=1 (Random Forest) is used here: fastest and dependency-light (no Perl/MERCI),
    appropriate for per-generation scoring inside the evolutionary loop. Switch to
    model=3 (ESM2-t6) for slightly higher accuracy at extra cost.

    Falls back per-sequence to calculate_hemo_proxy if HemoPI2 is not installed
    or fails.
    """
    import subprocess, tempfile
    if os.environ.get("PFES_SKIP_HEMO") == "1":
        # Hemolysis is attribute-only; skip the per-generation HemoPI2 subprocess
        # (avoids re-importing transformers/torch every generation). Score the
        # final candidates with HemoPI2 once, post-hoc. hemo_prob logged as 0.0.
        return {seq: 0.0 for seq in sequences}
    fallback = {seq: calculate_hemo_proxy(seq) for seq in sequences}
    if not sequences:
        return {}
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, 'hemo_in.fa')
            out_name = 'hemo_out.csv'          # HemoPI2 writes to <wd>/<-o>, so -o must be a bare name
            out_path = os.path.join(tmpdir, out_name)
            order = []
            with open(fasta_path, 'w') as fh:
                for i, seq in enumerate(sequences):
                    fh.write(f'>seq{i}\n{seq}\n')
                    order.append(seq)
            try:
                result = subprocess.run(
                    ['hemopi2_classification', '-i', fasta_path,
                     '-o', out_name, '-m', str(model), '-j', '1',
                     '-wd', tmpdir],
                    capture_output=True, text=True, timeout=600, cwd=tmpdir,
                )
            except FileNotFoundError:
                sys.stderr.write(
                    '  Warning: hemopi2 not installed — pip install hemopi2\n'
                    '           falling back to biophysical hemo proxy\n'
                )
                return fallback
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip() or result.stdout.strip())
            if not os.path.isfile(out_path):
                raise FileNotFoundError(f'no hemopi2 output at {out_path}')

            df = pd.read_csv(out_path)
            # HemoPI2 output column names are not contractually fixed across
            # versions, so locate the hemolytic score robustly: prefer a column
            # whose name mentions score/prob whose values fall in [0, 1], else
            # the first such numeric column.
            score_col = None
            for prefer_named in (True, False):
                for col in df.columns:
                    lc = str(col).strip().lower()
                    if prefer_named and not ('score' in lc or 'prob' in lc):
                        continue
                    vals = pd.to_numeric(df[col], errors='coerce')
                    if vals.notna().mean() > 0.8 and vals.dropna().between(0, 1).mean() > 0.8:
                        score_col = col
                        break
                if score_col is not None:
                    break
            if score_col is None:
                raise ValueError(f'no hemolytic score column in {list(df.columns)}')

            scores = pd.to_numeric(df[score_col], errors='coerce').tolist()
            if len(scores) != len(order):
                raise ValueError(
                    f'hemopi2 returned {len(scores)} rows for {len(order)} sequences')
            return {
                seq: (round(float(min(max(s, 0.0), 1.0)), 4)
                      if pd.notna(s) else fallback[seq])
                for seq, s in zip(order, scores)
            }
    except Exception as e:
        sys.stderr.write(
            f'  Warning: HemoPI2 failed ({type(e).__name__}: {e})'
            ' — using biophysical hemo proxy fallback\n'
        )
        return fallback


def macrel_score_batch(sequences):
    """
    Score sequences for AMP activity (MACREL) and hemolysis (HemoPI2).
    Returns dict {sequence: (amp_probability, hemolytic_probability)}.

    AMP probability: from MACREL (falls back to calculate_samp for sequences
    outside 10-100 AA range or on MACREL failure).
    Hemolytic probability: from HemoPI2 (falls back to calculate_hemo_proxy if
    HemoPI2 is unavailable). HemoPI2 replaces MACREL's own hemolytic output,
    which saturated to 0.000 on the evolved sequences and gave no gradient.
    """
    import subprocess, tempfile, glob as _glob
    if not sequences:
        return {}
    # Hemolytic probability for the whole batch (one HemoPI2 call).
    hemo_scores = hemopi2_score_batch(sequences)
    fallback = {seq: (calculate_samp(seq), hemo_scores[seq]) for seq in sequences}
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, 'batch.faa')
            with open(fasta_path, 'w') as fh:
                for i, seq in enumerate(sequences):
                    fh.write(f'>seq{i}\n{seq}\n')
            try:
                clean_env = {k: v for k, v in os.environ.items()
                             if not k.startswith('Malloc')}
                result = subprocess.run(
                    # --keep-negatives is REQUIRED here. Macrel's default is to
                    # emit only the sequences it classifies as AMP:
                    #     if not keep_negatives:
                    #         final = final.query('is_AMP').drop('is_AMP', axis=1)
                    # (macrel/AMP_predict.py). Without the flag every sequence
                    # below the classifier's decision threshold is absent from
                    # the output, gets no score, and silently falls through to
                    # the biophysical surrogate. That is fatal for an
                    # optimisation target: the whole point is to read a
                    # probability for candidates that are NOT yet antimicrobial
                    # and let selection climb the gradient.
                    ['macrel', 'peptides', '--fasta', fasta_path,
                     '--output', tmpdir, '--force', '--keep-negatives'],
                    capture_output=True, text=True, timeout=300,
                    env=clean_env
                )
            except FileNotFoundError:
                sys.stderr.write(
                    '  Warning: macrel not installed — '
                    'conda install -c bioconda macrel\n'
                    '           falling back to biophysical s_amp\n'
                )
                return fallback
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip())
            pred_files = _glob.glob(os.path.join(tmpdir, '*.prediction.gz'))
            if not pred_files:
                raise FileNotFoundError('no .prediction.gz in macrel output')
            df = pd.read_csv(pred_files[0], sep='\t', comment='#',
                             compression='gzip')
            seq_map = {
                str(row['Sequence']): float(row['AMP_probability'])
                for _, row in df.iterrows()
            }
            # MACREL is defined for 10-100 residues and silently omits anything
            # outside that window, in which case the biophysical proxy is
            # substituted below. That substitution changes what the objective
            # measures, so it must not pass unnoticed in a long run. Reported
            # here only; no behaviour or log-schema change.
            missing = [s for s in sequences if s not in seq_map]
            if missing:
                lens = sorted({len(s) for s in missing})
                sys.stderr.write(
                    f'  Warning: MACREL returned no score for {len(missing)}/'
                    f'{len(sequences)} sequence(s), lengths {lens[0]}-{lens[-1]}. '
                    'Biophysical proxy substituted for those. Check that '
                    '--keep-negatives is passed and that lengths are within '
                    'MACREL\'s 10-100 residue range.\n'
                )
            return {
                seq: (seq_map.get(seq, fallback[seq][0]),
                      hemo_scores[seq])
                for seq in sequences
            }
    except Exception as e:
        sys.stderr.write(
            f'  Warning: MACREL failed ({type(e).__name__}: {e})'
            ' — using fallback\n'
        )
        return fallback
