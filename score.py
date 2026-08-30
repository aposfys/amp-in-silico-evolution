import numpy as np
import sys, os
import tempfile
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
 


def _cbeta(r, res_ca, res_cb, res_n, res_c):
    """
    CB position for residue r: the real atom when present, otherwise the
    standard virtual CB built from N, CA and C (Yang et al. 2020, PNAS 117,
    1496-1503), falling back to CA if the backbone is incomplete.
    """
    if r in res_cb:
        return res_cb[r]
    ca = res_ca[r]
    if r in res_n and r in res_c:
        b = ca - res_n[r]
        c = res_c[r] - ca
        a = np.cross(b, c)
        return -0.58273431 * a + 0.56802827 * b - 0.54067466 * c + ca
    return ca


def get_nconts(pdb_txt, chain="A", distance_cutoff=6.0, plddt_cutoff=0):
    """
    Number of tertiary contacts, per Eq. 5 of Sahakyan et al. 2025 (PNAS 122,
    e2509015122):

        c(i,j) = 1 if  d(Cbeta_i, Cbeta_j) < 6 A
                   and |i - j| > 5 residues
                   and pLDDT of both residues > 50 (0.5 on a 0-1 scale)

    The |i - j| > 5 rule is what makes this a measure of TERTIARY compactness:
    it excludes the i -> i+4 alpha-helical register, so a long helix does not
    score as a compact fold. That distinction matters biologically as well as
    arithmetically for this fork's target class: membrane-active antimicrobial
    peptides are disordered in aqueous solution and fold to an amphipathic helix
    only on contact with a lipid bilayer (Cardoso et al. 2021, Biophys Rev 13,
    35-69; Zhang et al. 2021, Mil Med Res 8, 48), so they never occupy the
    compact state this term scores. Expect ~1.0 for a 25-30 aa helix; that is
    correct, not a regression. Upstream PFES loops `range(i + 4, n_atoms)`, which
    admits |i - j| = 4 and 5 and therefore counts exactly those helical
    contacts; at its 6 A cutoff the geometry mostly masks the error (Ca(i)-
    Ca(i+4) is ~6.2 A in an ideal helix), but at a wider cutoff it dominates
    the count. Measured on this project's own output at 8 A, 95-98% of counted
    contacts were |i - j| = 4. Hence: Cbeta, 6 A, and range(i + 6, ...).

    ESM3 emits BACKBONE-ONLY structures (N, CA, C, O -- no CB at all), so a
    literal "read the CB atom" implementation would silently degrade to CA for
    every residue. Where CB is absent it is therefore reconstructed from the
    backbone with the standard virtual-CB formula (Yang et al. 2020, trRosetta),
    which is exact for ideal tetrahedral geometry. A real CB is used when the
    model provides one (e.g. ESMFold full-atom output). Glycine gets the same
    virtual CB, which is the usual convention for contact definitions.

    pLDDT is read scale-agnostically. ESMFold writes 0-100 into the B-factor
    column, ESM3 writes 0-1; `plddt_cutoff` may be given on either scale and is
    converted to match. Returns (n_contacts, mean_plddt) with mean_plddt on a
    0-1 scale, averaged PER RESIDUE (not per atom, which would weight
    tryptophan 14x against glycine 4x).
    """

    # One entry per residue: prefer CB, fall back to CA (glycine, or a model
    # that omits side chains).
    res_ca, res_cb, res_plddt = {}, {}, {}
    res_n, res_c = {}, {}
    for line in pdb_txt.splitlines():
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        if line[21] != chain:
            continue
        try:
            res_seq = int(line[22:26].strip())
            b_factor = float(line[60:66].strip())
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        atom_name = line[12:16].strip()
        if atom_name == 'CA':
            res_ca[res_seq] = xyz
            res_plddt[res_seq] = b_factor
        elif atom_name == 'CB':
            res_cb[res_seq] = xyz
            res_plddt.setdefault(res_seq, b_factor)
        elif atom_name == 'N':
            res_n[res_seq] = xyz
        elif atom_name == 'C':
            res_c[res_seq] = xyz

    if not res_plddt:
        return (0, 0.0)

    # Detect the pLDDT scale from the data and normalise both it and the cutoff
    # to 0-1, so a caller passing 50 or 0.5 gets the same behaviour.
    raw = np.array(list(res_plddt.values()), dtype=float)
    scale = 100.0 if raw.max() > 1.0 else 1.0
    plddt = {k: v / scale for k, v in res_plddt.items()}
    cutoff = plddt_cutoff / 100.0 if plddt_cutoff > 1.0 else float(plddt_cutoff)
    mean_plddt = float(np.mean(list(plddt.values())))

    keep = sorted(r for r in res_ca if plddt.get(r, 0.0) > cutoff)
    if len(keep) < 2:
        return (0, round(mean_plddt, 3))

    coords = np.array([_cbeta(r, res_ca, res_cb, res_n, res_c) for r in keep])
    resnum = np.array(keep)
    d = np.linalg.norm(coords[:, None] - coords, axis=2)
    sep = np.abs(resnum[:, None] - resnum)          # true residue separation,
                                                    # not index separation: a
                                                    # pLDDT-filtered gap must
                                                    # not shrink |i - j|.
    n_conts = int(np.count_nonzero(np.triu((d < distance_cutoff) & (sep > 5), k=1)))
    return (n_conts, round(mean_plddt, 3))


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
    Interface contacts and interface pLDDT between two chains.

    Per Sahakyan et al. 2025 (Methods): ipLDDT is "pLDDT on the interface of
    interacting proteins, that includes residues between chains where C-beta
    atoms are closer than 6 A", and iCD is "the number of contacts between
    chains calculated in the same way".

    Rewritten during the Eq. 5 audit. The previous version:

      * was named and documented for C-beta but filtered `== "CA"`, the same
        misnomer get_nconts carried;
      * returned `round(iPLDDT * 0.01, 2)`, which assumes ESMFold's 0-100
        B-factor scale. ESM3 writes 0-1, so this returned ~0.0098 instead of
        ~0.98 -- and unlike the copy in get_nconts, this value is NOT
        discarded: it multiplies straight into the dimer score, scaling every
        inter_chain fitness down by about 100x;
      * averaged chain A's interface pLDDT only, ignoring chain B, though the
        paper defines the interface as residues *between* chains.

    Cbeta is reconstructed from the backbone where absent (ESM3 emits
    backbone-only structures) and the real atom used where present, matching
    get_nconts. Returns (n_interface_contacts, mean_interface_pLDDT) with
    pLDDT on a 0-1 scale.
    """
    ca, cb, n_at, c_at, bf, ch = {}, {}, {}, {}, {}, {}
    for line in pdb_txt.splitlines():
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        try:
            key = (line[21], int(line[22:26].strip()))
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            b = float(line[60:66].strip())
        except ValueError:
            continue
        name = line[12:16].strip()
        if name == 'CA':
            ca[key] = xyz; bf[key] = b; ch[key] = line[21]
        elif name == 'CB':
            cb[key] = xyz; bf.setdefault(key, b)
        elif name == 'N':
            n_at[key] = xyz
        elif name == 'C':
            c_at[key] = xyz

    if not bf:
        return (1, 0.01)
    raw = np.array(list(bf.values()), dtype=float)
    scale = 100.0 if raw.max() > 1.0 else 1.0
    cut = plddt_cutoff / 100.0 if plddt_cutoff > 1.0 else float(plddt_cutoff)

    def side(chain):
        keys = [k for k in ca if k[0] == chain and bf.get(k, 0.0) / scale > cut]
        keys.sort(key=lambda k: k[1])
        return keys

    A, B = side(chainA), side(chainB)
    if not A or not B:
        return (1, 0.01)

    Ac = np.array([_cbeta(k, ca, cb, n_at, c_at) for k in A])
    Bc = np.array([_cbeta(k, ca, cb, n_at, c_at) for k in B])
    d = np.linalg.norm(Ac[:, None] - Bc, axis=2)
    mask = d <= distance_cutoff
    n_contacts = int(mask.sum())
    if n_contacts == 0:
        return (1, 0.01)

    ai, bi = np.where(mask)
    # Interface pLDDT over residues on BOTH sides of the interface.
    vals = ([bf[A[i]] / scale for i in np.unique(ai)] +
            [bf[B[j]] / scale for j in np.unique(bi)])
    return (n_contacts, round(float(np.mean(vals)), 3))


def iplddt_all_atom(pdb_txt, chainA='A', chainB='B', distance_cutoff=6.0,):
    iplddt_all_atom = 'not ready yet'
    return iplddt_all_atom



if __name__ == '__main__':
    input_pdb_path = str(sys.argv[1])

    if os.path.isfile(input_pdb_path):
        with open(input_pdb_path, 'r') as file:
            pdb_txt = file.read()

        print("inner contacts, plddt:" + str(get_nconts(pdb_txt, "A", 6.0, 0.5)))
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


HEMOPI2_MODEL = int(os.environ.get('PFES_HEMO_MODEL', '3'))


def hemopi2_score_batch(sequences, model=None):
    """
    Predict hemolytic probability with HemoPI2 (Rathore et al. 2025).
    Returns dict {sequence: hemo_probability in [0, 1]}, higher = more hemolytic.

    HemoPI2 replaces the old biophysical calculate_hemo_proxy. It is an ML model
    trained on 1,926 experimentally validated hemolytic peptides (AUROC ~0.92)
    and produces a real, varying signal on the cationic amphipathic peptides PFES
    evolves — unlike MACREL's hemolytic RF, which saturates to ~0.000 out of its
    training distribution and gave no selection gradient.

    Install:  pip install hemopi2
    CLI:      hemopi2_classification -i in.fa -o out.csv -m {1=RF,2=RF+MERCI,3=ESM2,4=ESM+MERCI}

    **model=3 (ESM2-t6), not the package default.** Measured on magainin 2,
    melittin and poly-alanine, three of the four models are unusable:

        m=1 RF          12.000  44.535   2.710   -- not in [0,1] at all
        m=2 RF+MERCI    MERCI -1.0 for every sequence, Hybrid pinned to 1.0
        m=3 ESM2-t6      0.229   0.764   0.638   -- in range, correctly ordered
        m=4 ESM2+MERCI  MERCI -1.0 for every sequence, Hybrid pinned to 0.0

    HemoPI2 documents a decision threshold of 0.46 (RF) or 0.55 (ESM), so the
    score is a probability by construction; models 1 and 2 emit uncalibrated
    values, and against a 0.46 threshold every peptide including poly-alanine is
    then called hemolytic. MERCI returns its -1.0 sentinel for everything, which
    makes both hybrids degenerate in opposite directions -- model 4, the package
    default, calls melittin non-hemolytic. Only the ESM2 path is sound here.

    This is the MACREL onnxruntime failure in a second tool: a classifier
    returning raw decision values rather than calibrated probabilities, with the
    threshold comparison then marking everything positive. The [0,1] range check
    below is what catches it; clamping instead would have logged 1.0000 for every
    candidate of every generation and looked entirely plausible.

    Override with PFES_HEMO_MODEL if a future release repairs the others.

    Falls back per-sequence to calculate_hemo_proxy if HemoPI2 is not installed
    or fails.
    """
    if model is None:
        model = HEMOPI2_MODEL
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


def hemopi2_agrees():
    """True if HemoPI2 is installed, in range, and oriented the right way.

    Range alone is not enough. Model 4 returns a Hybrid Score of exactly 0.0 for
    every sequence -- perfectly inside [0,1] and perfectly useless, and it calls
    melittin non-hemolytic. So check the ordering against peptides whose relative
    hemolytic activity is not in question: melittin is the bee-venom lytic
    peptide and must outrank magainin 2, which is the standard example of an
    antimicrobial peptide with minimal hemolysis.

    Mirrors macrel_inproc_agrees(). Used by preflight.sh; a run that fails this
    is logging something other than what hemo_prob claims to be.
    """
    MELITTIN = 'GIGAVLKVLTTGLPALISWIKRKRQQ'
    MAGAININ = 'GIGKFLHSAKKFGKAFVGEIMNS'
    d = hemopi2_score_batch([MAGAININ, MELITTIN])
    mag, mel = d.get(MAGAININ), d.get(MELITTIN)
    if mag is None or mel is None:
        return False, 'HemoPI2 returned no score'
    # An exact match against the surrogate means HemoPI2 never answered.
    if (abs(mag - calculate_hemo_proxy(MAGAININ)) < 1e-9
            and abs(mel - calculate_hemo_proxy(MELITTIN)) < 1e-9):
        return False, 'returned the biophysical surrogate, not HemoPI2'
    if mel <= mag:
        return False, (f'orientation wrong: melittin {mel:.3f} does not exceed '
                       f'magainin 2 {mag:.3f}')
    if mel == mag:
        return False, 'degenerate: both peptides scored identically'
    return True, f'magainin 2 {mag:.3f} < melittin {mel:.3f}'


def _macrel_subprocess(sequences):
    """The documented CLI path. Returns {normalised_sequence: amp_probability}."""
    import subprocess, glob as _glob
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, 'batch.faa')
        with open(fasta_path, 'w') as fh:
            for i, seq in enumerate(sequences):
                fh.write(f'>seq{i}\n{seq}\n')
        clean_env = {k: v for k, v in os.environ.items()
                     if not k.startswith('Malloc')}
        # --keep-negatives is REQUIRED: macrel's default emits only sequences it
        # classifies as AMP, so every candidate below the decision threshold
        # would be absent and fall through to the surrogate. Fatal for an
        # optimisation target, whose whole point is to score candidates that are
        # NOT yet active so selection has a gradient to climb.
        result = subprocess.run(
            ['macrel', 'peptides', '--fasta', fasta_path,
             '--output', tmpdir, '--force', '--keep-negatives'],
            capture_output=True, text=True, timeout=300, env=clean_env)
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip())
        pred_files = _glob.glob(os.path.join(tmpdir, '*.prediction.gz'))
        if not pred_files:
            raise FileNotFoundError('no .prediction.gz in macrel output')
        df = pd.read_csv(pred_files[0], sep='\t', comment='#', compression='gzip')
        return {str(r['Sequence']): float(r['AMP_probability'])
                for _, r in df.iterrows()}


# --------------------------------------------------------------------------- #
# In-process classifiers
#
# MACREL and HemoPI2 are invoked once per generation as subprocesses, and each
# invocation re-loads its model. Measured against the v2 series -- which ran the
# identical pipeline with neither installed, falling back to in-process proxies
# -- the two subprocesses account for ~60 s of an 88.6 s generation at pop 100,
# against 5.8 s for ESM3 folding. They are the pipeline's bottleneck, and the
# cost is per-generation reload rather than inference, so it does not amortise
# over a larger population.
#
# Loading each model ONCE per process removes that. The risk is that it means
# calling library internals rather than the documented CLI, and this project has
# already been burned twice by a classifier that silently returned the wrong
# numbers. So the fast path is self-validating: on first use it is checked
# against the subprocess on a known peptide, and any disagreement beyond
# tolerance disables it permanently for the life of the process. It can never
# quietly produce different numbers than the path it replaces.
#
# Set PFES_NO_INPROC=1 to force the subprocess path.
# --------------------------------------------------------------------------- #

_INPROC = {}          # name -> callable | False (tried and unavailable)
_INPROC_TOL = 1e-3


def _macrel_inproc():
    """
    A callable(list[str]) -> {normalised_sequence: amp_prob}, or None.

    The saving is NOT simply "call the library instead of the CLI".
    macrel.AMP_predict.predict() opens and deserialises both ONNX models on
    every invocation:

        with gzip.open(model1, 'rb') as f:
            model1 = rt.InferenceSession(f.read(), ...)

    so calling it in-process would reload the models exactly as often as the
    subprocess does. The reload IS the cost. This holds one InferenceSession
    for the lifetime of the process and runs inference directly against it.

    Only the AMP model is loaded. macrel's second model is its hemolytic
    classifier, which this pipeline does not use -- it saturated to 0.000 on
    evolved cationic sequences and HemoPI2 replaced it -- so loading it would
    double the startup for a column that is thrown away.

    Feature construction and the ['AMP'] extraction deliberately mirror
    AMP_predict.predict() exactly, including reading features from column 2
    onward. That extraction is the line onnxruntime 1.26 broke by changing the
    shape of output_probability, which is why requirements.txt pins
    onnxruntime<=1.25.1; behaving identically to macrel here means the pin
    protects this path too.
    """
    if 'macrel' in _INPROC:
        return _INPROC['macrel'] or None
    if os.environ.get('PFES_NO_INPROC') == '1':
        _INPROC['macrel'] = False
        return None
    try:
        import gzip as _gz
        import onnxruntime as _rt
        import macrel
        from macrel import AMP_features as _af

        _model = os.path.join(os.path.dirname(macrel.__file__),
                              'data', 'models', 'AMP.onnx.gz')
        with _gz.open(_model, 'rb') as fh:
            _sess = _rt.InferenceSession(
                fh.read(), providers=["CPUExecutionProvider"])

        def _run(seqs):
            with tempfile.TemporaryDirectory() as td:
                fa = os.path.join(td, 'b.faa')
                with open(fa, 'w') as fh:
                    for i, q in enumerate(seqs):
                        fh.write(f'>s{i}\n{q}\n')
                data = _af.fasta_features(fa)
            if not len(data):
                return {}
            feats = data.iloc[:, 2:].values.astype(np.float32)
            [prob] = _sess.run(['output_probability'],
                               {'input_features': feats})
            amp = [float(x['AMP']) for x in prob]
            # data['sequence'] is already macrel-normalised (leading M stripped)
            return dict(zip(data['sequence'].astype(str), amp))

        _INPROC['macrel'] = _run
        return _run
    except Exception as e:
        sys.stderr.write(f'  in-process MACREL unavailable ({type(e).__name__}: {e}); '
                         'using the subprocess path\n')
        _INPROC['macrel'] = False
        return None


def macrel_inproc_agrees(probe='GIGKFLHSAKKFGKAFVGEIMNS'):
    """
    True if the in-process path reproduces the subprocess path on `probe`.

    Called by preflight.sh and on first use. A False here disables the fast
    path rather than accepting its numbers: an objective that silently changes
    identity is the failure this whole guard exists to prevent.
    """
    fast = _macrel_inproc()
    if fast is None:
        return False
    try:
        a = fast([probe]).get(_macrel_key(probe))
        b = _macrel_subprocess([probe]).get(_macrel_key(probe))
    except Exception:
        _INPROC['macrel'] = False
        return False
    if a is None or b is None or abs(a - b) > _INPROC_TOL:
        sys.stderr.write(f'  in-process MACREL disagrees with the subprocess '
                         f'({a} vs {b}); disabling the fast path\n')
        _INPROC['macrel'] = False
        return False
    return True


def _macrel_key(s):
    """MACREL reports the normalised sequence: leading M and trailing * removed."""
    if s and s[0] == 'M':
        s = s[1:]
    if s and s[-1] == '*':
        s = s[:-1]
    return s


def macrel_score_batch_src(sequences, with_hemo=True):
    """
    As macrel_score_batch, but returns 3-tuples
    {sequence: (amp_probability, hemolytic_probability, source)} where
    source is 'macrel' or 'proxy'.

    The distinction matters because MACREL is defined for 10-100 residues
    and the biophysical calculate_samp surrogate is substituted silently
    outside that window. On an arm with no length term the population
    drifts past 100 aa, so the amp_prob column changes identity mid-run
    with only a line on stderr to say so. Recording the source per row
    puts that in the data instead of the console, where it survives into
    the analysis.

    Score sequences for AMP activity (MACREL) and hemolysis (HemoPI2).
    Returns dict {sequence: (amp_probability, hemolytic_probability, source)}.

    AMP probability: from MACREL (falls back to calculate_samp for sequences
    outside 10-100 AA range or on MACREL failure).
    Hemolytic probability: from HemoPI2 (falls back to calculate_hemo_proxy if
    HemoPI2 is unavailable). HemoPI2 replaces MACREL's own hemolytic output,
    which saturated to 0.000 on the evolved sequences and gave no gradient.
    """
    import subprocess, glob as _glob
    if not sequences:
        return {}
    # Hemolytic probability for the whole batch (one HemoPI2 call).
    #
    # with_hemo=False skips it entirely and records NaN. HemoPI2 is measured at
    # 57.7 s per generation on the production node, 65% of an 88.6 s generation,
    # for a column that never enters the fitness -- 16 h of a 1000-generation
    # run. NaN rather than 0.0 because 0.0 is a valid hemolysis probability and
    # cannot mark "not measured"; PFES_SKIP_HEMO keeps writing 0.0 for
    # compatibility with the existing series.
    hemo_scores = (hemopi2_score_batch(sequences) if with_hemo
                   else {seq: float('nan') for seq in sequences})

    # Fast path: model loaded once per process rather than once per generation.
    # Used only after it has been shown to reproduce the subprocess on a known
    # peptide; see macrel_inproc_agrees().
    if _INPROC.get('macrel') is None:
        macrel_inproc_agrees()
    _fast = _INPROC.get('macrel') or None
    if _fast is not None:
        try:
            seq_map = _fast(sequences)
            miss = [q for q in sequences if _macrel_key(q) not in seq_map]
            if miss:
                lens = sorted({len(q) for q in miss})
                sys.stderr.write(
                    f'  Warning: MACREL returned no score for {len(miss)}/'
                    f'{len(sequences)} sequence(s), lengths {lens[0]}-{lens[-1]}. '
                    'Biophysical proxy substituted for those.\n')
            return {
                q: (seq_map.get(_macrel_key(q), calculate_samp(q)),
                    hemo_scores[q],
                    'macrel' if _macrel_key(q) in seq_map else 'proxy')
                for q in sequences
            }
        except Exception as e:
            sys.stderr.write(f'  in-process MACREL failed mid-run '
                             f'({type(e).__name__}: {e}); reverting to the '
                             'subprocess for the rest of this process\n')
            _INPROC['macrel'] = False
    fallback = {seq: (calculate_samp(seq), hemo_scores[seq], 'proxy')
                for seq in sequences}
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
            # MACREL reports the NORMALISED sequence, not the one it was given.
            # macrel_features.normalize_seq strips a leading methionine (it was
            # written for ORFs, where M is the start codon) and a trailing stop
            # character:
            #     if seq[0] == 'M':  seq = seq[1:]
            #     if seq[-1] == '*': seq = seq[:-1]
            # Looking up the original string therefore misses every candidate
            # that begins with M, which is about one in twenty, and each of
            # those silently falls through to the biophysical surrogate.
            def macrel_key(s):
                if s and s[0] == 'M':
                    s = s[1:]
                if s and s[-1] == '*':
                    s = s[:-1]
                return s

            # MACREL is defined for 10-100 residues and silently omits anything
            # outside that window, in which case the biophysical proxy is
            # substituted below. That substitution changes what the objective
            # measures, so it must not pass unnoticed in a long run. Reported
            # here only; no behaviour or log-schema change.
            missing = [s for s in sequences if macrel_key(s) not in seq_map]
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
                seq: (seq_map.get(macrel_key(seq), fallback[seq][0]),
                      hemo_scores[seq],
                      'macrel' if macrel_key(seq) in seq_map else 'proxy')
                for seq in sequences
            }
    except Exception as e:
        sys.stderr.write(
            f'  Warning: MACREL failed ({type(e).__name__}: {e})'
            ' — using fallback\n'
        )
        return fallback


def macrel_score_batch(sequences, with_hemo=True):
    """
    Backward-compatible view of macrel_score_batch_src: returns
    {sequence: (amp_probability, hemolytic_probability)}, dropping the scorer
    provenance. Kept because several analysis entry points unpack 2-tuples.
    New code that cares which scorer answered should call
    macrel_score_batch_src directly.
    """
    return {k: (v[0], v[1])
            for k, v in macrel_score_batch_src(sequences, with_hemo).items()}
