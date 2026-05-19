Please use [AMES](https://github.com/sahakyanhk/ames). It is an updated version of PFES that supports more simulations.
#

# PFES: protein fold evolution simulation

Code for [In silico evolution of globular protein folds from random sequences](https://www.pnas.org/doi/10.1073/pnas.2509015122)


### Installation and usage examples 
```
wget https://github.com/sahakyanhk/PFES/archive/refs/heads/amps-esm3.zip -O pfes-amps-esm3.zip; unzip pfes-amps-esm3.zip

python pfes-amps-esm3/pfes.py -h

# A HuggingFace token with access to esm3-sm-open-v1 is required:
export HF_TOKEN=your_huggingface_token

# run a simulation evolving AMPs from random peptides and analyse results
python pfes-amps-esm3/pfes.py  -ng 100 -ps 50 -sm weak -em single_chain -iseq random --random_seq_len 24 -o pfes_test_random 
python pfes-amps-esm3/visual_pfes.py -l pfes_test_random/progress.log -s pfes_test_random/structures/ -o pfes_test_random/ 

# run a simulation starting from polyalanine 
python pfes-amps-esm3/pfes.py  -ng 100 -ps 50 -sm weak -em single_chain -iseq AAAAAAAAAAAAAAAAAAAAAAAA -o pfes_test_polyA

```
<p align="center">
  <img src="examples/example2.gif" width="350" height="350"/>
  <img src="examples/Summary.png" width="350" height="350"/>
</p>

### Extended data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14061036.svg)](https://doi.org/10.5281/zenodo.14061036)


### Hardware requirements 
This branch uses [ESM3](https://github.com/evolutionaryscale/esm) (`esm3-sm-open-v1`) for structure prediction. Access is gated on HuggingFace — set `HF_TOKEN` before running.

Supported devices: NVIDIA CUDA GPUs (V100, A100), Apple Silicon (MPS), and CPU.

PFES was originally tested on Rocky Linux 8.7 (Green Obsidian) with NVIDIA Tesla V100 and A100 GPUs.

### AMP score (`s_amp`)
The `amps-esm3` branch adds a biophysical AMP fitness term (`s_amp`) that steers evolution toward antimicrobial peptide-like sequences. It is a geometric mean of three components:
- **S_charge** — net positive charge reward (R, K favored over D, E)
- **S_hydro** — hydrophobic ratio reward (optimal ~50%)
- **S_amphi** — amphipathicity reward (Eisenberg scale, window N=11, angle 100°)

`s_amp` is multiplied directly into the fitness score and logged as a column in `progress.log`.




