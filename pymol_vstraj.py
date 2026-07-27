import os 
import re
import sys
import shutil
import pandas as pd
from pymol import cmd
from moviepy.editor import ImageClip, concatenate_videoclips

# Renders the winning lineage as an animation (the visual in the upstream README).
#
# usage: python pymol_vstraj.py <outdir>/analysis/trajpdb
#
# The argument is the trajpdb/ directory, NOT structures/. Run visual_pfes.py
# first WITHOUT --notraj: backbone_traj() selects the lineage, deduplicates it
# by sequence, and copies those structures into trajpdb/ renamed by generation
# index (gndx0.pdb.gz, gndx1.pdb.gz, ...). This script relies on that ordering
# and loads gndx0.pdb.gz by name.
#
# Writes frames/ and test_vid.mp4. The README shows a .gif; converting the mp4
# is a separate step and is not part of this repository.

def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(data, key=alphanum_key)



pdb_dir = os.path.abspath(str(sys.argv[1])) # get path from the first argument

# pdb_dir = wd + '/pdb/'
# selected_pdb = wd + '/pdb2/'
frames = 'frames/'
if os.path.exists(frames):
    shutil.rmtree(frames)
os.makedirs(frames, exist_ok=True)
wdname = 'test_vid' #os.path.basename(wd)

# os.chdir(wd)



# if os.path.exists(selected_pdb):
#     shutil.rmtree(selected_pdb)
# os.makedirs(selected_pdb, exist_ok=True)


# names=['id','len', 'len_p', 'helix_p', 'ptm', 'plddt', 'nconts', 'niconts', 'score', 'seq','dssp'] #
# log = pd.read_csv('output.log', sep='\t', header=None, names=names)


print('preparing strucutres with best score')

pdbs=sorted_alphanumeric(os.listdir(pdb_dir))
# for i in range(len(log[log.id.str.contains(f'seq0')])):
#     gen1 = log[log.id.str.contains(f'gen{i}_')]
#     pdbid = gen1[gen1.score == gen1.score.max()].head(1).id.item()
#     pdbs.append(pdbid +'.pdb')

# Upstream takes every 2nd structure, which suits the long lineages it evolved
# (500-800 steps, so 250-400 frames). A peptide run produces far fewer: a
# 35-step lineage halved gives 17 frames, i.e. 1.7 s of video.
#
# Keep upstream's step of 2 wherever it still yields a watchable animation, and
# fall back to every structure only when the lineage is too short to bear
# halving. Long runs therefore render exactly as they did before.
MIN_FRAMES = 60
step = 2 if len(pdbs) // 2 >= MIN_FRAMES else 1
pdbs = pdbs[0:len(pdbs):step]

n = len(pdbs)
print(f'lineage {len(sorted_alphanumeric(os.listdir(pdb_dir)))} structures, '
      f'step {step} -> {n} frames ({n * 0.1:.1f} s)')


print(f'{n} structures will be rendered')

cmd.load(f'{pdb_dir}/gndx0.pdb.gz', 'pdb_0')
cmd.orient()

# pLDDT lives in the B-factor column, but its scale depends on the predictor:
# ESMFold writes 0-100, ESM3 writes 0-1. A hardcoded 0-100 spectrum applied to
# ESM3 output puts every residue in the bottom one percent of the range, so the
# whole chain renders one colour and the confidence gradient disappears without
# any error. Read the scale off the data instead.
_bfactors = [a.b for a in cmd.get_model('pdb_0 and chain A').atom]
PLDDT_MAX = 100.0 if (_bfactors and max(_bfactors) > 1.0) else 1.0
print(f'pLDDT scale detected: 0-{PLDDT_MAX:g} '
      f'(max B-factor {max(_bfactors) if _bfactors else float("nan"):.3f})')
# cmd.color("grey", 'pdb_0' and "chain B")
# cmd.spectrum("b", "blue_white_red",  'pdb_0' and "chain A")
#view = cmd.get_view(0) # or set view from pymol

# Upstream pasted in a camera matrix captured by hand from a PyMOL session and
# tuned for the ~100-residue globular domains it evolved: the camera sits 232 A
# back. A 26-residue peptide spans about 31 A and would occupy 13% of the frame
# width, a speck in the middle of white space.
#
# Fit the view to the molecules instead. The first and last structures of the
# lineage are framed together, so the view accommodates the whole trajectory
# including any growth (the fold-only arm runs from 24 to 68 residues). The
# result is then held FIXED for every frame, exactly as upstream does, so the
# peptide appears to evolve rather than the camera to drift.
cmd.load(f'{pdb_dir}/{pdbs[-1]}', 'pdb_last')
cmd.align('pdb_last', 'pdb_0')
cmd.orient('pdb_0 or pdb_last')
cmd.zoom('pdb_0 or pdb_last', buffer=4.0)
view = cmd.get_view()
cmd.delete('pdb_last')


print(pdbs)

i = 0
for pdb in pdbs: 
    i+=1
    cmd.load(f'{pdb_dir}/{pdb}', f'pdb_{i}')
    q=f'pdb_{i}'
    t=f'pdb_{i-1}'
    cmd.align(q, t)
    cmd.delete(t)
    cmd.set_view(view)
    cmd.bg_color(color="white")
    cmd.color("white", q and "chain B")
    cmd.spectrum("b", "red_yellow_blue",  q and "chain A", minimum=0, maximum=PLDDT_MAX)
    cmd.show("sticks",  q and "chain A")
    cmd.set("stick_radius", 0.15)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_trace_mode", 3)
    cmd.set("ray_trace_color", "black")
    cmd.set("ray_trace_gain", 25)
    #cmd.ray(800, 800, -1, 0, 0, 0)
    cmd.png(f'{frames}/frame_{i}.png', width=1800, height=1800, dpi=300, ray=0, quiet=1)
    print(f'{i} of {n}')


clips = [ImageClip(frames + m).set_duration(0.1)
           for m in sorted_alphanumeric(os.listdir(frames))]

concat_clip = concatenate_videoclips(clips,
                                     method="compose",
                                     bg_color=(255, 255, 255))
concat_clip.write_videofile(f'{wdname}.mp4', 24)

# The upstream README embeds a .gif, but the script only ever wrote an .mp4 and
# the conversion was done outside the repository. Write both, so the animation
# that goes in a README or a slide is reproducible from this script alone.
concat_clip.write_gif(f'{wdname}.gif', fps=10)
print(f'wrote {wdname}.mp4 and {wdname}.gif')

