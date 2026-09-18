# HDOCK toolkit

Scripts for running [HDOCKlite](https://github.com/huang-laboratory/HDOCKlite)
and, more importantly, for **reading its output honestly**: which poses recur,
where on the receptor they land, and whether that site is a real pocket or noise.

Built while docking an intrinsically disordered peptide onto DNA and onto a
protein–DNA complex, so it handles protein–protein, protein–DNA and
protein–RNA cases the same way.

## Install

```bash
git clone https://github.com/huang-laboratory/HDOCKlite.git ~/.local/opt/HDOCKlite
chmod +x ~/.local/opt/HDOCKlite/{hdock,createpl}
ln -s ~/.local/opt/HDOCKlite/hdock    ~/.local/bin/hdock
ln -s ~/.local/opt/HDOCKlite/createpl ~/.local/bin/createpl
pip install numpy scipy            # MDAnalysis only for the optional DSSP check
```

HDOCKlite needs `libfftw3`. Check with `ldd $(which hdock)`.

## Use

```bash
# dock: receptor first (the larger molecule), then one or more ligands
bash scripts/run_docking.sh myjob receptor.pdb ligand1.pdb ligand2.pdb

# what the scores say
python3 scripts/score_summary.py docking/myjob

# which poses recur, and where -> pick a starting structure
python3 scripts/select_pose.py docking/myjob

# per-residue contact frequency across every pose
python3 scripts/contact_map.py docking/myjob --top 15

# which residues are helix / strand / loop?  (check before trusting helix_wheel)
python3 scripts/ss_table.py receptor.pdb --chain A --ranges
python3 scripts/ss_table.py 3tq6                 # deposited assignment, PDBe API
python3 scripts/ss_table.py 3tq6 --fetch         # download and run DSSP yourself

# is the contacted set one face of a helix, or scattered?
python3 scripts/helix_wheel.py docking/myjob/k1_model_1.pdb \
        --chain A --range 56-71 --contacts 60,63,64,67,68,71
```

Output goes to `./docking/<label>/`; override with `HDOCK_OUT=/some/dir`.

## Scripts

| Script | Does |
|---|---|
| `prepare_pdb.py` | strip hydrogens, set chain ID, fill the element column |
| `run_docking.sh` | prepare inputs, pick a free ligand chain, dock, build models |
| `score_summary.py` | docking score, confidence score, pose dispersion |
| `select_pose.py` | cluster poses by **contact fingerprint**, pick a representative |
| `contact_map.py` | per-residue contact frequency over all poses |
| `helix_wheel.py` | are the contacts on one face of a helix? |
| `ss_table.py` | which residue is in which secondary structure (DSSP) |

## Four things worth knowing

**1. The ligand chain must not collide with a receptor chain.** If it does, the
two merge in the model PDB and every downstream analysis silently breaks — the
ligand centroid gets pinned by receptor atoms that never move, so all poses look
tightly clustered. `run_docking.sh` picks a free letter automatically.

**2. "Ligand RMSD" is not a quality score.** It is the distance from the pose to
the *input* ligand coordinates. Meaningful only when the input was a known
complex (see `examples/1cgi.md`); meaningless when the ligand was placed
arbitrarily.

**3. Take the biggest cluster, not the best score.** One good score can be noise;
the same site found from several independent ligand conformers cannot be.

**4. Check the residue numbering before you quote it.** A file straight from the
PDB carries author numbering — what the literature uses. Anything that has been
through `pdb2gmx` is renumbered from 1 per chain, so residue 17 in the file may
be residue 60 in the paper. `ss_table.py` prints the file's range up front and
warns when it starts at 1; pass `--offset` to restore the real numbers. Getting
this wrong is silent — nothing errors, you just discuss the wrong residues.

Scoring in detail: [`docs/scoring.md`](docs/scoring.md) ·
Validation run: [`examples/1cgi.md`](examples/1cgi.md) ·
Turkish walkthrough: [`docs/rehber-TR.md`](docs/rehber-TR.md)

## Limits

Rigid-body docking: neither partner flexes. Solvent is implicit. The score is a
knowledge-based potential, not an energy — **never report it as ΔG**. The
confidence score is calibrated on protein–protein benchmarks; for protein–DNA
it is indicative, not evidence. Use the result as a starting structure for MD,
not as a conclusion.

## Credit

The docking engine is HDOCKlite by the Huang Lab (HUST). Cite:
Yan Y, Tao H, He J, Huang S-Y. *Nat Protoc* 2020;15:1829–1852.
These scripts are MIT-licensed and unaffiliated with the Huang Lab.
