# Validation: re-docking a known complex (1CGI)

HDOCKlite ships two PDB files as its own example. The `_b` suffix means **bound**:
these are the two partners as they appear in the crystal structure of the complex.
Docking them back together should recover the crystal pose — if the pipeline is
wired correctly.

```bash
bash scripts/run_docking.sh demo ~/.local/opt/HDOCKlite/1CGI_r_b.pdb \
                                 ~/.local/opt/HDOCKlite/1CGI_l_b.pdb
python3 scripts/score_summary.py docking/demo
```

```
pose   score  confidence  RMSD-to-input  distance-to-pose-1
   1  -445.23       1.00            0.5                 0.0
   2  -249.91       0.88            6.0                 4.3
   3  -243.15       0.87           18.3                 3.1
...
  8/9 poses within 10 A of the best
```

**Pose 1 lands 0.5 A from the crystal structure.** That is the pipeline working.

It also demonstrates the one case where the RMSD column means something: here the
input *was* the answer, so 0.5 A reads as "found it". In a real prediction the
ligand starts somewhere arbitrary and the column carries no information.

## Where does it bind?

```bash
python3 scripts/contact_map.py docking/demo --top 8
```

```
chain  residue   poses
    E      192   10/10
    E       99   10/10
    E      146    9/10
    E       39    9/10
    E       97    9/10
    E       57    9/10
    E      215    9/10
    E       41    9/10

hot spot (>=60% of poses): [39, 40, 41, 57, 58, 96, 97, 99, 143, 146, 149,
                            150, 192, 193, 195, 214, 215, 216, 217, 218, 219]
```

Residues **57**, **195**, **102** are the catalytic triad of a serine protease;
**192**, **214–219** line the substrate groove and oxyanion hole. The inhibitor
docks into the active site, which is where it binds in reality.

Two independent checks — geometric (0.5 A RMSD) and biochemical (the right
residues) — agreeing is what validation looks like. Run this after any change to
the scripts.
