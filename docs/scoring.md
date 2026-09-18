# Reading HDOCK scores

## Docking score

Unitless, in ITScore units. **More negative is better.** It is a knowledge-based
statistical potential, not an energy: it cannot be converted to kJ/mol and must
never be reported as a binding free energy. Use it for **ranking only**.

## Confidence score

The HDOCK server derives a calibrated probability from the docking score:

```
confidence = 1 / (1 + exp(0.02 * (score + 150)))
```

`score_summary.py` computes this locally.

| Confidence | Reading |
|---|---|
| > 0.7 | binding very likely |
| 0.5 – 0.7 | binding possible |
| < 0.5 | binding unlikely |

Calibrated on protein–protein benchmarks. For protein–DNA/RNA it is indicative,
not evidence.

What the scale looks like in practice:

| System | Best score | Confidence |
|---|---|---|
| 1CGI, bound crystal structures re-docked | −445 | **1.00** |
| Disordered peptide on a protein–DNA complex | −207 | 0.76 |
| Disordered peptide on free B-DNA | −179 | 0.64 |

A genuine, tight complex looks like the first row. Anything in the 0.4–0.7 band
is a weak or non-specific association — worth simulating, not worth claiming.

## Ligand RMSD

The distance from the pose to the **input** ligand coordinates — not a quality
measure. Model 1 does not come out at 0.

It is meaningful only when the input *was* a known complex: then a low RMSD means
docking recovered the known pose. That is exactly the 1CGI validation case
(`examples/1cgi.md`), where model 1 lands 0.5 A from the crystal pose.

When the ligand was placed arbitrarily — a conformer pulled out of an MD
trajectory, say — ignore the column.

## Pose dispersion — the number that actually matters

One good score can be noise. The same site found repeatedly, from several
independent ligand conformers, cannot be.

`score_summary.py` reports how many of the top 10 poses sit within 10 A of the
best one. `select_pose.py` goes further and clusters every pose by **contact
fingerprint**: the set of receptor residues within 4.5 A of the ligand, compared
between poses by Jaccard overlap. Poses that touch the same residues cluster
together regardless of how the ligand is rotated.

Take the representative of the **largest cluster**, not the single best score.

## A caution about geometry

Long or rod-shaped receptors bias contacts toward the middle simply because the
middle has more neighbouring surface. If your hot spot sits at the geometric
centre, check whether that is chemistry or arithmetic — a randomized-sequence or
shuffled-residue control settles it.
