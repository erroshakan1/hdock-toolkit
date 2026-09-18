#!/usr/bin/env python3
"""Docking scores, confidence scores, and how far apart the poses sit.

    python3 score_summary.py <docking_dir>

Confidence is the HDOCK server's published calibration of the docking score:

    confidence = 1 / (1 + exp(0.02 * (score + 150)))

    > 0.7   binding very likely
    0.5-0.7 binding possible
    < 0.5   binding unlikely

Calibrated on protein-protein benchmarks; for protein-DNA/RNA read it as
indicative, not as evidence.
"""
import math, re, sys
from pathlib import Path
import numpy as np


def ligand_chain(d):
    g = sorted(Path(d).glob("ligand_*.pdb"))
    if not g:
        return "A"
    for line in open(g[0]):
        if line.startswith("ATOM"):
            return line[21]
    return "A"


ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
LIG = ligand_chain(ROOT)


def confidence(score):
    return 1.0 / (1.0 + math.exp(0.02 * (score + 150.0)))


def read_model(path):
    score = rmsd = None
    xyz = []
    for line in open(path):
        if line.startswith("REMARK Score:"):
            score = float(line.split(":")[1])
        elif line.startswith("REMARK RMSD:"):
            rmsd = float(line.split(":")[1])
        elif line.startswith("ATOM") and line[21] == LIG:
            xyz.append([float(line[30 + 8 * k:38 + 8 * k]) for k in range(3)])
    return score, rmsd, np.array(xyz)


print("=" * 70)
print(f"HDOCK summary: {ROOT}    (ligand chain {LIG})")
print("=" * 70)

ligands = sorted({p.name.split("_")[0] for p in ROOT.glob("k*_model_*.pdb")})
for tag in ligands:
    models = sorted(ROOT.glob(f"{tag}_model_*.pdb"),
                    key=lambda p: int(re.search(r"_(\d+)\.pdb$", p.name).group(1)))
    if not models:
        continue
    print(f"\n--- ligand {tag[1:]} ---")
    print(f"{'pose':>4s} {'score':>9s} {'conf':>6s} {'RMSD-to-input':>14s} {'to-pose-1':>10s}")
    first = None
    for p in models:
        s, r, x = read_model(p)
        c = x.mean(0)
        if first is None:
            first = c
        n = re.search(r"_(\d+)\.pdb$", p.name).group(1)
        print(f"{n:>4s} {s:9.2f} {confidence(s):6.2f} {r:14.1f} "
              f"{np.linalg.norm(c - first):10.1f}")
    near = sum(1 for p in models[1:]
               if np.linalg.norm(read_model(p)[2].mean(0) - first) < 10.0)
    print(f"  {near}/{len(models)-1} poses within 10 A of the best")

print("""
How to read this
  score          lower (more negative) is better. Unitless ITScore, not an
                 energy - never report it as a binding free energy.
  conf           see the formula in the module docstring.
  RMSD-to-input  distance from the pose to the INPUT ligand coordinates. Only
                 meaningful when the input was a known complex; if the ligand
                 was placed arbitrarily, ignore this column entirely.
  to-pose-1      the number that matters. Poses landing in the same place mean
                 a defined site; scattered poses mean no preferred site.""")
