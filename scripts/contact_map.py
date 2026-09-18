#!/usr/bin/env python3
"""Per-residue contact frequency across all docking poses.

    python3 contact_map.py <docking_dir> [--cutoff 4.5] [--offset 0] [--chain X]

For every receptor residue, counts in how many of the N poses it lies within
`cutoff` A of the ligand. A residue contacted in many independent poses is a
real hot spot; one contacted once is noise.

`--offset` is added to residue numbers on output, for when the structure was
renumbered from 1 but you want to report the numbering used in the literature
(e.g. GROMACS renumbered Ser44 to 1 -> --offset 43).
"""
import argparse, re, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("docking_dir", type=Path)
ap.add_argument("--cutoff", type=float, default=4.5)
ap.add_argument("--offset", type=int, default=0)
ap.add_argument("--chain", help="restrict to this receptor chain")
ap.add_argument("--top", type=int, default=0, help="print only the N hottest")
a = ap.parse_args()


def ligand_chain(d):
    g = sorted(d.glob("ligand_*.pdb"))
    if not g:
        sys.exit("no ligand_*.pdb found - was this produced by run_docking.sh?")
    for l in open(g[0]):
        if l.startswith("ATOM"):
            return l[21]
    sys.exit("ligand file has no ATOM records")


LIG = ligand_chain(a.docking_dir)
models = sorted(a.docking_dir.glob("k*_model_*.pdb"),
                key=lambda p: (p.name.split("_")[0], int(re.search(r"_(\d+)\.pdb$", p.name).group(1))))
if not models:
    sys.exit(f"no models in {a.docking_dir}")

count, n = {}, 0
for p in models:
    lig, rec, rid = [], [], []
    for l in open(p):
        if not l.startswith("ATOM"):
            continue
        xyz = [float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)]
        if l[21] == LIG:
            lig.append(xyz)
        elif a.chain is None or l[21] == a.chain:
            rec.append(xyz)
            rid.append((l[21], int(l[22:26])))
    if not lig or not rec:
        continue
    n += 1
    t = cKDTree(rec)
    # count each residue ONCE per pose, not once per atom pair
    for r in {rid[j] for i in t.query_ball_point(lig, a.cutoff) for j in i}:
        count[r] = count.get(r, 0) + 1

print(f"{len(models)} models, {n} usable | ligand chain {LIG} | cutoff {a.cutoff} A")
print(f"{'chain':>5s} {'residue':>8s} {'poses':>7s}  {'frequency':<22s}")
items = sorted(count.items(), key=lambda kv: -kv[1])
if a.top:
    items = items[:a.top]
else:
    items = sorted(items, key=lambda kv: (kv[0][0], kv[0][1]))
top = max(count.values())
for (ch, r), v in items:
    bar = "#" * int(22 * v / top)
    print(f"{ch:>5s} {r + a.offset:8d} {v:4d}/{n}  {bar}")

hot = sorted(r + a.offset for (ch, r), v in count.items() if v >= 0.6 * n)
print(f"\nhot spot (>=60% of poses): {hot}")
