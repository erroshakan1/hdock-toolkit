#!/usr/bin/env python3
"""Cluster poses by contact fingerprint and pick a starting structure.

    python3 select_pose.py <docking_dir> [output.pdb]

Every pose is reduced to the set of receptor residues within 4.5 A of the
ligand - its contact fingerprint. Two poses whose fingerprints overlap by more
than half (Jaccard) are the same binding mode, however the ligand is rotated.

The representative of the LARGEST cluster is chosen, not the single best score:
one good score can be noise, but the same site found independently from several
ligand conformers cannot be.
"""
import re, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

CONTACT = 4.5     # A
OVERLAP = 0.5     # Jaccard threshold


def ligand_chain(d):
    g = sorted(Path(d).glob("ligand_*.pdb"))
    if not g:
        return "A"
    for line in open(g[0]):
        if line.startswith("ATOM"):
            return line[21]
    return "A"


ROOT = Path(sys.argv[1])
DEST = Path(sys.argv[2]) if len(sys.argv) > 2 else None
LIG = ligand_chain(ROOT)


def read(path):
    score = None
    lig, rec, rid = [], [], []
    for line in open(path):
        if line.startswith("REMARK Score:"):
            score = float(line.split(":")[1])
        elif line.startswith("ATOM"):
            xyz = [float(line[30 + 8 * k:38 + 8 * k]) for k in range(3)]
            if line[21] == LIG:
                lig.append(xyz)
            else:
                rec.append(xyz)
                rid.append((line[21], int(line[22:26])))
    return score, np.array(lig), np.array(rec), rid


models = sorted(ROOT.glob("k*_model_*.pdb"),
                key=lambda p: (p.name.split("_")[0],
                               int(re.search(r"_(\d+)\.pdb$", p.name).group(1))))
if not models:
    sys.exit(f"no models found in {ROOT}")

poses = []
for p in models:
    s, lig, rec, rid = read(p)
    tree = cKDTree(rec)
    touched = {rid[j] for i in tree.query_ball_point(lig, CONTACT) for j in i}
    poses.append({"name": p.name, "score": s, "print": touched})

remaining = sorted(range(len(poses)), key=lambda i: poses[i]["score"])
clusters = []
while remaining:
    seed = remaining.pop(0)
    group = [seed]
    for i in list(remaining):
        a, b = poses[seed]["print"], poses[i]["print"]
        if a and b and len(a & b) / len(a | b) >= OVERLAP:
            group.append(i)
            remaining.remove(i)
    clusters.append(group)

print("=" * 74)
print(f"POSE CLUSTERING: {ROOT}   ({len(poses)} poses, {CONTACT} A contacts)")
print("=" * 74)
for n, group in enumerate(clusters[:6], 1):
    best = poses[group[0]]
    by_chain = {}
    for ch, num in sorted(best["print"]):
        by_chain.setdefault(ch, []).append(num)
    where = " | ".join(f"{ch}{min(v)}-{max(v)} ({len(v)} res)"
                       for ch, v in by_chain.items())
    print(f"\ncluster {n}: {len(group):2d} poses | best score {best['score']:8.2f}"
          f" | {best['name']}")
    print(f"   contacts: {where}")
    if len(group) > 1:
        members = ", ".join(poses[i]["name"].replace("_model", "").replace(".pdb", "")
                            for i in group)
        print(f"   members: {members}")

biggest = max(clusters, key=len)
chosen = poses[biggest[0]]
print("\n" + "=" * 74)
print(f"SELECTED: {chosen['name']}   "
      f"(largest cluster, {len(biggest)}/{len(poses)} poses, score {chosen['score']:.2f})")
if len(biggest) == 1:
    print("WARNING: no pose recurs. There is no defined binding site here and")
    print("         the choice is essentially arbitrary.")
print("=" * 74)

if DEST:
    DEST.write_text((ROOT / chosen["name"]).read_text())
    print(f"\nwrote {DEST}")
