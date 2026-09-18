#!/usr/bin/env python3
"""Which residue is in which secondary structure?

    python3 ss_table.py structure.pdb [--chain A] [--offset 0] [--ranges]
    python3 ss_table.py 3tq6                       # fetch from the PDBe API

Runs DSSP (Kabsch & Sander 1983) on the coordinates and prints, per residue,
whether it sits in a helix, a strand or a loop.

DSSP does not look at dihedral angles: it computes an electrostatic energy for
every backbone N-H...O=C hydrogen bond and reads the PATTERN of those bonds.
An i -> i+4 ladder is an alpha-helix, i -> i+3 a 3-10 helix, bonds between
distant segments a beta-sheet.

`--offset` is added to the file's numbering, for structures renumbered from 1.
"""
import argparse, json, sys, urllib.request
from pathlib import Path

AD = {"H": "alpha-helix", "G": "3-10 helix", "I": "pi helix", "P": "PPII",
      "E": "beta-strand", "B": "beta-bridge", "T": "turn", "S": "bend",
      "-": "loop/coil", " ": "loop/coil"}

ap = argparse.ArgumentParser()
ap.add_argument("target", help="PDB file, or a 4-character PDB code to fetch")
ap.add_argument("--chain", default=None)
ap.add_argument("--offset", type=int, default=0)
ap.add_argument("--ranges", action="store_true", help="only print the ranges")
a = ap.parse_args()


def from_pdbe(code):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/{code.lower()}"
    with urllib.request.urlopen(url, timeout=30) as f:
        d = json.load(f)[code.lower()]["molecules"]
    print(f"PDBe deposited assignment for {code.upper()}")
    for m in d:
        for ch in m["chains"]:
            if a.chain and ch["chain_id"] != a.chain:
                continue
            print(f"\n  chain {ch['chain_id']}")
            for kind, lst in ch["secondary_structure"].items():
                for s in lst:
                    ad = {"helices": "helix", "strands": "strand"}.get(kind, kind)
                    print(f"    {s['start']['residue_number']:>4d}-"
                          f"{s['end']['residue_number']:<4d} {ad}")


if not Path(a.target).exists():
    if len(a.target) == 4:
        from_pdbe(a.target)
        sys.exit()
    sys.exit(f"no such file: {a.target}")

try:
    import warnings; warnings.filterwarnings("ignore")
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP
except ImportError:
    sys.exit("needs MDAnalysis >= 2.4  (pip install MDAnalysis)")

u = mda.Universe(a.target)
sel = "protein" + (f" and chainID {a.chain}" if a.chain else "")
prot = u.select_atoms(sel)
if not len(prot):
    sys.exit(f"no atoms matched '{sel}'")

code = DSSP(prot).run().results.dssp[0]
resid = prot.residues.resids + a.offset

if not a.ranges:
    print(f"{len(code)} residues | H alpha-helix  E beta-strand  - loop\n")
    for b in range(0, len(code), 60):
        print(f"  {resid[b]:4d} {''.join(code[b:b+60])} {resid[min(b+59, len(code)-1)]:4d}")
    print()

# collapse to ranges, treating every helix type as helix
def kaba(c):
    return "H" if c in "HGI" else ("E" if c in "EB" else "-")

start, prev, son, out = resid[0], kaba(code[0]), resid[0], []
for r, c in zip(resid[1:], code[1:]):
    k = kaba(c)
    if k != prev:
        out.append((start, son, prev))
        start, prev = r, k
    son = r
out.append((start, resid[-1], prev))

print("ranges:")
for b, e, k in out:
    n = e - b + 1
    print(f"  {b:4d}-{e:<4d} ({n:3d})  {AD[k]}")

from collections import Counter
c = Counter(kaba(x) for x in code)
print("\ncomposition: " + "  ".join(
    f"{AD[k]} {v} ({100*v/len(code):.0f}%)" for k, v in c.most_common()))
