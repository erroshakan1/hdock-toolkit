#!/usr/bin/env python3
"""Are the contacted residues on ONE FACE of a helix, or scattered?

    python3 helix_wheel.py <model.pdb> --chain A --range 58-76 \
            --contacts 60,63,64,67,68,71,72 [--offset 43]

Why this matters
----------------
An alpha-helix makes a full turn every 3.6 residues, so consecutive residues
are rotated 360/3.6 = 100 degrees from one another. Residues i, i+3, i+4 and
i+7 therefore end up pointing the SAME way; i+1 and i+2 point elsewhere.

So a contact list like 60, 63, 64, 67, 68, 71 (spacing 3,1,3,1,3) is the
signature of a binding partner lying against one face of a helix, i.e. a real
surface patch. A list of consecutive residues is geometrically impossible for
a surface contact, and a scattered list means no defined site.

This script fits the helix axis, projects every CA onto the plane normal to it,
and reports the angular position of each residue plus how tightly the contacted
ones cluster.
"""
import argparse, sys
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("model", type=Path)
ap.add_argument("--chain", default="A")
ap.add_argument("--range", required=True, help="e.g. 58-76 (in OUTPUT numbering)")
ap.add_argument("--contacts", required=True, help="comma-separated, OUTPUT numbering")
ap.add_argument("--offset", type=int, default=0, help="added to file numbering")
a = ap.parse_args()

lo, hi = (int(x) for x in a.range.split("-"))
contacts = {int(x) for x in a.contacts.split(",")}

ca = {}
for l in open(a.model):
    if l.startswith("ATOM") and l[21] == a.chain and l[12:16].strip() == "CA":
        ca[int(l[22:26]) + a.offset] = np.array(
            [float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)])
idx = [i for i in range(lo, hi + 1) if i in ca]
if len(idx) < 5:
    sys.exit(f"only {len(idx)} CA atoms in range - check --chain/--range/--offset")

P = np.array([ca[i] for i in idx])
centre = P.mean(0)
axis = np.linalg.svd(P - centre)[2][0]
if np.dot(axis, P[-1] - P[0]) < 0:
    axis = -axis

u = P - centre
u = u - np.outer(u @ axis, axis)               # drop the along-axis component
ref = u[0] / np.linalg.norm(u[0])
perp = np.cross(axis, ref)
ang = np.degrees(np.arctan2(u @ perp, u @ ref)) % 360

print(f"{'residue':>8s} {'angle':>7s}   contact")
for i, r in enumerate(idx):
    print(f"{r:8d} {ang[i]:6.0f}°   {'<<<' if r in contacts else ''}")

t = np.radians([ang[i] for i, r in enumerate(idx) if r in contacts])
if len(t) < 2:
    sys.exit("\nneed at least 2 contacted residues in range")
v = np.exp(1j * t)
spread = np.degrees(np.sqrt(-2 * np.log(abs(v.mean()))))   # circular std
# rotation per residue: wrap each step into (0,360) then take the MEDIAN.
# A circular mean of the steps is wrong here - steps near 0/360 cancel.
steps = np.diff(ang) % 360
gaps = np.diff(idx)                      # skip over missing residues
step = float(np.median(steps[gaps == 1] / 1.0)) if (gaps == 1).any() else float("nan")

print(f"\ncontacted residues span {spread:.0f}° (circular std)")
print(f"   < 60°  one tight face      -> defined binding patch")
print(f"  60-120° one broad face      -> plausible patch")
print(f"   >120°  spread all around   -> no defined site, likely noise")
print(f"\nmean rotation per residue: {step:.0f}°"
      f"   (alpha-helix 100°, 3-10 helix 120°, beta-strand ~180°)")
print("NOTE: verify the secondary structure with DSSP - a poorly chosen range")
print("      that straddles a loop will give a meaningless axis fit.")
