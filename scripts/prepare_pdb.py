#!/usr/bin/env python3
"""Clean a PDB so HDOCK can read it.

    python3 prepare_pdb.py input.pdb output.pdb [chain]

Three fixes, all of which matter:

  * strip hydrogens   - docking works on heavy-atom geometry; hydrogens are
                        wasted cycles and a parsing hazard
  * set the chain ID  - GROMACS writes column 22 blank, and without a chain ID
                        you cannot tell ligand from receptor in the output
  * fill the element  - columns 77-78; some tools read the atom type from here
                        and misinterpret the record when it is empty

If no chain is given the existing one is kept, defaulting to A. Pass a chain
explicitly for the ligand, and make sure it is one the receptor does not use.
"""
import sys
from pathlib import Path

if len(sys.argv) < 3:
    sys.exit(__doc__)

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
forced = sys.argv[3] if len(sys.argv) > 3 else None


def element_of(line):
    e = line[76:78].strip()
    if e:
        return e
    name = line[12:16]
    # PDB convention: a name starting in column 13 may be a two-letter element
    return (name[0] if name[0] != " " else name[1]).upper()


kept = dropped = 0
with open(dst, "w") as out:
    for line in open(src):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        e = element_of(line)
        if e == "H":
            dropped += 1
            continue
        kept += 1
        chain = forced if forced else (line[21] if line[21] != " " else "A")
        rebuilt = ("ATOM  " + f"{kept:5d}" + line[11:21] + chain + line[22:72]
                   + "    " + f"{e:>2s}" + "  ")
        out.write(rebuilt.rstrip() + "\n")
    out.write("END\n")

print(f"{src.name} -> {dst.name}: {kept} heavy atoms, {dropped} hydrogens dropped")
