#!/usr/bin/env python3
"""PDB'yi HDOCK icin temizle: hidrojenleri at, zincir kimligi ve element sutunu yaz.

  python3 scripts/pdb_hdock_hazirla.py giris.pdb cikis.pdb [zincir]

Zincir verilmezse dosyadaki mevcut zincir korunur; bosa A yazilir.
"""
import sys
from pathlib import Path

giris, cikis = Path(sys.argv[1]), Path(sys.argv[2])
zorla = sys.argv[3] if len(sys.argv) > 3 else None

def element(satir):
    e = satir[76:78].strip()
    if e:
        return e
    ad = satir[12:16]
    # PDB kurali: 13. sutun doluysa element iki harfli olabilir
    s = ad[0] if ad[0] != " " else ad[1]
    return s.upper()

n = 0
atlanan = 0
with open(cikis, "w") as f:
    for l in open(giris):
        if l[:6] not in ("ATOM  ", "HETATM"):
            continue
        e = element(l)
        if e == "H":
            atlanan += 1
            continue
        n += 1
        zin = zorla if zorla else (l[21] if l[21] != " " else "A")
        yeni = ("ATOM  " + f"{n:5d}" + l[11:21] + zin + l[22:72]
                + "    " + f"{e:>2s}" + "  ")
        f.write(yeni.rstrip() + "\n")
    f.write("END\n")

print(f"{giris.name} -> {cikis.name} : {n} agir atom yazildi, {atlanan} hidrojen atildi")
