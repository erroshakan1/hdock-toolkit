#!/usr/bin/env python3
"""HDOCK pozlarini temas parmak iziyle kumele ve baslangic yapisini sec.

  python3 scripts/hdock_poz_sec.py analiz/docking/S4 [cikis.pdb]

Uc konformerin 30 pozunu BIRLIKTE degerlendirir. Ayni DNA bolgesine oturan
pozlar ayni kumeye duser; en kalabalik kumenin en iyi skorlu uyesi secilir.
"""
import sys, re
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def ligand_zinciri(kok):
    import glob
    g = sorted(Path(kok).glob("ligand_*.pdb"))
    if not g:
        return "A"
    for l in open(g[0]):
        if l.startswith("ATOM"):
            return l[21]
    return "A"


KOK = Path(sys.argv[1])
CIKIS = Path(sys.argv[2]) if len(sys.argv) > 2 else None
LIGZ = ligand_zinciri(KOK)
TEMAS = 4.5     # A
ORTAK = 0.5     # parmak izi ortusme esigi (Jaccard)


def oku(p):
    skor = None
    lig, res, res_id = [], [], []
    for l in open(p):
        if l.startswith("REMARK Score:"):
            skor = float(l.split(":")[1])
        elif l.startswith("ATOM"):
            xyz = [float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)]
            if l[21] == LIGZ:
                lig.append(xyz)
            else:
                res.append(xyz)
                res_id.append((l[21], int(l[22:26])))
    return skor, np.array(lig), np.array(res), res_id


modeller = sorted(KOK.glob("k?_model_*.pdb"),
                  key=lambda p: (p.name[1], int(re.search(r"_(\d+)\.pdb$", p.name).group(1))))
if not modeller:
    sys.exit(f"HATA: {KOK} icinde model bulunamadi.")

kayit = []
for p in modeller:
    s, lig, res, rid = oku(p)
    agac = cKDTree(res)
    dokunan = set()
    for i in agac.query_ball_point(lig, TEMAS):
        for j in i:
            dokunan.add(rid[j])
    kayit.append({"ad": p.name, "skor": s, "izi": dokunan, "n": len(dokunan)})

# Jaccard ortusmesine gore acgozlu kumeleme
kalan = sorted(range(len(kayit)), key=lambda i: kayit[i]["skor"])
kumeler = []
while kalan:
    c = kalan.pop(0)
    grup = [c]
    for i in list(kalan):
        a, b = kayit[c]["izi"], kayit[i]["izi"]
        if a and b and len(a & b) / len(a | b) >= ORTAK:
            grup.append(i)
            kalan.remove(i)
    kumeler.append(grup)

print("=" * 74)
print(f"POZ KUMELEME : {KOK}   ({len(kayit)} poz, temas esigi {TEMAS} A)")
print("=" * 74)
for k, g in enumerate(kumeler[:6], 1):
    en = kayit[g[0]]
    zin = {}
    for ch, no in sorted(en["izi"], key=lambda x: (x[0], x[1])):
        zin.setdefault(ch, []).append(no)
    yer = " · ".join(f"{ch}{min(v)}-{max(v)} ({len(v)} nt)" for ch, v in zin.items())
    print(f"\nkume {k}: {len(g):2d} poz | en iyi skor {en['skor']:8.2f} | {en['ad']}")
    print(f"   temas: {yer}")
    if len(g) > 1:
        print(f"   uyeler: {', '.join(kayit[i]['ad'].replace('_model','').replace('.pdb','') for i in g)}")

büyük = max(kumeler, key=len)
sec = kayit[büyük[0]]
print("\n" + "=" * 74)
print(f"SECIM: {sec['ad']}   ({len(büyük)}/{len(kayit)} pozluk en kalabalik kume,"
      f" skor {sec['skor']:.2f})")
if len(büyük) == 1:
    print("UYARI: hicbir poz tekrarlanmiyor, kume yok. Bu, belirgin bir baglanma")
    print("       bolgesi olmadigi anlamina gelir; secim buyuk olcude keyfi.")
print("=" * 74)

if CIKIS:
    CIKIS.write_text((KOK / sec["ad"]).read_text())
    print(f"\n{CIKIS} yazildi.")
