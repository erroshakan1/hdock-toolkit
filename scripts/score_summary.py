#!/usr/bin/env python3
"""HDOCKlite sonuclarini ozetle: skor, guven skoru, pozlarin ayni bolgede toplanip
toplanmadigi.

  python3 scripts/hdock_ozet.py analiz/docking/S4
"""
import sys, math, re
from pathlib import Path
import numpy as np

def ligand_zinciri(kok):
    import glob
    g = sorted(Path(kok).glob("ligand_*.pdb"))
    if not g:
        return "A"
    for l in open(g[0]):
        if l.startswith("ATOM"):
            return l[21]
    return "A"


KOK = Path(sys.argv[1] if len(sys.argv) > 1 else "analiz/docking/S4")
LIGZ = ligand_zinciri(KOK)


def guven(skor):
    """HDOCK sunucusunun yayinladigi kalibrasyon:
       guven = 1 / (1 + exp(0.02 * (skor + 150)))"""
    return 1.0 / (1.0 + math.exp(0.02 * (skor + 150.0)))


def model_oku(p):
    skor = rmsd = None
    xyz = []
    for l in open(p):
        if l.startswith("REMARK Score:"):
            skor = float(l.split(":")[1])
        elif l.startswith("REMARK RMSD:"):
            rmsd = float(l.split(":")[1])
        elif l.startswith("ATOM") and l[21] == LIGZ:     # ligand zinciri
            xyz.append([float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)])
    return skor, rmsd, np.array(xyz)


print("=" * 70)
print(f"HDOCK ozeti : {KOK}")
print("=" * 70)

for k in (1, 2, 3):
    modeller = sorted(KOK.glob(f"k{k}_model_*.pdb"),
                      key=lambda p: int(re.search(r"_(\d+)\.pdb$", p.name).group(1)))
    if not modeller:
        continue
    print(f"\n--- konformer {k} ---")
    print(f"{'poz':>4s} {'skor':>9s} {'guven':>7s} {'girdiye RMSD':>13s} {'poz1e uzaklik':>14s}")
    ilk = None
    merkezler = []
    for p in modeller:
        s, r, x = model_oku(p)
        m = x.mean(0)
        merkezler.append(m)
        if ilk is None:
            ilk = m
        d = np.linalg.norm(m - ilk)
        print(f"{p.name.split('_')[-1][:-4]:>4s} {s:9.2f} {guven(s):7.2f}"
              f" {r:13.1f} {d:14.1f}")
    M = np.array(merkezler)
    yakin = sum(1 for m in M[1:] if np.linalg.norm(m - M[0]) < 10.0)
    print(f"  en iyi pozun 10 A'sinda toplanan diger poz sayisi: {yakin}/9")
print("""
OKUMA
  skor    : dusuk (daha negatif) = daha iyi. Birimi yok, ITScore biriminde.
  guven   : >0,7 cok muhtemel baglanma | 0,5-0,7 olasi | <0,5 zayif.
            Formul protein-protein kiyaslamalarinda kalibre edildi,
            protein-DNA icin yon gostericidir, kanit degildir.
  girdiye RMSD : pozun GIRDI ligand konumundan uzakligi. Girdi rastgele
            yerlestirildigi icin kalite olcusu DEGIL, sadece bilgi.
  poz1e uzaklik: asil bakilacak sey. Ilk 10 poz ayni bolgede toplaniyorsa
            (kucuk degerler) docking tutarli; dagilmissa belirgin bir
            baglanma bolgesi yok demektir.""")
