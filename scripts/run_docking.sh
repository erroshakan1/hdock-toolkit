#!/bin/bash
# HDOCKlite ile bir reseptoru bir veya daha fazla ligandla yerlestir.
#
# KULLANIM:
#   bash scripts/hdock_kos.sh <etiket> <reseptor.pdb> [ligand1.pdb ligand2.pdb ...]
#
# Ligand verilmezse bu projenin uc Abeta konformeri kullanilir:
#   analiz/S1/konformer_1.pdb  2  3
#
# Ornekler:
#   bash scripts/hdock_kos.sh S4 yapilar/hazir/siteY_duz.pdb
#   bash scripts/hdock_kos.sh S6 yapilar/hazir/S6_reseptor.pdb
#   bash scripts/hdock_kos.sh test ~/.local/opt/HDOCKlite/1CGI_r_b.pdb \
#                                  ~/.local/opt/HDOCKlite/1CGI_l_b.pdb
#
# Reseptor BUYUK olan, ligand KUCUK olan molekul olmali.
set -e
ETIKET=$1; RESEPTOR=$2; shift 2 2>/dev/null || true
[ -z "$ETIKET" ] || [ -z "$RESEPTOR" ] && {
    echo "KULLANIM: $0 <etiket> <reseptor.pdb> [ligand.pdb ...]"; exit 1; }
[ -f "$RESEPTOR" ] || { echo "HATA: $RESEPTOR yok."; exit 1; }
command -v hdock >/dev/null || { echo "HATA: hdock PATH'te yok."; exit 1; }

KOK=$(pwd)
BURADA=$(cd "$(dirname "$0")" && pwd)
LIGANDLAR=("$@")
if [ ${#LIGANDLAR[@]} -eq 0 ]; then
    echo "HATA: en az bir ligand PDB'si ver."; exit 1
fi
for L in "${LIGANDLAR[@]}"; do
    [ -f "$L" ] || { echo "HATA: ligand $L yok."; exit 1; }
done

CIK=${HDOCK_OUT:-$KOK/docking}/$ETIKET
mkdir -p "$CIK"

python3 "$BURADA/prepare_pdb.py" "$RESEPTOR" "$CIK/reseptor.pdb"

# Ligand zinciri reseptorunkiyle CAKISMAMALI. Cakisirsa model dosyasinda ikisi
# birlesir, ligand ayirt edilemez ve kumeleme analizi sahte sonuc verir.
KULLANILAN=$(awk '/^ATOM/{print substr($0,22,1)}' "$CIK/reseptor.pdb" | sort -u | tr -d '\n')
LZ=""
for h in L B E F G H I J K M N; do
    case "$KULLANILAN" in *"$h"*) ;; *) LZ=$h; break;; esac
done
[ -z "$LZ" ] && { echo "HATA: bos zincir harfi bulunamadi ($KULLANILAN)"; exit 1; }
echo "reseptor zincirleri: $KULLANILAN  ->  ligand zinciri: $LZ"

i=0
for L in "${LIGANDLAR[@]}"; do
    i=$((i+1))
    python3 "$BURADA/prepare_pdb.py" "$L" "$CIK/ligand_$i.pdb" $LZ
    if [ ! -f "$CIK/k$i.out" ]; then
        echo ">>> $ETIKET ligand $i : docking basliyor  $(date +%H:%M:%S)"
        ( cd "$CIK" && hdock reseptor.pdb ligand_$i.pdb -out k$i.out > k$i.log 2>&1 )
    fi
    ( cd "$CIK" && createpl k$i.out top10_k$i.pdb -nmax 10 -rmsd 5.0 -complex -models \
        > createpl_k$i.log 2>&1 )
    ( cd "$CIK" && for m in model_*.pdb; do [ -f "$m" ] && mv "$m" "k${i}_$m"; done )
    echo "    bitti: $(ls "$CIK"/k${i}_model_*.pdb 2>/dev/null | wc -l) model"
done

echo
echo "Cikti: $CIK   (ligand zinciri $LZ)"
echo "Skorlar     :  python3 $BURADA/score_summary.py $CIK"
echo "Poz secimi  :  python3 $BURADA/select_pose.py $CIK"
