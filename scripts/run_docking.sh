#!/bin/bash
# Dock one receptor against one or more ligands with HDOCKlite.
#
#   bash run_docking.sh <label> <receptor.pdb> <ligand.pdb> [ligand2.pdb ...]
#
# The receptor is the LARGER molecule: HDOCK builds a grid around it and slides
# the ligand through that grid.
#
# Output goes to ./docking/<label>/ ; override the parent with HDOCK_OUT.
#
# Examples:
#   bash run_docking.sh trial receptor.pdb ligand.pdb
#   bash run_docking.sh ensemble dna.pdb conf1.pdb conf2.pdb conf3.pdb
set -e

LABEL=$1; RECEPTOR=$2; shift 2 2>/dev/null || true
[ -z "$LABEL" ] || [ -z "$RECEPTOR" ] && {
    echo "usage: $0 <label> <receptor.pdb> <ligand.pdb> [ligand2.pdb ...]"; exit 1; }
[ -f "$RECEPTOR" ] || { echo "error: no such receptor: $RECEPTOR"; exit 1; }
command -v hdock >/dev/null || { echo "error: hdock not on PATH"; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)
LIGANDS=("$@")
[ ${#LIGANDS[@]} -eq 0 ] && { echo "error: give at least one ligand PDB"; exit 1; }
for L in "${LIGANDS[@]}"; do
    [ -f "$L" ] || { echo "error: no such ligand: $L"; exit 1; }
done

OUT=${HDOCK_OUT:-$(pwd)/docking}/$LABEL
mkdir -p "$OUT"

python3 "$HERE/prepare_pdb.py" "$RECEPTOR" "$OUT/receptor.pdb"

# The ligand chain must not collide with a receptor chain. If it does, the two
# merge in the model PDB: the ligand can no longer be separated, and every
# downstream analysis silently gives wrong answers, because receptor atoms that
# never move end up pinning what you think is the ligand centroid.
USED=$(awk '/^ATOM/{print substr($0,22,1)}' "$OUT/receptor.pdb" | sort -u | tr -d '\n')
LCHAIN=""
for c in L B E F G H I J K M N; do
    case "$USED" in *"$c"*) ;; *) LCHAIN=$c; break;; esac
done
[ -z "$LCHAIN" ] && { echo "error: no free chain letter (receptor uses $USED)"; exit 1; }
echo "receptor chains: $USED  ->  ligand chain: $LCHAIN"

i=0
for L in "${LIGANDS[@]}"; do
    i=$((i+1))
    python3 "$HERE/prepare_pdb.py" "$L" "$OUT/ligand_$i.pdb" $LCHAIN
    if [ ! -f "$OUT/k$i.out" ]; then
        echo ">>> $LABEL ligand $i: docking  $(date +%H:%M:%S)"
        ( cd "$OUT" && hdock receptor.pdb ligand_$i.pdb -out k$i.out > k$i.log 2>&1 )
    fi
    ( cd "$OUT" && createpl k$i.out top10_k$i.pdb -nmax 10 -rmsd 5.0 -complex -models \
        > createpl_k$i.log 2>&1 )
    ( cd "$OUT" && for m in model_*.pdb; do [ -f "$m" ] && mv "$m" "k${i}_$m"; done )
    echo "    $(ls "$OUT"/k${i}_model_*.pdb 2>/dev/null | wc -l) models written"
done

echo
echo "output: $OUT   (ligand chain $LCHAIN)"
echo "  scores : python3 $HERE/score_summary.py $OUT"
echo "  poses  : python3 $HERE/select_pose.py $OUT"
