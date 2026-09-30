#!/usr/bin/env bash
# All-vs-all structural comparison with Foldseek, following the Domainome
# protocol: build a database of all structures, then search every structure
# against it exhaustively with no E-value cut-off.
#
# Output: results/foldseek/foldseek_easy_allvsall  (TSV, no header)
#   query  target  alntmscore  qtmscore  ttmscore  prob
set -euo pipefail
cd "$(dirname "$0")/.."

STRUCT_DIR=${MP_DATA:-data}/structures
OUT=${MP_RESULTS:-results}/foldseek
FOLDSEEK=${FOLDSEEK:-foldseek}
mkdir -p "$OUT"

"$FOLDSEEK" createdb "$STRUCT_DIR" "$OUT/foldseek_domainome_db"
"$FOLDSEEK" easy-search "$STRUCT_DIR"/*.pdb "$OUT/foldseek_domainome_db" \
    "$OUT/foldseek_easy_allvsall" "$OUT/tmp" \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,prob" \
    --exhaustive-search 1 -e inf
rm -rf "$OUT/tmp"
echo "hits: $(wc -l < "$OUT/foldseek_easy_allvsall")"
