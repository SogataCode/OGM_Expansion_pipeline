#!/usr/bin/env bash
set -euo pipefail
usage() {
  echo 'Usage: bash expansion_distribution.sh --input RAW_DIR --output OUTPUT_DIR [--prefix NAME] --locus chr4:3934999-3953000'
}
INPUT='' OUTPUT='' PREFIX='' LOCUS=''
while (($#)); do
  case "$1" in
    --input|--output|--prefix|--locus)
      if (($# < 2)); then usage >&2; exit 2; fi
      case "$1" in
        --input) INPUT=$2;; --output) OUTPUT=$2;; --prefix) PREFIX=$2;; --locus) LOCUS=$2;;
      esac
      shift 2;;
    --help|-h) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2;;
  esac
done
if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$LOCUS" ]]; then usage >&2; exit 2; fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"
"$PYTHON" "$ROOT/ogm_pipeline.py" --input "$INPUT" --output "$OUTPUT" --prefix "$PREFIX" --locus "$LOCUS"
"$PYTHON" "$ROOT/finish_gmm.py" --output "$OUTPUT" --prefix "$PREFIX" --locus "$LOCUS"
