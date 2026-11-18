#!/bin/sh
# Runs every check of this directory and stores the output next to it.
# usage: sh run_all.sh <python> <repository root>
PY="$1"; REPO="$2"
cd "$(dirname "$0")" || exit 1
for f in check_*.py; do
  n="${f%.py}"
  "$PY" -B "$f" "$REPO" > "output/$n.txt" 2>&1
  echo "$f exit $?" >> output/_summary.txt
done
