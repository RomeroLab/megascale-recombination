#!/usr/bin/env bash
# Map one FASTQ in a preconfigured Python/Biopython + MAFFT environment.
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
python3 "$script_dir/fastq2chimera.py" "$1"
