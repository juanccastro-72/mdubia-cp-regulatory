#!/usr/bin/env bash
# End-to-end reproduction of the Myrciaria dubia chloroplast regulatory-element analysis.
# Usage: bash scripts/run_all.sh [EMAIL] [NCBI_API_KEY]
set -euo pipefail
EMAIL="${1:-you@example.com}"
APIKEY="${2:-}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
export PYTHONPATH="$HERE/src:${PYTHONPATH:-}"

EXCEL="data/ChloroplastGenomes_ncbi_dataset_Myrtaceae_29102026.xlsx"
INV="data/Chloroplast_regulatory_elements_inventory.xlsx"
THREADS="${THREADS:-4}"

echo "[1/4] Core pipeline (downloads plastomes, predicts & compares elements)"
python src/cpreg/cpreg_pipeline.py --excel "$EXCEL" --inventory "$INV" \
    --genomes-dir genomes --outdir results --threads "$THREADS" \
    --email "$EMAIL" ${APIKEY:+--api-key "$APIKEY"}

echo "[2/4] Self-repeats and GC windows for the Circos figures"
python src/cpreg/make_circos_inputs.py --genomes-dir genomes --out circos_work

echo "[3/4] Figures (main Fig1-7, supplementary FigS1-5)"
python src/cpreg/cpreg_figures.py --results results --genomes-dir genomes --out figures

echo "[3b] Synteny + Circos (Fig8-12)"
python src/cpreg/cpreg_synteny.py       --results results --genomes-dir genomes --out figures
python src/cpreg/cpreg_circos_detail.py --results results --genomes-dir genomes --out figures --which all

echo "[4/4] Done. Tables in results/, figures in figures/."
