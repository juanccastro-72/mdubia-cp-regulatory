# mdubia-cp-regulatory

Reproducible, statistically calibrated **in silico identification of chloroplast
regulatory elements** in *Myrciaria dubia* (camu-camu; GenBank **PX251231.1**) and
comparative analysis across **121 Myrtaceae** plastomes plus **12 reference**
plastomes carrying experimentally validated elements.

This repository contains the full analysis software (v2 workflow), a one-command
run script, a conda environment specification, and offline regression tests that
accompany the article:

> Castro JC et al. *In silico identification of plastid regulatory elements in the
> chloroplast genome of Myrciaria dubia (camu-camu) and comparative analysis across
> Myrtaceae.* Plant Molecular Biology (2026).

> **All predictions produced by this pipeline are in silico hypotheses** that
> require experimental validation (TSS mapping, reporter assays, transplastomics).

---

## What the pipeline does

1. **Reference library** — maps experimentally validated plastid regulatory
   elements (promoters, 5′UTRs, Shine–Dalgarno motifs, 3′UTRs/terminators, the
   intercistronic expression element, the PPR10 site) onto concrete plastome loci.
2. **Homology** — tests each element against *M. dubia* by synteny-aware BLASTn and
   anchored global alignment versus an empirical null (BH-FDR).
3. **De novo prediction** — PEP (−35/−10 PSSM) and NEP (YRTA/GAA-box) promoters,
   Shine–Dalgarno motifs (anti-SD hybridisation energy), translation-initiation-region
   structure and 3′ stem-loops, each against dinucleotide-shuffle null models with
   negative and positive controls.
4. **Comparative genomics** — annotation-independent orthologue mapping across the
   121 Myrtaceae, motif-vs-flank conservation, stem-loop covariation, and a
   maximum-likelihood CDS phylogeny.
5. **Synteny** — whole-plastome alignment (MUMmer4) and Circos visualisations
   (pyCirclize): family-wide, tribe-representative and reference-genome comparisons.
6. **Evidence integration** into explicit confidence tiers, plus publication
   figures (editable SVG + 600 dpi JPEG) and result tables.

## Repository layout

```
src/cpreg/
  cpreg_pipeline.py        core analysis (homology, de novo prediction, comparative, phylogeny)
  cpreg_style.py           shared figure style (Times New Roman, SVG export)
  cpreg_figures.py         main/supplementary figures (Fig1–7, FigS1–5)
  cpreg_synteny.py         MUMmer4 synteny + family-wide Circos (Fig8, Fig9, FigS6)
  make_circos_inputs.py    self-repeats + GC windows for the detailed Circos
  cpreg_circos_detail.py   detailed focal-genome Circos (Fig10–12)
scripts/run_all.sh         end-to-end reproduction
data/                      input tables (Myrtaceae NCBI dataset; validated-element inventory)
tests/                     offline regression tests (pytest)
environment.yml            conda environment (Python 3.13 + tools)
requirements.txt           pip fallback
```

## Installation

```bash
conda env create -f environment.yml
conda activate mdubia-cp
```

This installs Python 3.13 and the required Python packages plus the external tools
**BLAST+**, **MAFFT**, **FastTree** and **MUMmer4** (`nucmer`, `show-coords`,
`delta-filter`). With `pip`, install `requirements.txt` and provide those four tools
separately.

## Usage

```bash
# full reproduction (downloads plastomes from NCBI; give your e-mail, optionally an API key)
bash scripts/run_all.sh you@domain.org YOUR_NCBI_API_KEY
```

or step by step:

```bash
export PYTHONPATH=src/cpreg
python src/cpreg/cpreg_pipeline.py \
    --excel data/ChloroplastGenomes_ncbi_dataset_Myrtaceae_29102026.xlsx \
    --inventory data/Chloroplast_regulatory_elements_inventory.xlsx \
    --genomes-dir genomes --outdir results --threads 4 --email you@domain.org
python src/cpreg/make_circos_inputs.py   --genomes-dir genomes --out circos_work
python src/cpreg/cpreg_figures.py        --results results --genomes-dir genomes --out figures
python src/cpreg/cpreg_synteny.py        --results results --genomes-dir genomes --out figures
python src/cpreg/cpreg_circos_detail.py  --results results --genomes-dir genomes --out figures --which all
```

Outputs: result tables in `results/tables/` (TSV) and `results/*.xlsx`, figures in
`figures/` (SVG + JPEG), the ML tree (`results/myrtaceae_cds_fasttree.nwk`), and a
reproducibility record (`results/run_metadata.json`, `statistics.json`).

## Tests

```bash
pip install pytest
python -m pytest -q
```

The tests are offline (no network, no downloaded genomes) and check the sequence
utilities, the dinucleotide-preserving shuffle, BH-FDR, the PEP PSSM (recovery of a
planted TTGACA…TATAAT promoter), the anti-SD hybridisation energy and the RNA-folding
helpers. Continuous integration runs them on every push (`.github/workflows/ci.yml`).

## Reproducibility

All randomised procedures use a fixed seed; software versions, parameters and input
SHA-256 checksums are written to `results/run_metadata.json`. Plastome records are
cached under `genomes/`.

## Data

`data/` holds the two project inputs: the NCBI Myrtaceae chloroplast-genome table and
the curated inventory of experimentally validated regulatory elements. Plastome
sequences are downloaded from NCBI GenBank at run time and are not redistributed here.

## Citation

See `CITATION.cff`. Please cite both the software (this repository, archived on Zenodo
upon acceptance) and the associated article.

## License

MIT — see `LICENSE`.
