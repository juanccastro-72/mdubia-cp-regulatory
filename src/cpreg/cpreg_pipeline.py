#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpreg_pipeline.py  (version 2.0.0)
==================================

Reproducible in silico identification of chloroplast regulatory elements in
*Myrciaria dubia* (camu-camu; GenBank PX251231.1, NCBI Nuccore 3055957251) and
comparative analysis across the Myrtaceae plastomes listed in
``ChloroplastGenomes_ncbi_dataset_Myrtaceae_29102026.xlsx``.

This is a full re-implementation of ``chloroplast_regulatory_element_pipeline.py``
(v1). The main scientific changes relative to v1 are:

1.  Experimentally validated reference library. The 48 records of
    ``Chloroplast_regulatory_elements_inventory.xlsx`` (plus additional
    literature-supported elements: the tobacco intercistronic expression element,
    PPR10 binding site, psbD light-responsive promoter, PclpP and PrpoB NEP
    promoters) are mapped to concrete loci in the plastomes of the species in
    which they were validated, and their sequences are extracted directly from
    the RefSeq records (no hand-typed sequences except the heterologous phage T7g10
    leader).
2.  Homology evidence is two-fold and statistically calibrated: (i) BLASTn with a
    positional-orthology (synteny) check and (ii) gene-anchored global alignment
    whose identity is compared with an empirical null distribution obtained by
    aligning the same element to random, non-orthologous M. dubia intergenic
    windows of identical length (empirical P, Benjamini-Hochberg FDR).
3.  De novo predictions are tested against composition-matched null models
    (Altschul-Erickson dinucleotide-preserving shuffles of each tested window,
    reproducing the complete search procedure, i.e. the maximum-over-positions
    statistic), followed by FDR control. Degenerate regexes of v1 (e.g. "TANNNT",
    "TTGNNN", "GGA") that produce thousands of hits by chance were removed.
       * PEP (sigma70-type) promoters: position-specific scoring matrices for the
         -35 (TTGACA) and -10 (TATAAT) hexamers with a 15-19 nt spacer model.
       * NEP promoters: YRTA / GAA-box composite motif enrichment, plus homology
         to experimentally mapped NEP promoters (PaccD-129, PatpB-289, PclpP-53,
         PrpoB-345, Prrn16-Pc).
       * Shine-Dalgarno (SD): RNA-RNA hybridisation free energy (ViennaRNA
         duplexfold) between the 3' end of the M. dubia 16S rRNA (anti-SD) and
         positions -20..-1 of every CDS, threshold dG <= -4.4 kcal/mol
         (Ma et al. 2002) plus shuffle-based empirical P.
       * Translation-initiation-region (TIR) structure: MFE z-scores.
       * 3' stem-loops (terminator / RNA-stability elements): sliding-window MFE
         with a null model that repeats the sliding-window minimisation on
         shuffled sequences; z-score and empirical P.
4.  Comparative genomics across all Myrtaceae plastomes is annotation-independent
    (orthologous regions are located by BLASTn of M. dubia loci), and includes
    position-level conservation of predicted motifs relative to flanking spacer
    sequence, cross-species re-scoring of each element, compensatory base-pair
    changes (covariation) in predicted stem-loops, and a maximum-likelihood
    phylogeny inferred from concatenated protein-coding genes (MAFFT + FastTree).
5.  Evidence integration into explicit tiers:
       Tier 1  homologous to an experimentally validated element (syntenic and
               significant) AND/OR significant de novo prediction conserved in
               >= 95 % of Myrtaceae;
       Tier 2  significant de novo prediction (FDR < 0.05) with partial
               conservation;
       Tier 3  nominal prediction only (exploratory).
6.  The untrained CNN of v1 was removed: no experimentally labelled
    positive/negative training set of sufficient size exists for plastid
    regulatory elements, and an untrained or overfitted network cannot provide
    calibrated evidence.
7.  Full reproducibility: fixed random seed, cached NCBI records with SHA-256
    checksums, recorded software versions and parameters (run_metadata.json).

All predictions are hypotheses (in silico evidence) and must be reported as
"putative"/"candidate" until confirmed by TSS mapping (dRNA-seq/5'RACE), reporter
assays or transplastomic analysis.

Usage
-----
python cpreg_pipeline.py \
    --excel ChloroplastGenomes_ncbi_dataset_Myrtaceae_29102026.xlsx \
    --inventory Chloroplast_regulatory_elements_inventory.xlsx \
    --target PX251231.1 --outdir results --email you@domain --threads 4

Dependencies: python>=3.9, biopython>=1.80, pandas, numpy, scipy, openpyxl,
requests, ViennaRNA (python bindings, RNA module) and the executables blastn /
makeblastdb (BLAST+ >= 2.12), mafft (>= 7.4) and FastTree (>= 2.1).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    import RNA  # ViennaRNA
except ImportError as exc:  # pragma: no cover
    sys.exit("ViennaRNA python bindings are required: pip install ViennaRNA")

__version__ = "2.0.0"

# ---------------------------------------------------------------------------
# Global parameters (all recorded in run_metadata.json)
# ---------------------------------------------------------------------------
PARAMS = dict(
    seed=20261029,
    promoter_window_max=300,       # nt upstream of start codon / 5' end scanned
    promoter_window_min=40,
    pep_spacer_min=15, pep_spacer_max=19, pep_spacer_opt=17,
    pep_spacer_penalty=0.5,        # bits per nt deviation from 17
    pssm_consensus_prob=0.70,      # P(consensus base) in -35/-10 PSSMs
    n_shuffles_motif=10000,
    n_shuffles_rna=1000,
    n_shuffles_stemloop=500,
    sd_window=(-20, -1),           # relative to first nt of start codon
    sd_dG_threshold=-4.4,          # kcal/mol, Ma et al. (2002)
    antisd_length=13,
    tir_window=(-60, 30),
    stemloop_region_max=150,
    stemloop_window=60, stemloop_step=5,
    homology_null_n=500,
    blast_evalue=1e-5,
    fdr_alpha=0.05,
    conservation_tier1=0.95,
    temperature_C=37.0,
)

TARGET_ACCESSION = "PX251231.1"
TARGET_SPECIES = "Myrciaria dubia"

REFERENCE_GENOMES = {
    "NC_001879.2": "Nicotiana tabacum",
    "NC_000932.1": "Arabidopsis thaliana",
    "NC_008590.1": "Hordeum vulgare",
    "NC_005353.1": "Chlamydomonas reinhardtii",
    "NC_037507.1": "Marchantia polymorpha",
    "NC_003119.8": "Medicago truncatula",
    "NC_009259.1": "Phaseolus vulgaris",
    "NC_008096.2": "Solanum tuberosum",
    "NC_007578.1": "Lactuca sativa",
    "NC_043882.1": "Erigeron breviscapus",
    "NC_001666.2": "Zea mays",
    "NC_001320.1": "Oryza sativa",
}

log = logging.getLogger("cpreg")

# ---------------------------------------------------------------------------
# Basic sequence utilities
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
ENC = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    ENC[ord(_b)] = _i


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def clean(s: str) -> str:
    return re.sub(r"[^ACGT]", "N", str(s).upper())


def gc(s: str) -> float:
    s = s.upper()
    n = sum(s.count(b) for b in "ACGT")
    return (s.count("G") + s.count("C")) / n if n else float("nan")


def circ_slice(seq: str, start: int, end: int) -> str:
    """0-based half-open slice on a circular sequence (start may be <0, end > len)."""
    L = len(seq)
    if end - start >= L:
        raise ValueError("slice longer than genome")
    s, e = start % L, end % L
    if start >= 0 and end <= L:
        return seq[start:end]
    if s < e:
        return seq[s:e]
    return seq[s:] + seq[:e]


def encode(s: str) -> np.ndarray:
    return ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


def dinuc_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Altschul & Erickson (1985) dinucleotide-preserving shuffle (Euler walk)."""
    if len(seq) < 4:
        return seq
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        edges[a].append(b)
    last = seq[-1]
    verts = [v for v in edges if v != last]
    while True:  # random arborescence of last-exit edges rooted at `last`
        lastedge = {v: int(rng.integers(len(edges[v]))) for v in verts}
        ok = True
        for v in verts:
            seen, u = set(), v
            while u != last:
                if u in seen:
                    ok = False
                    break
                seen.add(u)
                u = edges[u][lastedge[u]]
            if not ok:
                break
        if ok:
            break
    lists = {}
    for v, lst in edges.items():
        if v in lastedge:
            k = lastedge[v]
            others = lst[:k] + lst[k + 1:]
            rng.shuffle(others)
            lists[v] = others + [lst[k]]
        else:
            others = list(lst)
            rng.shuffle(others)
            lists[v] = others
    idx = {v: 0 for v in lists}
    out, u = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        nxt = lists[u][idx[u]]
        idx[u] += 1
        out.append(nxt)
        u = nxt
    return "".join(out)


def bh_fdr(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def empirical_p(obs: float, null: np.ndarray, greater: bool = True) -> float:
    null = np.asarray(null)
    k = (null >= obs).sum() if greater else (null <= obs).sum()
    return (k + 1) / (len(null) + 1)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# NCBI retrieval (cached)
# ---------------------------------------------------------------------------
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_genbank(accessions: list[str], outdir: Path, email=None, api_key=None,
                  batch: int = 10) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accessions if not (outdir / f"{a}.gb").exists()]
    for i in range(0, len(todo), batch):
        ids = todo[i:i + batch]
        params = dict(db="nuccore", id=",".join(ids), rettype="gbwithparts", retmode="text")
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key
        for attempt in range(5):
            try:
                r = requests.get(EUTILS, params=params, timeout=300)
                r.raise_for_status()
                break
            except Exception as exc:  # pragma: no cover
                log.warning("NCBI retry %d (%s)", attempt + 1, exc)
                time.sleep(3 * (attempt + 1))
        else:
            raise RuntimeError(f"NCBI download failed for {ids}")
        for rec in r.text.split("//\n"):
            if not rec.strip().startswith("LOCUS"):
                continue
            ver = [l for l in rec.splitlines() if l.startswith("VERSION")][0].split()[1]
            (outdir / f"{ver}.gb").write_text(rec + "//\n")
        time.sleep(0.4 if not api_key else 0.12)


# ---------------------------------------------------------------------------
# Genome model
# ---------------------------------------------------------------------------
RRNA_PRODUCTS = {"16s": "rrn16", "23s": "rrn23", "4.5s": "rrn4.5", "5s": "rrn5"}
GENE_ALIASES = {"rrns": "rrn16", "rrn16s": "rrn16", "rrnl": "rrn23", "rrn23s": "rrn23",
                "rrn4.5s": "rrn4.5", "rrn5s": "rrn5", "ycf3": "pafI", "ycf4": "pafII", "clpp1": "clpP",
                "pafi": "pafI", "pafii": "pafII", "pbf1": "psbN"}


def gene_name(feat) -> str:
    q = feat.qualifiers
    g = q.get("gene", [""])[0].strip()
    if not g and feat.type == "rRNA":
        prod = q.get("product", [""])[0].lower()
        for k, v in RRNA_PRODUCTS.items():
            if prod.startswith(k):
                return v
    if not g:
        g = q.get("locus_tag", ["unknown"])[0]
    g = g.replace(" ", "")
    return GENE_ALIASES.get(g.lower(), g)


@dataclass
class Locus:
    gene: str
    ftype: str
    strand: int
    five: int        # 0-based index of first transcribed/translated nt (start codon 1st nt)
    three: int       # 0-based index of last nt (stop codon last nt)
    length: int
    start_codon: str = ""
    copies: int = 1


class Plastome:
    """Thin wrapper around a GenBank plastome with circular-aware utilities."""

    def __init__(self, gb_path: Path):
        self.path = Path(gb_path)
        self.record = SeqIO.read(str(gb_path), "genbank")
        self.accession = self.record.id
        self.organism = self.record.annotations.get("organism", "")
        self.seq = clean(self.record.seq)
        self.L = len(self.seq)
        self._build()

    def _build(self):
        L = self.L
        self.mask = np.zeros(L, dtype=bool)
        self.strand_at = np.zeros(L, dtype=np.int8)
        loci = []
        for f in self.record.features:
            if f.type not in ("CDS", "tRNA", "rRNA"):
                continue
            g = gene_name(f)
            parts = f.location.parts
            if g.startswith("rps12") or len(parts) > 1 and len({p.strand for p in parts}) > 1:
                spans = [(int(p.start), int(p.end), p.strand) for p in parts]
            elif int(f.location.end) - int(f.location.start) > 20000:
                spans = [(int(p.start), int(p.end), p.strand) for p in parts]
            else:
                spans = [(int(f.location.start), int(f.location.end), parts[0].strand)]
            for s, e, st in spans:
                self.mask[s:e] = True
                self.strand_at[s:e] = st if st else 0
            p0, p1 = parts[0], parts[-1]
            st5 = p0.strand
            five = int(p0.start) if st5 == 1 else int(p0.end) - 1
            st3 = p1.strand
            three = int(p1.end) - 1 if st3 == 1 else int(p1.start)
            try:
                flen = len(f.location)
            except Exception:
                flen = 0
            sc = ""
            if f.type == "CDS":
                try:
                    sc = str(f.location.extract(self.record.seq)[:3]).upper()
                except Exception:
                    sc = ""
            loci.append(Locus(g, f.type, st5 if st5 else 1, five, three, flen, sc))
        # de-duplicate inverted-repeat copies: keep longest copy per (gene, type)
        best: dict[tuple, Locus] = {}
        cnt = Counter((l.gene, l.ftype) for l in loci)
        for l in loci:
            k = (l.gene, l.ftype)
            if k not in best or l.length > best[k].length:
                best[k] = l
        for k, l in best.items():
            l.copies = cnt[k]
        self.loci = sorted(best.values(), key=lambda x: x.five)
        self.all_loci = loci

    # -- windows ---------------------------------------------------------
    def upstream(self, five: int, strand: int, n: int) -> str:
        if n <= 0:
            return ""
        if strand == 1:
            return circ_slice(self.seq, five - n, five)
        return revcomp(circ_slice(self.seq, five + 1, five + 1 + n))

    def downstream(self, three: int, strand: int, n: int) -> str:
        if n <= 0:
            return ""
        if strand == 1:
            return circ_slice(self.seq, three + 1, three + 1 + n)
        return revcomp(circ_slice(self.seq, three - n, three))

    def around_start(self, five: int, strand: int, a: int, b: int) -> str:
        """Sequence from a to b (relative, start codon 1st nt = 0, b exclusive)."""
        if strand == 1:
            return circ_slice(self.seq, five + a, five + b)
        return revcomp(circ_slice(self.seq, five - b + 1, five - a + 1))

    def free_upstream(self, five: int, strand: int, limit: int) -> int:
        n = 0
        pos = five
        while n < limit:
            pos = (pos - 1) % self.L if strand == 1 else (pos + 1) % self.L
            if self.mask[pos]:
                break
            n += 1
        return n

    def free_downstream(self, three: int, strand: int, limit: int) -> int:
        n = 0
        pos = three
        while n < limit:
            pos = (pos + 1) % self.L if strand == 1 else (pos - 1) % self.L
            if self.mask[pos]:
                break
            n += 1
        return n

    def genome_coord(self, five: int, strand: int, rel: int) -> int:
        """1-based genome coordinate of a position relative to the start codon."""
        return ((five + rel) % self.L) + 1 if strand == 1 else ((five - rel) % self.L) + 1

    def intergenic_spacers(self, min_len: int = 1) -> pd.DataFrame:
        L, m = self.L, self.mask
        if m.all():
            return pd.DataFrame()
        # rotate so that index 0 is genic, to avoid a spacer spanning the origin split
        first_genic = int(np.argmax(m))
        rows = []
        i = 0
        while i < L:
            pos = (first_genic + i) % L
            if not m[pos]:
                j = i
                while j < L and not m[(first_genic + j) % L]:
                    j += 1
                s = (first_genic + i) % L
                e = (first_genic + j - 1) % L
                n = j - i
                left = self.strand_at[(s - 1) % L]
                right = self.strand_at[(e + 1) % L]
                if left == 1 and right == 1:
                    ori = "tandem(+)"
                elif left == -1 and right == -1:
                    ori = "tandem(-)"
                elif left == -1 and right == 1:
                    ori = "divergent"
                elif left == 1 and right == -1:
                    ori = "convergent"
                else:
                    ori = "undetermined"
                if n >= min_len:
                    rows.append(dict(start=s + 1, end=e + 1, length=n, orientation=ori,
                                     seq=circ_slice(self.seq, s, s + n)))
                i = j
            else:
                i += 1
        df = pd.DataFrame(rows)
        # flanking gene names
        if not df.empty:
            df["left_gene"] = [self._nearest_gene(r.start - 2) for r in df.itertuples()]
            df["right_gene"] = [self._nearest_gene(r.end) for r in df.itertuples()]
        return df

    def _nearest_gene(self, pos0: int) -> str:
        best, bd = "", 10 ** 9
        for l in self.all_loci:
            for p in (l.five, l.three):
                d = min(abs(p - pos0), self.L - abs(p - pos0))
                if d < bd:
                    bd, best = d, l.gene
        return best

    def locus(self, gene: str, ftype: str | None = None) -> Locus | None:
        for l in self.loci:
            if l.gene == gene and (ftype is None or l.ftype == ftype):
                return l
        return None


# ---------------------------------------------------------------------------
# Reference library of experimentally validated plastid regulatory elements
# ---------------------------------------------------------------------------
# Each element is defined by an anchor gene in the plastome of the species in
# which the element was characterised, the region relative to that anchor
# ('up' = upstream of the start codon / mature 5' end; 'down' = downstream of the
# stop codon / 3' end) and the window length. Windows are deliberately
# generous so that they contain the experimentally mapped TSS / processing site
# and core motifs reported in the cited studies. `inv_rows` link to the
# 0-based row index of Chloroplast_regulatory_elements_inventory.xlsx.
T7G10_LEADER = ("GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT")

REF_ELEMENTS = [
    # --- Nicotiana tabacum -------------------------------------------------
    dict(id="Nt_Prrn16", name="Prrn (rrn16 operon promoter; P1 PEP & Pc NEP)", acc="NC_001879.2", gene="rrn16", ftype="rRNA", region="up", length=250, cls="promoter", evidence="Reporter assay; transplastomic expression; TSS mapping", inv=[0], cite="Svab and Maliga 1993; Allison et al. 1996; Hajdukiewicz et al. 1997"),
    dict(id="Nt_PpsbA_5UTR", name="PpsbA promoter + psbA 5'UTR", acc="NC_001879.2", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Reporter assay; mutagenesis; transplastomic expression", inv=[1, 7], cite="Staub and Maliga 1993; Eibl et al. 1999; Zou et al. 2003"),
    dict(id="Nt_PrbcL_5UTR", name="PrbcL promoter + rbcL 5'UTR", acc="NC_001879.2", gene="rbcL", ftype="CDS", region="up", length=250, cls="promoter+5UTR", evidence="Reporter assay; transcript mapping; transplastomic expression", inv=[2, 8], cite="Eibl et al. 1999; Occhialini et al. 2019"),
    dict(id="Nt_PatpB_5UTR", name="PatpB (NEP -289/PEP -255) + atpB 5'UTR", acc="NC_001879.2", gene="atpB", ftype="CDS", region="up", length=350, cls="promoter+5UTR", evidence="Transcript mapping; in vitro translation; mutagenesis", inv=[3, 9], cite="Hajdukiewicz et al. 1997; Hirose and Sugiura 2004"),
    dict(id="Nt_PaccD", name="PaccD-129 (NEP promoter)", acc="NC_001879.2", gene="accD", ftype="CDS", region="up", length=250, cls="promoter", evidence="Transcript mapping; functional characterization", inv=[4], cite="Hajdukiewicz et al. 1997; Occhialini et al. 2019"),
    dict(id="Nt_PatpI", name="PatpI (NEP-associated promoter)", acc="NC_001879.2", gene="atpI", ftype="CDS", region="up", length=300, cls="promoter", evidence="Transcript mapping; functional characterization", inv=[5], cite="Occhialini et al. 2019; Liere et al. 2011"),
    dict(id="Nt_PpsbN_5UTR", name="PpsbN promoter + psbN 5'UTR", acc="NC_001879.2", gene="psbN", ftype="CDS", region="up", length=150, cls="promoter+5UTR", evidence="Transcript mapping; functional characterization", inv=[6, 11], cite="Occhialini et al. 2019"),
    dict(id="Nt_rps2_5UTR", name="rps2 5'UTR", acc="NC_001879.2", gene="rps2", ftype="CDS", region="up", length=150, cls="5UTR", evidence="Transcript mapping; functional characterization", inv=[10], cite="Occhialini et al. 2019"),
    dict(id="Nt_PclpP", name="PclpP-53 (NEP promoter)", acc="NC_001879.2", gene="clpP", ftype="CDS", region="up", length=200, cls="promoter", evidence="Transcript mapping (NEP-only rpoB-deletion plastids)", inv=[], cite="Hajdukiewicz et al. 1997"),
    dict(id="Nt_PrpoB", name="PrpoB-345 (type I NEP promoter)", acc="NC_001879.2", gene="rpoB", ftype="CDS", region="up", length=400, cls="promoter", evidence="In vitro transcription; mutagenesis", inv=[], cite="Liere and Maliga 1999"),
    dict(id="Nt_PpsbD_BLRP", name="psbD blue-light-responsive promoter (BLRP)", acc="NC_001879.2", gene="psbD", ftype="CDS", region="up", length=1000, cls="promoter", evidence="Transplastomic promoter deletion analysis", inv=[], cite="Allison and Maliga 1995"),
    dict(id="Nt_IEE", name="Intercistronic expression element (psbT-psbH spacer, psbH 5' end)", acc="NC_001879.2", gene="psbH", ftype="CDS", region="up", length=60, cls="IEE", evidence="Transplastomic operon expression; RNA processing", inv=[], cite="Zhou et al. 2007"),
    dict(id="Nt_TrbcL", name="TrbcL (rbcL 3'UTR)", acc="NC_001879.2", gene="rbcL", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay; transplastomic expression", inv=[12], cite="Tangphatsornruang et al. 2011; Occhialini et al. 2019"),
    dict(id="Nt_TrpoA", name="TrpoA (rpoA 3'UTR)", acc="NC_001879.2", gene="rpoA", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay; transplastomic expression", inv=[13], cite="Tangphatsornruang et al. 2011; Occhialini et al. 2019"),
    dict(id="Nt_Trps16", name="Trps16 (rps16 3'UTR)", acc="NC_001879.2", gene="rps16", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay; transplastomic expression", inv=[14], cite="Occhialini et al. 2019"),
    dict(id="Nt_TpsbA", name="TpsbA (psbA 3'UTR)", acc="NC_001879.2", gene="psbA", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay; transplastomic expression", inv=[40], cite="Staub and Maliga 1993; Ruhlman et al. 2010; Tangphatsornruang et al. 2011"),
    # --- Arabidopsis / barley: genome-wide TSS catalogues -----------------
    dict(id="At_PpsbA_TSS", name="psbA promoter region with mapped primary TSS", acc="NC_000932.1", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Terminome-seq; 5' RACE (TSS mapping)", inv=[15, 18], cite="Castandet et al. 2019"),
    dict(id="At_PrbcL_TSS", name="rbcL promoter region with mapped primary TSS", acc="NC_000932.1", gene="rbcL", ftype="CDS", region="up", length=250, cls="promoter+5UTR", evidence="Terminome-seq; 5' RACE (TSS mapping)", inv=[15, 16], cite="Castandet et al. 2019"),
    dict(id="At_Prrn16_TSS", name="rrn16 promoter region with mapped TSS", acc="NC_000932.1", gene="rrn16", ftype="rRNA", region="up", length=250, cls="promoter", evidence="Terminome-seq; 5' RACE (TSS mapping)", inv=[15, 19], cite="Castandet et al. 2019"),
    dict(id="Hv_PpsbA_TSS", name="psbA promoter region with mapped TSS", acc="NC_008590.1", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Differential RNA-seq; 5' RACE", inv=[20, 21], cite="Zhelyazkova et al. 2012"),
    dict(id="Hv_PrbcL_TSS", name="rbcL promoter region with mapped TSS", acc="NC_008590.1", gene="rbcL", ftype="CDS", region="up", length=250, cls="promoter+5UTR", evidence="Differential RNA-seq; 5' RACE", inv=[20, 21], cite="Zhelyazkova et al. 2012"),
    # --- Chlamydomonas reinhardtii -----------------------------------------
    dict(id="Cr_Prrn16", name="rrn16 (rrnS) promoter", acc="NC_005353.1", gene="rrn16", ftype="rRNA", region="up", length=250, cls="promoter", evidence="Reporter assay; transplastomic expression", inv=[22, 29], cite="Inckemann et al. 2025"),
    dict(id="Cr_psaA_5UTR", name="psaA exon-1 5'UTR", acc="NC_005353.1", gene="psaA", ftype="CDS", region="up", length=200, cls="5UTR", evidence="Reporter assay; mutagenesis", inv=[23, 29], cite="Inckemann et al. 2025"),
    dict(id="Cr_psbA_5UTR", name="psbA 5'UTR", acc="NC_005353.1", gene="psbA", ftype="CDS", region="up", length=200, cls="5UTR", evidence="Reporter assay; mutagenesis", inv=[24, 27, 29], cite="Mayfield et al. 1994; Inckemann et al. 2025"),
    dict(id="Cr_TrbcL", name="rbcL 3'UTR", acc="NC_005353.1", gene="rbcL", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay", inv=[25, 29], cite="Inckemann et al. 2025"),
    dict(id="Cr_TatpB", name="atpB 3'UTR", acc="NC_005353.1", gene="atpB", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Reporter assay", inv=[26, 29], cite="Stern et al. 1991; Inckemann et al. 2025"),
    # --- Marchantia polymorpha ---------------------------------------------
    dict(id="Mp_PrbcL_5UTR", name="MprbcL promoter + 5'UTR", acc="NC_037507.1", gene="rbcL", ftype="CDS", region="up", length=250, cls="promoter+5UTR", evidence="Transplastomic reporter assay", inv=[30, 32], cite="Sauret-Güeto et al. 2020"),
    dict(id="Mp_petB_leader", name="psbH-petB intergenic 5' leader", acc="NC_037507.1", gene="petB", ftype="CDS", region="up", length=150, cls="5UTR", evidence="Transplastomic reporter assay", inv=[31], cite="Sauret-Güeto et al. 2020"),
    # --- Heterologous elements tested in tobacco / other hosts -------------
    dict(id="Mt_PpsbA_5UTR", name="Medicago psbA promoter + 5'UTR", acc="NC_003119.8", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Transplastomic reporter assay (in tobacco)", inv=[33], cite="Gerasymenko et al. 2017"),
    dict(id="Pv_PrbcL_5UTR", name="Phaseolus rbcL promoter + 5'UTR", acc="NC_009259.1", gene="rbcL", ftype="CDS", region="up", length=250, cls="promoter+5UTR", evidence="Transplastomic reporter assay (in tobacco)", inv=[34], cite="Gerasymenko et al. 2017"),
    dict(id="St_Prrn16", name="Potato Prrn", acc="NC_008096.2", gene="rrn16", ftype="rRNA", region="up", length=250, cls="promoter", evidence="Reporter assay; transplastomic expression", inv=[35, 38], cite="Occhialini et al. 2019"),
    dict(id="St_PaccD", name="Potato PaccD", acc="NC_008096.2", gene="accD", ftype="CDS", region="up", length=250, cls="promoter", evidence="Reporter assay; transplastomic expression", inv=[36], cite="Occhialini et al. 2019"),
    dict(id="St_PpsbA_5UTR", name="Potato PpsbA + 5'UTR", acc="NC_008096.2", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Reporter assay; transplastomic expression", inv=[37], cite="Occhialini et al. 2019"),
    dict(id="Ls_PpsbA_5UTR", name="Lettuce psbA promoter + 5'UTR", acc="NC_007578.1", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Transplastomic reporter assay", inv=[39], cite="Ruhlman et al. 2010"),
    dict(id="Ls_TpsbA", name="Lettuce psbA 3'UTR", acc="NC_007578.1", gene="psbA", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Transplastomic reporter assay", inv=[39], cite="Ruhlman et al. 2010"),
    # --- Erigeron breviscapus (candidate elements) -------------------------
    dict(id="Eb_Prrn16", name="Eb-Prrn", acc="NC_043882.1", gene="rrn16", ftype="rRNA", region="up", length=250, cls="promoter", evidence="Sequence identification; candidate validation", inv=[41], cite="Yu et al. 2021"),
    dict(id="Eb_PpsbA", name="Eb-PpsbA", acc="NC_043882.1", gene="psbA", ftype="CDS", region="up", length=200, cls="promoter+5UTR", evidence="Sequence identification; candidate validation", inv=[42], cite="Yu et al. 2021"),
    dict(id="Eb_TpsbA", name="Eb-TpsbA", acc="NC_043882.1", gene="psbA", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Sequence identification; candidate validation", inv=[43], cite="Yu et al. 2021"),
    dict(id="Eb_TrbcL", name="Eb-TrbcL", acc="NC_043882.1", gene="rbcL", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Sequence identification; candidate validation", inv=[43], cite="Yu et al. 2021"),
    dict(id="Eb_Trps16", name="Eb-Trps16", acc="NC_043882.1", gene="rps16", ftype="CDS", region="down", length=200, cls="terminator/3UTR", evidence="Sequence identification; candidate validation", inv=[43], cite="Yu et al. 2021"),
    dict(id="Eb_rbcL_SD", name="rbcL-derived SD region", acc="NC_043882.1", gene="rbcL", ftype="CDS", region="up", length=30, cls="SD", evidence="Sequence identification; candidate validation", inv=[44], cite="Yu et al. 2021"),
    dict(id="Eb_IEE", name="Eb-IEE (psbH 5' processing element)", acc="NC_043882.1", gene="psbH", ftype="CDS", region="up", length=60, cls="IEE", evidence="Sequence identification; candidate validation", inv=[45], cite="Yu et al. 2021"),
    # --- Maize: PPR-controlled UTRs ----------------------------------------
    dict(id="Zm_atpB_5UTR", name="atpB/E 5'UTR (PPR target)", acc="NC_001666.2", gene="atpB", ftype="CDS", region="up", length=200, cls="5UTR", evidence="Protein-RNA interaction; genetic analysis", inv=[46, 47], cite="Barkan and Small 2014"),
    dict(id="Zm_atpH_PPR10", name="PPR10 binding site (atpI-atpH spacer, atpH 5'UTR)", acc="NC_001666.2", gene="atpH", ftype="CDS", region="up", length=60, cls="5UTR/PPR site", evidence="Protein-RNA binding; genetic analysis; sRNA footprint", inv=[46], cite="Pfalz et al. 2009; Prikryl et al. 2011"),
    # --- heterologous phage leader -----------------------------------------
    dict(id="T7g10_leader", name="Bacteriophage T7 gene 10 5'UTR", acc=None, gene=None, ftype=None, region="seq", length=len(T7G10_LEADER), cls="5UTR", evidence="Reporter assay; transplastomic expression (heterologous)", inv=[28], cite="Kuroda and Maliga 2001"),
]

# inventory rows that are class-level / dataset-level and therefore inform motif
# models rather than sequence homology
INVENTORY_CLASS_LEVEL = {
    15: "Genome-wide TSS catalogue; represented by At_* TSS loci and by PEP/NEP motif models",
    16: "Processed 5' ends; represented by At_PrbcL_TSS and TIR-structure analysis",
    17: "Genome-wide 3' ends; represented by 3' stem-loop analysis",
    18: "PEP -10/-35 consensus; implemented as PSSM model (PEP scan)",
    19: "NEP type Ia/Ib consensus; implemented as YRTA/GAA-box composite scan",
    20: "Barley TSS catalogue; represented by Hv_* TSS loci",
    21: "Barley 5'UTR stem-loops; represented by TIR-structure analysis",
    27: "SD-dependent translation; implemented as anti-SD hybridisation model",
    29: "Chlamydomonas part collection; represented by Cr_* loci",
    32: "Marchantia TSS map; represented by Mp_* loci",
    38: "Synthetic RBS operons; informs SD model (no native plastome locus)",
    47: "Maize/rice TSS catalogues; represented by Zm_* loci",
}


def build_reference_library(refdir: Path, inventory: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genomes = {}
    rows = []
    for el in REF_ELEMENTS:
        r = dict(el)
        r["species"] = REFERENCE_GENOMES.get(el["acc"], "Escherichia phage T7") if el["acc"] else "Escherichia phage T7"
        if el["region"] == "seq":
            r.update(sequence=T7G10_LEADER, mapped=True, ref_coord="synthetic/heterologous", anchor_start_codon="")
            rows.append(r)
            continue
        if el["acc"] not in genomes:
            genomes[el["acc"]] = Plastome(refdir / f"{el['acc']}.gb")
        pg = genomes[el["acc"]]
        loc = pg.locus(el["gene"], el["ftype"])
        if loc is None:
            r.update(sequence="", mapped=False, ref_coord="anchor gene not annotated")
            rows.append(r)
            continue
        if el["region"] == "up":
            seq = pg.upstream(loc.five, loc.strand, el["length"])
            a, b = pg.genome_coord(loc.five, loc.strand, -el["length"]), pg.genome_coord(loc.five, loc.strand, -1)
        else:
            seq = pg.downstream(loc.three, loc.strand, el["length"])
            a = pg.genome_coord(loc.three, loc.strand, 1)
            b = pg.genome_coord(loc.three, loc.strand, el["length"])
        r.update(sequence=seq, mapped=True, ref_coord=f"{min(a, b)}-{max(a, b)} ({'+' if loc.strand == 1 else '-'})",
                 anchor_start_codon=loc.start_codon)
        rows.append(r)
    lib = pd.DataFrame(rows)
    lib["gc"] = lib["sequence"].map(lambda s: round(gc(s), 4) if s else np.nan)
    # inventory mapping table
    inv_rows = []
    for i, row in inventory.reset_index(drop=True).iterrows():
        ids = [e["id"] for e in REF_ELEMENTS if i in e["inv"]]
        status = "mapped to plastome locus" if ids else ("class-level model" if i in INVENTORY_CLASS_LEVEL else "not mappable")
        if i == 28:
            status = "heterologous sequence (T7 phage)"
        inv_rows.append(dict(inventory_row=i + 1, species=row.get("species"), element_name=row.get("element_name"),
                             element_class=row.get("element_class"), evidence_type=row.get("evidence_type"),
                             mapping_status=status, reference_element_ids="; ".join(ids),
                             how_used=INVENTORY_CLASS_LEVEL.get(i, "sequence homology (BLASTn + anchored global alignment)" if ids else "")))
    return lib, pd.DataFrame(inv_rows)


# ---------------------------------------------------------------------------
# PEP (sigma70-type) promoter model
# ---------------------------------------------------------------------------
def consensus_pssm(consensus: str, bg: np.ndarray, p_cons: float) -> np.ndarray:
    """Log2-odds PSSM (len x 4) from a consensus, P(cons)=p_cons, others equal."""
    m = np.zeros((len(consensus), 4))
    for i, b in enumerate(consensus):
        probs = np.full(4, (1 - p_cons) / 3)
        probs["ACGT".index(b)] = p_cons
        m[i] = np.log2(probs / bg)
    return m


class PEPModel:
    """-35 TTGACA / -10 TATAAT PSSMs with a spacer-length penalty."""

    def __init__(self, bg: np.ndarray, p=PARAMS):
        self.m35 = consensus_pssm("TTGACA", bg, p["pssm_consensus_prob"])
        self.m10 = consensus_pssm("TATAAT", bg, p["pssm_consensus_prob"])
        self.smin, self.smax, self.sopt = p["pep_spacer_min"], p["pep_spacer_max"], p["pep_spacer_opt"]
        self.pen = p["pep_spacer_penalty"]
        self.max_score = self.m35.max(1).sum() + self.m10.max(1).sum()

    @staticmethod
    def _motif_scores(X: np.ndarray, M: np.ndarray) -> np.ndarray:
        # X: (B, L) int8 with -1 for N; returns (B, L-w+1)
        B, L = X.shape
        w = M.shape[0]
        n = L - w + 1
        out = np.zeros((B, n))
        Xc = np.where(X < 0, 0, X)
        for k in range(w):
            sc = M[k][Xc[:, k:k + n]]
            sc[X[:, k:k + n] < 0] = -2.0
            out += sc
        return out

    def scan(self, X: np.ndarray):
        """Return best total score and (pos35, spacer) for each row of X."""
        B, L = X.shape
        s35 = self._motif_scores(X, self.m35)
        s10 = self._motif_scores(X, self.m10)
        best = np.full(B, -np.inf)
        bpos = np.zeros(B, dtype=int)
        bsp = np.zeros(B, dtype=int)
        for sp in range(self.smin, self.smax + 1):
            off = 6 + sp
            n = s10.shape[1] - off
            if n <= 0:
                continue
            tot = s35[:, :n] + s10[:, off:off + n] - self.pen * abs(sp - self.sopt)
            idx = tot.argmax(1)
            val = tot[np.arange(B), idx]
            upd = val > best
            best[upd], bpos[upd], bsp[upd] = val[upd], idx[upd], sp
        return best, bpos, bsp


def pep_scan_windows(windows: pd.DataFrame, model: PEPModel, rng, n_shuf: int) -> pd.DataFrame:
    """windows: DataFrame with columns id, seq (5'->3' sense, ending at start codon)."""
    out = []
    for r in windows.itertuples():
        s = r.seq
        if len(s) < 6 + model.smin + 6:
            continue
        X = encode(s)[None, :]
        sc, pos, sp = model.scan(X)
        shuf = np.stack([encode(dinuc_shuffle(s, rng)) for _ in range(n_shuf)])
        null, _, _ = model.scan(shuf)
        p = empirical_p(sc[0], null)
        p35 = int(pos[0])
        p10 = p35 + 6 + int(sp[0])
        L = len(s)
        out.append(dict(id=r.id, window_len=L, pep_score=round(float(sc[0]), 3),
                        pep_score_rel=round(float(sc[0] / model.max_score), 3),
                        minus35=s[p35:p35 + 6], minus10=s[p10:p10 + 6], spacer=int(sp[0]),
                        minus35_pos=p35 - L, minus10_pos=p10 - L,
                        predicted_tss_pos=p10 + 6 + 7 - L,  # TSS ~7 nt downstream of -10 hexamer
                        null_mean=round(float(null.mean()), 3), null_sd=round(float(null.std()), 3),
                        z=round(float((sc[0] - null.mean()) / (null.std() + 1e-9)), 3),
                        p_empirical=p))
    df = pd.DataFrame(out)
    if not df.empty:
        df["q_BH"] = bh_fdr(df["p_empirical"].values)
    return df


# ---------------------------------------------------------------------------
# NEP promoter motif enrichment (YRTA core; GAA-box + YRTA composite, type Ib)
# ---------------------------------------------------------------------------
RE_YRTA = re.compile(r"(?=([CT][AG]TA))")
RE_GAA = re.compile(r"(?=(GAA))")


def nep_counts(s: str) -> tuple[int, int]:
    y = [m.start() for m in RE_YRTA.finditer(s)]
    g = [m.start() for m in RE_GAA.finditer(s)]
    comp = 0
    gs = np.array(g)
    for j in y:
        if len(gs) and np.any(((j - (gs + 3)) >= 14) & ((j - (gs + 3)) <= 20)):
            comp += 1
    return len(y), comp


def nep_scan_windows(windows: pd.DataFrame, rng, n_shuf: int) -> pd.DataFrame:
    out = []
    for r in windows.itertuples():
        s = r.seq
        if len(s) < 30:
            continue
        ny, nc = nep_counts(s)
        null = np.array([nep_counts(dinuc_shuffle(s, rng))[1] for _ in range(n_shuf)])
        out.append(dict(id=r.id, window_len=len(s), YRTA_count=ny,
                        YRTA_per_100nt=round(100 * ny / len(s), 2), GAA_YRTA_composite=nc,
                        composite_null_mean=round(float(null.mean()), 3),
                        p_empirical=empirical_p(nc, null)))
    df = pd.DataFrame(out)
    if not df.empty:
        df["q_BH"] = bh_fdr(df["p_empirical"].values)
    return df


# ---------------------------------------------------------------------------
# Shine-Dalgarno: anti-SD hybridisation free energy
# ---------------------------------------------------------------------------
def anti_sd_from_rrn16(pg: Plastome, n: int) -> str:
    loc = pg.locus("rrn16", "rRNA")
    if loc is None:
        return "GATCACCTCCTTT"[-n:]
    seq = pg.around_start(loc.five, loc.strand, 0, loc.length)
    return seq[-n:]


def duplex_dG(a: str, b: str) -> float:
    """Minimum free energy of intermolecular hybridisation (kcal/mol), RNA."""
    if not a or not b:
        return 0.0
    d = RNA.duplexfold(a.replace("T", "U"), b.replace("T", "U"))
    return min(0.0, float(d.energy))


def sd_analysis(pg: Plastome, antisd: str, rng, n_shuf: int, species: str) -> pd.DataFrame:
    a, b = PARAMS["sd_window"]
    rows = []
    for l in pg.loci:
        if l.ftype != "CDS":
            continue
        w = pg.around_start(l.five, l.strand, a, b + 1)
        dG = duplex_dG(w, antisd)
        d = RNA.duplexfold(w.replace("T", "U"), antisd.replace("T", "U"))
        # i = end position (1-based) of the duplex on the mRNA window
        core_end = int(d.i)
        core_len = len(d.structure.split("&")[0])
        core = w[max(0, core_end - core_len):core_end]
        dist = len(w) - core_end  # nt between SD core 3' end and start codon
        null = np.array([duplex_dG(dinuc_shuffle(w, rng), antisd) for _ in range(n_shuf)])
        rows.append(dict(species=species, gene=l.gene, start_codon=l.start_codon,
                         sd_window=w, sd_core=core, sd_core_to_start_nt=dist,
                         dG_kcal_mol=round(dG, 2), null_mean_dG=round(float(null.mean()), 2),
                         p_empirical=empirical_p(dG, null, greater=False),
                         SD_present=dG <= PARAMS["sd_dG_threshold"]))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_BH"] = bh_fdr(df["p_empirical"].values)
    return df


def sd_positional_profile(pg: Plastome, antisd: str, offsets=range(-35, 11), w: int = 10) -> pd.DataFrame:
    rows = []
    for l in pg.loci:
        if l.ftype != "CDS":
            continue
        for o in offsets:
            s = pg.around_start(l.five, l.strand, o, o + w)
            rows.append(dict(gene=l.gene, offset=o, dG=duplex_dG(s, antisd)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# RNA secondary structure: TIR accessibility and 3' stem-loops
# ---------------------------------------------------------------------------
def mfe(s: str) -> tuple[str, float]:
    st, e = RNA.fold(s.replace("T", "U"))
    return st, float(e)


def _tir_one(args):
    gene, s, n_shuf, seed = args
    rng = np.random.default_rng(seed)
    st, e = mfe(s)
    null = np.array([mfe(dinuc_shuffle(s, rng))[1] for _ in range(n_shuf)])
    return dict(gene=gene, tir_seq=s, structure=st, mfe=round(e, 2),
                null_mean=round(float(null.mean()), 2), z=round(float((e - null.mean()) / (null.std() + 1e-9)), 2),
                p_more_structured=empirical_p(e, null, greater=False),
                p_less_structured=empirical_p(e, null, greater=True))


def _pmap(fn, jobs, procs):
    if procs <= 1:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(procs) as pool:
        return pool.map(fn, jobs, chunksize=1)


def tir_structure(pg: Plastome, rng, n_shuf: int, procs: int = 1) -> pd.DataFrame:
    a, b = PARAMS["tir_window"]
    jobs = [(l.gene, pg.around_start(l.five, l.strand, a, b), n_shuf, PARAMS["seed"] + 1000 + k)
            for k, l in enumerate(pg.loci) if l.ftype == "CDS"]
    df = pd.DataFrame(_pmap(_tir_one, jobs, procs))
    if not df.empty:
        df["q_more_structured"] = bh_fdr(df["p_more_structured"])
        df["q_less_structured"] = bh_fdr(df["p_less_structured"])
    return df


def best_window_mfe(s: str, w: int, step: int):
    best = (0.0, "", 0, "")
    if len(s) <= w:
        st, e = mfe(s)
        return e, st, 0, s
    for i in range(0, len(s) - w + 1, step):
        sub = s[i:i + w]
        st, e = mfe(sub)
        if e < best[0]:
            best = (e, st, i, sub)
    if best[3] == "":
        sub = s[:w]
        st, e = mfe(sub)
        best = (e, st, 0, sub)
    return best


def hairpin_stats(structure: str) -> tuple[int, int]:
    """Longest contiguous helix (bp) and number of hairpin loops."""
    stack, pairs = [], {}
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")" and stack:
            j = stack.pop()
            pairs[j] = i
    longest, cur = 0, 0
    for i in sorted(pairs):
        if i - 1 in pairs and pairs[i - 1] == pairs[i] + 1:
            cur += 1
        else:
            cur = 1
        longest = max(longest, cur)
    loops = len(re.findall(r"\(\.+\)", structure))
    return longest, loops


def _stemloop_one(args):
    gene, region, free, n_shuf, seed, gstart = args
    w, step = PARAMS["stemloop_window"], PARAMS["stemloop_step"]
    rng = np.random.default_rng(seed)
    e, st, off, sub = best_window_mfe(region, w, step)
    null = np.array([best_window_mfe(dinuc_shuffle(region, rng), w, step)[0] for _ in range(n_shuf)])
    helix, loops = hairpin_stats(st)
    down = region[off + len(sub): off + len(sub) + 10]
    return dict(gene=gene, free_3prime_spacer=free, region_len=len(region), window_offset=off + 1,
                stemloop_seq=sub, structure=st, mfe=round(e, 2), longest_helix_bp=helix, hairpin_loops=loops,
                downstream_U_fraction=round(down.count("T") / max(1, len(down)), 2),
                null_mean=round(float(null.mean()), 2), z=round(float((e - null.mean()) / (null.std() + 1e-9)), 2),
                p_empirical=empirical_p(e, null, greater=False), region_genome_start=gstart)


def stemloop_3prime(pg: Plastome, rng, n_shuf: int, procs: int = 1) -> pd.DataFrame:
    Lmax = PARAMS["stemloop_region_max"]
    jobs = []
    for k, l in enumerate(pg.loci):
        if l.ftype != "CDS":
            continue
        free = pg.free_downstream(l.three, l.strand, Lmax)
        n = max(free, 40)  # at least 40 nt analysed (may overlap the downstream gene)
        region = pg.downstream(l.three, l.strand, min(n, Lmax))
        jobs.append((l.gene, region, free, n_shuf, PARAMS["seed"] + 5000 + k, pg.genome_coord(l.three, l.strand, 1)))
    df = pd.DataFrame(_pmap(_stemloop_one, jobs, procs))
    if not df.empty:
        df["q_BH"] = bh_fdr(df["p_empirical"])
    return df


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------
def make_aligner(free_target_end_gaps: bool = False) -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2, -3
    al.open_gap_score, al.extend_gap_score = -5, -2
    if free_target_end_gaps:
        al.end_insertion_score = 0.0
        al.end_deletion_score = 0.0
    return al


ALN = make_aligner()
ALN_FREE = make_aligner(True)


def global_identity(a: str, b: str) -> tuple[float, object]:
    """Identity = identical aligned columns / len(a) (a = reference element)."""
    if not a or not b:
        return 0.0, None
    aln = ALN.align(a, b)[0]
    ident = 0
    for (ts, te), (qs, qe) in zip(*aln.aligned):
        ident += sum(1 for x, y in zip(a[ts:te], b[qs:qe]) if x == y)
    return ident / len(a), aln


def position_map(query: str, subject: str):
    """Map every query position to a subject position (or -1); free end gaps."""
    aln = ALN_FREE.align(query, subject)[0]
    qmap = np.full(len(query), -1, dtype=int)
    for (ts, te), (ss, se) in zip(*aln.aligned):
        qmap[ts:te] = np.arange(ss, se)
    return qmap


def subject_segment(subject: str, qmap: np.ndarray, qa: int, qb: int) -> str:
    idx = qmap[qa:qb]
    idx = idx[idx >= 0]
    if len(idx) == 0:
        return ""
    return subject[idx.min(): idx.max() + 1]


# ---------------------------------------------------------------------------
# BLAST wrappers
# ---------------------------------------------------------------------------
BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


def write_fasta(seqs: dict[str, str], path: Path) -> Path:
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")
    return path


def make_blastdb(fasta: Path) -> None:
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids"],
                   check=True, stdout=subprocess.DEVNULL)


def run_blastn(query: Path, db: Path, out: Path, task="blastn", evalue=10.0, word_size=7,
               threads=4, max_target=5000, extra=None) -> pd.DataFrame:
    cmd = ["blastn", "-task", task, "-query", str(query), "-db", str(db), "-evalue", str(evalue),
           "-word_size", str(word_size), "-outfmt", "6 " + BLAST_COLS, "-num_threads", str(threads),
           "-max_target_seqs", str(max_target), "-dust", "no", "-out", str(out)]
    if extra:
        cmd += extra
    subprocess.run(cmd, check=True)
    df = pd.read_csv(out, sep="\t", header=None, names=BLAST_COLS.split())
    return df


# ---------------------------------------------------------------------------
# Homology of validated elements to M. dubia (and Myrtaceae)
# ---------------------------------------------------------------------------
def homology_to_target(lib: pd.DataFrame, tg: Plastome, model: PEPModel, workdir: Path, rng,
                       threads: int) -> pd.DataFrame:
    # 1) BLASTn against the M. dubia plastome
    qfa = write_fasta({r.id: r.sequence for r in lib.itertuples() if r.mapped and r.sequence}, workdir / "ref_elements.fa")
    tfa = write_fasta({tg.accession: tg.seq}, workdir / "target.fa")
    make_blastdb(tfa)
    bl = run_blastn(qfa, tfa, workdir / "ref_vs_target.tsv", evalue=10, word_size=7, threads=threads)
    # null windows for anchored identity: intergenic-rich random windows
    nongenic = ~tg.mask
    rows = []
    for r in lib.itertuples():
        if not r.mapped or not r.sequence:
            continue
        L = len(r.sequence)
        loc = tg.locus(r.gene, r.ftype) if r.gene else None
        ortho, ortho_coord = "", ""
        if loc is not None:
            if r.region == "up":
                ortho = tg.upstream(loc.five, loc.strand, L)
                a = tg.genome_coord(loc.five, loc.strand, -L); b = tg.genome_coord(loc.five, loc.strand, -1)
            else:
                ortho = tg.downstream(loc.three, loc.strand, L)
                a = tg.genome_coord(loc.three, loc.strand, 1); b = tg.genome_coord(loc.three, loc.strand, L)
            ortho_coord = f"{min(a, b)}-{max(a, b)} ({'+' if loc.strand == 1 else '-'})"
        ident, aln = global_identity(r.sequence, ortho) if ortho else (np.nan, None)
        # empirical null: random windows of length L with >= 50 % non-genic sequence,
        # excluding +/- L around the orthologous window
        null = []
        tries = 0
        excl_center = None
        if loc is not None:
            excl_center = loc.five if r.region == "up" else loc.three
        while len(null) < PARAMS["homology_null_n"] and tries < 50 * PARAMS["homology_null_n"]:
            tries += 1
            st = int(rng.integers(0, tg.L))
            if excl_center is not None and min(abs(st - excl_center), tg.L - abs(st - excl_center)) < 2 * L:
                continue
            idx = np.arange(st, st + L) % tg.L
            if nongenic[idx].mean() < 0.5:
                continue
            w = circ_slice(tg.seq, st, st + L)
            if rng.random() < 0.5:
                w = revcomp(w)
            null.append(global_identity(r.sequence, w)[0])
        null = np.array(null)
        p = empirical_p(ident, null) if not np.isnan(ident) and len(null) else np.nan
        # best BLAST hit and synteny
        h = bl[bl.qseqid == r.id].sort_values("bitscore", ascending=False)
        best = h.iloc[0] if len(h) else None
        syntenic = False
        if best is not None and loc is not None:
            s0, s1 = sorted((int(best.sstart), int(best.send)))
            o0, o1 = [int(x) for x in ortho_coord.split()[0].split("-")]
            syntenic = not (s1 < o0 - 100 or s0 > o1 + 100) and (o1 - o0) < tg.L / 2
        # project PEP motifs (promoter classes)
        proj = {}
        if "promoter" in r.cls and aln is not None and len(r.sequence) >= 30:
            sc, pos, sp = model.scan(encode(r.sequence)[None, :])
            p35, p10 = int(pos[0]), int(pos[0]) + 6 + int(sp[0])
            qmap = np.full(L, -1)
            for (ts, te), (qs, qe) in zip(*aln.aligned):
                qmap[ts:te] = np.arange(qs, qe)
            def proj_hex(p0):
                idx = qmap[p0:p0 + 6]
                return "".join(ortho[i] if i >= 0 else "-" for i in idx)
            m35, m10 = proj_hex(p35), proj_hex(p10)
            proj = dict(ref_minus35=r.sequence[p35:p35 + 6], ref_minus10=r.sequence[p10:p10 + 6],
                        ref_spacer=int(sp[0]), ref_minus35_pos=p35 - L, ref_pep_score=round(float(sc[0]), 2),
                        Mdubia_minus35=m35, Mdubia_minus10=m10,
                        minus35_identical=m35 == r.sequence[p35:p35 + 6],
                        minus10_identical=m10 == r.sequence[p10:p10 + 6])
        rows.append(dict(element_id=r.id, element=r.name, reference_species=r.species, element_class=r.cls,
                         evidence=r.evidence, citation=r.cite, ref_length=L,
                         Mdubia_ortholog_gene=r.gene or "", Mdubia_ortholog_coord=ortho_coord,
                         anchored_identity=round(float(ident), 4) if not np.isnan(ident) else np.nan,
                         null_identity_mean=round(float(null.mean()), 4) if len(null) else np.nan,
                         null_identity_p99=round(float(np.quantile(null, 0.99)), 4) if len(null) else np.nan,
                         p_empirical=p,
                         blast_best_pident=float(best.pident) if best is not None else np.nan,
                         blast_best_alnlen=int(best.length) if best is not None else 0,
                         blast_best_qcov=round(float((best.qend - best.qstart + 1) / L), 3) if best is not None else 0,
                         blast_best_evalue=float(best.evalue) if best is not None else np.nan,
                         blast_best_sstart=int(best.sstart) if best is not None else np.nan,
                         blast_best_send=int(best.send) if best is not None else np.nan,
                         blast_syntenic=syntenic, **proj))
    df = pd.DataFrame(rows)
    df["q_BH"] = bh_fdr(df["p_empirical"].values)
    blast_ok = (df.blast_best_evalue <= PARAMS["blast_evalue"]) & df.blast_syntenic
    anch_ok = df.q_BH < PARAMS["fdr_alpha"]
    df["anchored_significant"] = anch_ok
    df["blast_syntenic_significant"] = blast_ok
    df["homology_support"] = np.select([blast_ok & anch_ok, blast_ok, anch_ok],
                                       ["strong (BLAST + anchored alignment)", "BLAST only", "anchored alignment only"], "none")
    # conservative criterion used for tiering: significant, syntenic BLASTn hit
    df["homology_supported"] = blast_ok
    return df


# ---------------------------------------------------------------------------
# Comparative genomics across Myrtaceae (annotation-independent)
# ---------------------------------------------------------------------------
UP_LEN, ANCHOR, DN_LEN = 300, 60, 150


def target_queries(tg: Plastome) -> dict[str, str]:
    q = {}
    for l in tg.loci:
        if l.ftype in ("CDS", "rRNA", "tRNA"):
            q[f"{l.gene}|UP"] = tg.around_start(l.five, l.strand, -UP_LEN, ANCHOR)
        if l.ftype == "CDS":
            down = tg.downstream(l.three, l.strand, DN_LEN)
            if l.strand == 1:
                anchor = circ_slice(tg.seq, l.three - ANCHOR + 1, l.three + 1)
            else:
                anchor = revcomp(circ_slice(tg.seq, l.three, l.three + ANCHOR))
            q[f"{l.gene}|DN"] = anchor + down
    return q


def map_queries_to_genomes(queries: dict[str, str], genomes: dict[str, Plastome], workdir: Path,
                           threads: int) -> dict[tuple[str, str], dict]:
    dbfa = write_fasta({acc: pg.seq for acc, pg in genomes.items()}, workdir / "myrtaceae_db.fa")
    make_blastdb(dbfa)
    qfa = write_fasta(queries, workdir / "target_loci.fa")
    bl = run_blastn(qfa, dbfa, workdir / "loci_vs_myrtaceae.tsv", task="blastn", evalue=1e-10,
                    word_size=11, threads=threads, max_target=len(genomes) + 5, extra=["-max_hsps", "1"])
    bl["sseqid"] = bl["sseqid"].astype(str).str.extract(r"([A-Z]{1,2}_?\d+\.\d+)", expand=False)
    bl = bl.sort_values("bitscore", ascending=False).drop_duplicates(["qseqid", "sseqid"])
    out = {}
    for h in bl.itertuples():
        pg = genomes.get(h.sseqid)
        if pg is None:
            continue
        q = queries[h.qseqid]
        pad = 20
        if h.sstart <= h.send:
            a = h.sstart - 1 - (h.qstart - 1) - pad
            b = h.send + (h.qlen - h.qend) + pad
            seg = circ_slice(pg.seq, a, b)
        else:
            a = h.send - 1 - (h.qlen - h.qend) - pad
            b = h.sstart + (h.qstart - 1) + pad
            seg = revcomp(circ_slice(pg.seq, a, b))
        qmap = position_map(q, seg)
        out[(h.qseqid, h.sseqid)] = dict(seg=seg, qmap=qmap, pident=h.pident,
                                         qcov=(h.qend - h.qstart + 1) / h.qlen, evalue=h.evalue)
    return out


def identity_vector(query: str, m: dict) -> np.ndarray:
    seg, qmap = m["seg"], m["qmap"]
    return np.array([(qmap[i] >= 0 and seg[qmap[i]] == query[i]) for i in range(len(query))], dtype=bool)


def phylogeny(tg: Plastome, genomes: dict[str, Plastome], workdir: Path, threads: int,
              min_len: int = 200) -> tuple[Path, pd.DataFrame]:
    """ML tree (FastTree GTR+G) from concatenated MAFFT alignments of CDS."""
    cds = {}
    for l in tg.loci:
        if l.ftype == "CDS" and l.length >= min_len and not l.gene.startswith("rps12") and l.gene != "ycf1":
            cds[l.gene] = tg.around_start(l.five, l.strand, 0, l.length) if l.length < 20000 else None
    cds = {k: v for k, v in cds.items() if v}
    qfa = write_fasta(cds, workdir / "cds_queries.fa")
    dbfa = workdir / "myrtaceae_db.fa"
    bl = run_blastn(qfa, dbfa, workdir / "cds_vs_myrtaceae.tsv", task="megablast", evalue=1e-20,
                    word_size=28, threads=threads, max_target=len(genomes) + 5, extra=["-max_hsps", "1"])
    bl["sseqid"] = bl["sseqid"].astype(str).str.extract(r"([A-Z]{1,2}_?\d+\.\d+)", expand=False)
    bl = bl.sort_values("bitscore", ascending=False).drop_duplicates(["qseqid", "sseqid"])
    aln_dir = workdir / "cds_alignments"
    aln_dir.mkdir(exist_ok=True)
    concat = defaultdict(str)
    part_rows = []
    accs = list(genomes.keys())
    for gname in cds:
        h = bl[bl.qseqid == gname]
        if h.sseqid.nunique() < 0.9 * len(accs):
            continue
        seqs = {}
        for r in h.itertuples():
            pg = genomes[r.sseqid]
            if r.sstart <= r.send:
                s = circ_slice(pg.seq, r.sstart - 1 - (r.qstart - 1), r.send + (r.qlen - r.qend))
            else:
                s = revcomp(circ_slice(pg.seq, r.send - 1 - (r.qlen - r.qend), r.sstart + (r.qstart - 1)))
            seqs[r.sseqid] = s
        fa = write_fasta(seqs, aln_dir / f"{gname}.fa")
        res = subprocess.run(["mafft", "--auto", "--thread", str(threads), "--quiet", str(fa)],
                             capture_output=True, text=True, check=True)
        aligned = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(__import__("io").StringIO(res.stdout), "fasta")}
        alen = len(next(iter(aligned.values())))
        for acc in accs:
            concat[acc] += aligned.get(acc, "-" * alen)
        part_rows.append(dict(gene=gname, aligned_length=alen, n_taxa=len(aligned)))
    sup = workdir / "concatenated_cds.fasta"
    write_fasta(dict(concat), sup)
    tree = workdir / "myrtaceae_cds_fasttree.nwk"
    with open(tree, "w") as fh:
        subprocess.run(["FastTree", "-nt", "-gtr", "-gamma", "-quiet", str(sup)], stdout=fh,
                       stderr=subprocess.DEVNULL, check=True)
    return tree, pd.DataFrame(part_rows)


def covariation(query_struct: str, q_off: int, query: str, m: dict) -> dict:
    """Classify base pairs of a M. dubia stem-loop in an orthologous sequence."""
    canon = {"AT", "TA", "GC", "CG", "GT", "TG"}
    stack, pairs = [], []
    for i, c in enumerate(query_struct):
        if c == "(":
            stack.append(i)
        elif c == ")" and stack:
            pairs.append((stack.pop(), i))
    res = Counter()
    seg, qmap = m["seg"], m["qmap"]
    for i, j in pairs:
        qi, qj = q_off + i, q_off + j
        si, sj = qmap[qi], qmap[qj]
        if si < 0 or sj < 0:
            res["gapped"] += 1
            continue
        x, y = seg[si], seg[sj]
        qx, qy = query[qi], query[qj]
        if x == qx and y == qy:
            res["identical"] += 1
        elif (x + y) in canon and x != qx and y != qy:
            res["compensatory"] += 1
        elif (x + y) in canon:
            res["consistent_single"] += 1
        else:
            res["disrupted"] += 1
    return dict(res)


# ---------------------------------------------------------------------------
# Window builders for de novo scans
# ---------------------------------------------------------------------------
def promoter_windows(pg: Plastome) -> pd.DataFrame:
    rows = []
    for l in pg.loci:
        free = pg.free_upstream(l.five, l.strand, PARAMS["promoter_window_max"])
        if free < PARAMS["promoter_window_min"]:
            continue
        rows.append(dict(id=l.gene, ftype=l.ftype, strand="+" if l.strand == 1 else "-",
                         seq=pg.upstream(l.five, l.strand, free), free_upstream=free,
                         window_genome_start=pg.genome_coord(l.five, l.strand, -free),
                         window_genome_end=pg.genome_coord(l.five, l.strand, -1)))
    return pd.DataFrame(rows)


def control_windows(pg: Plastome, rng, n: int = 100, L: int = 300) -> pd.DataFrame:
    """Negative controls: convergent spacers (both strands) and random in-frame CDS windows."""
    rows = []
    igs = pg.intergenic_spacers(min_len=60).drop_duplicates("seq")
    for r in igs[igs.orientation == "convergent"].itertuples():
        s = r.seq[-L:]
        rows.append(dict(id=f"conv_{r.left_gene}_{r.right_gene}_+", seq=s, control="convergent spacer"))
        rows.append(dict(id=f"conv_{r.left_gene}_{r.right_gene}_-", seq=revcomp(r.seq)[-L:], control="convergent spacer"))
    cds = [l for l in pg.loci if l.ftype == "CDS" and l.length >= L + 60 and l.length < 20000]
    for k in range(n):
        l = cds[int(rng.integers(len(cds)))]
        o = int(rng.integers(30, l.length - L - 30))
        rows.append(dict(id=f"cds_{l.gene}_{o}", seq=pg.around_start(l.five, l.strand, o, o + L), control="coding sequence"))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
def software_versions() -> dict:
    def ver(cmd):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            return (r.stdout or r.stderr).strip().splitlines()[0]
        except Exception:
            return "not found"
    import Bio, scipy
    return dict(python=platform.python_version(), biopython=Bio.__version__, numpy=np.__version__,
                pandas=pd.__version__, scipy=scipy.__version__, viennarna=RNA.__version__,
                blastn=ver(["blastn", "-version"]), mafft=ver(["mafft", "--version"]),
                fasttree=ver(["FastTree", "-expert"]) if shutil.which("FastTree") else "not found",
                platform=platform.platform(), pipeline=__version__)


def load_species_table(excel: Path) -> pd.DataFrame:
    df = pd.read_excel(excel)
    df["GenBank accession"] = df["GenBank accession"].astype(str).str.strip()
    df = df[df["GenBank accession"].str.match(r"^[A-Z]{1,2}_?\d+\.\d+$")].drop_duplicates("GenBank accession")
    return df.reset_index(drop=True)


def genus_to_tribe(genus: str) -> str:
    tribes = {"Myrciaria": "Myrteae", "Eugenia": "Myrteae", "Myrcia": "Myrteae", "Myrtus": "Myrteae",
              "Plinia": "Myrteae", "Psidium": "Myrteae", "Rhodomyrtus": "Myrteae", "Syzygium": "Syzygieae",
              "Eucalyptus": "Eucalypteae", "Corymbia": "Eucalypteae", "Angophora": "Eucalypteae",
              "Allosyncarpia": "Eucalypteae", "Stockwellia": "Eucalypteae", "Melaleuca": "Melaleuceae",
              "Baeckea": "Chamelaucieae", "Backhousia": "Backhousieae", "Lophostemon": "Lophostemoneae",
              "Xanthostemon": "Xanthostemoneae", "Heteropyxis": "Heteropyxideae (Psiloxyloideae)"}
    return tribes.get(genus, "other")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--excel", required=True)
    ap.add_argument("--inventory", required=True)
    ap.add_argument("--target", default=TARGET_ACCESSION)
    ap.add_argument("--genomes-dir", default="genomes")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--email", default=None)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--n-shuffles", type=int, default=PARAMS["n_shuffles_motif"])
    ap.add_argument("--n-shuffles-rna", type=int, default=PARAMS["n_shuffles_rna"])
    ap.add_argument("--seed", type=int, default=PARAMS["seed"])
    args = ap.parse_args(argv)
    PARAMS.update(seed=args.seed, n_shuffles_motif=args.n_shuffles, n_shuffles_rna=args.n_shuffles_rna)

    out = Path(args.outdir); tabs = out / "tables"; work = out / "work"
    for d in (out, tabs, work):
        d.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(out / "pipeline.log", mode="w"), logging.StreamHandler()])
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    gdir = Path(args.genomes_dir); refdir = gdir / "ref"

    # 1. inputs -----------------------------------------------------------
    log.info("[1] Loading species table and downloading/caching plastomes")
    sp = load_species_table(Path(args.excel))
    accs = [args.target] + [a for a in sp["GenBank accession"] if a != args.target]
    fetch_genbank(accs, gdir, args.email, args.api_key)
    fetch_genbank(list(REFERENCE_GENOMES), refdir, args.email, args.api_key)
    genomes = {a: Plastome(gdir / f"{a}.gb") for a in accs}
    tg = genomes[args.target]
    meta_rows = []
    for a, pg in genomes.items():
        row = sp[sp["GenBank accession"] == a]
        genus = pg.organism.split()[0]
        meta_rows.append(dict(accession=a, refseq=row["RefSeq accession"].iloc[0] if len(row) else "",
                              species=pg.organism, genus=genus, tribe=genus_to_tribe(genus), length_bp=pg.L,
                              GC=round(gc(pg.seq), 4), annotated_unique_CDS=sum(l.ftype == "CDS" for l in pg.loci),
                              annotated_gene_count_ncbi=int(row["Gene count"].iloc[0]) if len(row) else np.nan,
                              role="target" if a == args.target else "Myrtaceae comparison",
                              sha256=sha256(gdir / f"{a}.gb")))
    for a, s in REFERENCE_GENOMES.items():
        pg = Plastome(refdir / f"{a}.gb")
        meta_rows.append(dict(accession=a, refseq=a, species=pg.organism, genus=pg.organism.split()[0], tribe="reference",
                              length_bp=pg.L, GC=round(gc(pg.seq), 4),
                              annotated_unique_CDS=sum(l.ftype == "CDS" for l in pg.loci),
                              role="reference (validated elements)", sha256=sha256(refdir / f"{a}.gb")))
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(tabs / "S1_genomes.tsv", sep="\t", index=False)
    log.info("    %d Myrtaceae plastomes (incl. target) + %d reference plastomes", len(genomes), len(REFERENCE_GENOMES))

    # 2. target genome features --------------------------------------------
    log.info("[2] Target plastome features")
    igs = tg.intergenic_spacers(min_len=1)
    igs["IR_duplicate"] = igs.duplicated("seq", keep="first") & (igs.length > 20)
    igs["GC"] = igs["seq"].map(gc)
    igs.drop(columns=["seq"]).to_csv(tabs / "S2_target_spacers.tsv", sep="\t", index=False)
    loci_df = pd.DataFrame([dict(gene=l.gene, type=l.ftype, strand="+" if l.strand == 1 else "-", start_codon=l.start_codon,
                                 five_prime_coord=l.five + 1, three_prime_coord=l.three + 1, length=l.length, copies=l.copies,
                                 free_upstream_nt=tg.free_upstream(l.five, l.strand, 2000),
                                 free_downstream_nt=tg.free_downstream(l.three, l.strand, 2000)) for l in tg.loci])
    loci_df.to_csv(tabs / "S3_target_loci.tsv", sep="\t", index=False)
    bg = np.array([igs.loc[~igs.IR_duplicate, "seq"].str.count(b).sum() for b in "ACGT"], dtype=float)
    bg /= bg.sum()
    model = PEPModel(bg)
    antisd = anti_sd_from_rrn16(tg, PARAMS["antisd_length"])
    log.info("    background (IGS) ACGT = %s; anti-SD = %s", np.round(bg, 3), antisd)

    # 3. reference library and homology -----------------------------------
    log.info("[3] Reference library of validated elements and homology to the target")
    inv = pd.read_excel(args.inventory)
    lib, invmap = build_reference_library(refdir, inv)
    lib.drop(columns=["inv"]).to_csv(tabs / "S4_reference_library.tsv", sep="\t", index=False)
    invmap.to_csv(tabs / "S5_inventory_mapping.tsv", sep="\t", index=False)
    hom = homology_to_target(lib, tg, model, work, rng, args.threads)
    hom.to_csv(tabs / "S6_homology_validated_elements.tsv", sep="\t", index=False)

    # 4. de novo predictions in target --------------------------------------
    log.info("[4] De novo predictions (PEP, NEP, SD, TIR, 3' stem-loops) with shuffle null models")
    pw = promoter_windows(tg)
    pep = pep_scan_windows(pw, model, rng, args.n_shuffles).merge(
        pw[["id", "ftype", "strand", "free_upstream", "window_genome_start", "window_genome_end"]], on="id")
    pep.to_csv(tabs / "S7_PEP_promoters.tsv", sep="\t", index=False)
    ctrl = control_windows(tg, rng)
    pep_ctrl = pep_scan_windows(ctrl, model, rng, args.n_shuffles).merge(ctrl[["id", "control"]], on="id")
    pep_ctrl.to_csv(tabs / "S8_PEP_negative_controls.tsv", sep="\t", index=False)
    nep = nep_scan_windows(pw, rng, max(200, args.n_shuffles // 5))
    nep_ctrl = nep_scan_windows(ctrl, rng, max(200, args.n_shuffles // 5)).merge(ctrl[["id", "control"]], on="id")
    nep.to_csv(tabs / "S9_NEP_motifs.tsv", sep="\t", index=False)
    nep_ctrl.to_csv(tabs / "S9b_NEP_controls.tsv", sep="\t", index=False)
    sd = sd_analysis(tg, antisd, rng, args.n_shuffles_rna, tg.organism)
    sd.to_csv(tabs / "S10_SD_target.tsv", sep="\t", index=False)
    sdprof = sd_positional_profile(tg, antisd)
    sdprof.to_csv(work / "sd_profile_target.tsv", sep="\t", index=False)
    tir = tir_structure(tg, rng, args.n_shuffles_rna, args.threads)
    tir.to_csv(tabs / "S11_TIR_structure.tsv", sep="\t", index=False)
    sl = stemloop_3prime(tg, rng, PARAMS["n_shuffles_stemloop"], args.threads)
    sl.to_csv(tabs / "S12_3prime_stemloops.tsv", sep="\t", index=False)

    # 5. positive control: tobacco (validated PEP promoters) + model species SD
    log.info("[5] Positive controls on reference plastomes")
    nt = Plastome(refdir / "NC_001879.2.gb")
    nt_igs = nt.intergenic_spacers(min_len=1)
    bg_nt = np.array([nt_igs["seq"].str.count(b).sum() for b in "ACGT"], dtype=float); bg_nt /= bg_nt.sum()
    pep_nt = pep_scan_windows(promoter_windows(nt), PEPModel(bg_nt), rng, args.n_shuffles)
    pep_nt.to_csv(tabs / "S13_PEP_tobacco_positive_control.tsv", sep="\t", index=False)
    sd_ref = []
    for a in ("NC_001879.2", "NC_000932.1", "NC_008590.1", "NC_005353.1"):
        pg = Plastome(refdir / f"{a}.gb")
        sd_ref.append(sd_analysis(pg, anti_sd_from_rrn16(pg, PARAMS["antisd_length"]), rng, 100, pg.organism))
    sd_ref = pd.concat(sd_ref)
    sd_ref.to_csv(tabs / "S14_SD_reference_species.tsv", sep="\t", index=False)

    # 6. comparative genomics across Myrtaceae -------------------------------
    log.info("[6] Comparative mapping of target loci across %d Myrtaceae plastomes", len(genomes))
    queries = target_queries(tg)
    maps = map_queries_to_genomes(queries, genomes, work, args.threads)
    accs_all = list(genomes)
    # 6a region identity matrices
    ident_rows = []
    pos_cons = {}
    for qid, q in queries.items():
        gene, kind = qid.split("|")
        core = slice(0, UP_LEN) if kind == "UP" else slice(ANCHOR, ANCHOR + DN_LEN)
        mat = []
        for a in accs_all:
            m = maps.get((qid, a))
            if m is None:
                ident_rows.append(dict(gene=gene, region=kind, accession=a, identity=np.nan))
                continue
            v = identity_vector(q, m)
            mat.append(v)
            ident_rows.append(dict(gene=gene, region=kind, accession=a, identity=round(float(v[core].mean()), 4)))
        if mat:
            pos_cons[qid] = np.mean(np.stack(mat), axis=0)
    ident = pd.DataFrame(ident_rows)
    ident.to_csv(tabs / "S15_region_identity_long.tsv", sep="\t", index=False)
    np.savez_compressed(work / "position_conservation.npz", **{k.replace("|", "__"): v for k, v in pos_cons.items()})

    # 6b PEP motif conservation and re-scoring across species
    pep_cmp, mf_rows = [], []
    for r in pep.itertuples():
        qid = f"{r.id}|UP"
        if qid not in queries:
            continue
        q = queries[qid]
        p35 = UP_LEN + r.minus35_pos
        p10 = UP_LEN + r.minus10_pos
        wstart = UP_LEN - r.window_len
        for a in accs_all:
            m = maps.get((qid, a))
            if m is None:
                continue
            seg = subject_segment(m["seg"], m["qmap"], wstart, UP_LEN)
            sc = model.scan(encode(seg)[None, :])[0][0] if len(seg) > 30 else np.nan
            v = identity_vector(q, m)
            pep_cmp.append(dict(gene=r.id, accession=a, pep_score=round(float(sc), 3),
                                minus35_identical=bool(v[p35:p35 + 6].all()), minus10_identical=bool(v[p10:p10 + 6].all())))
        if qid in pos_cons:
            pc = pos_cons[qid]
            motif_idx = list(range(p35, p35 + 6)) + list(range(p10, p10 + 6))
            flank_idx = [i for i in range(max(wstart, p35 - 40), min(UP_LEN, p10 + 46)) if i not in motif_idx]
            mf_rows.append(dict(gene=r.id, q_BH=r.q_BH, motif_conservation=float(pc[motif_idx].mean()),
                                flank_conservation=float(pc[flank_idx].mean())))
    pep_cmp = pd.DataFrame(pep_cmp)
    pep_cmp.to_csv(tabs / "S16_PEP_conservation_long.tsv", sep="\t", index=False)
    mf = pd.DataFrame(mf_rows)
    mf.to_csv(tabs / "S17_motif_vs_flank_conservation.tsv", sep="\t", index=False)

    # 6c SD conservation
    sd_cmp = []
    for r in sd.itertuples():
        qid = f"{r.gene}|UP"
        if qid not in queries:
            continue
        for a in accs_all:
            m = maps.get((qid, a))
            if m is None:
                continue
            seg = subject_segment(m["seg"], m["qmap"], UP_LEN - 20, UP_LEN)
            sd_cmp.append(dict(gene=r.gene, accession=a, dG=round(duplex_dG(seg, antisd), 2)))
    sd_cmp = pd.DataFrame(sd_cmp)
    sd_cmp["SD_present"] = sd_cmp.dG <= PARAMS["sd_dG_threshold"]
    sd_cmp.to_csv(tabs / "S18_SD_conservation_long.tsv", sep="\t", index=False)

    # 6d stem-loop conservation + covariation
    sl_cmp, cov_rows = [], []
    for r in sl.itertuples():
        qid = f"{r.gene}|DN"
        if qid not in queries:
            continue
        q = queries[qid]
        qa = ANCHOR + r.window_offset - 1
        qb = qa + len(r.stemloop_seq)
        tot = Counter()
        for a in accs_all:
            if a == args.target:
                continue
            m = maps.get((qid, a))
            if m is None:
                continue
            seg = subject_segment(m["seg"], m["qmap"], qa, qb)
            e = mfe(seg)[1] if len(seg) >= 20 else np.nan
            c = covariation(r.structure, qa, q, m)
            tot.update(c)
            sl_cmp.append(dict(gene=r.gene, accession=a, mfe=round(e, 2), **c))
        npairs = sum(tot.values())
        cov_rows.append(dict(gene=r.gene, mfe_target=r.mfe, z_target=r.z, q_BH=r.q_BH,
                             pairs_evaluated=npairs, **{k: tot.get(k, 0) for k in
                             ("identical", "compensatory", "consistent_single", "disrupted", "gapped")},
                             fraction_pairing_maintained=round((tot.get("identical", 0) + tot.get("compensatory", 0) +
                                                               tot.get("consistent_single", 0)) / max(1, npairs - tot.get("gapped", 0)), 4)))
    sl_cmp = pd.DataFrame(sl_cmp)
    sl_cmp.to_csv(tabs / "S19_stemloop_conservation_long.tsv", sep="\t", index=False)
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(tabs / "S20_stemloop_covariation.tsv", sep="\t", index=False)

    # 6e homology of validated elements across all Myrtaceae (BLASTn)
    lib_fa = work / "ref_elements.fa"
    blm = run_blastn(lib_fa, work / "myrtaceae_db.fa", work / "ref_vs_myrtaceae.tsv", evalue=PARAMS["blast_evalue"],
                     word_size=7, threads=args.threads, max_target=len(genomes) + 5, extra=["-max_hsps", "1"])
    blm["sseqid"] = blm["sseqid"].astype(str).str.extract(r"([A-Z]{1,2}_?\d+\.\d+)", expand=False)
    blm = blm.sort_values("bitscore", ascending=False).drop_duplicates(["qseqid", "sseqid"])
    blm["qcov"] = (blm.qend - blm.qstart + 1) / blm.qlen
    blm.to_csv(tabs / "S21_validated_elements_vs_myrtaceae.tsv", sep="\t", index=False)

    # 7. phylogeny -----------------------------------------------------------
    log.info("[7] Phylogeny from concatenated CDS (MAFFT + FastTree)")
    tree, parts = phylogeny(tg, genomes, work, args.threads)
    parts.to_csv(tabs / "S22_phylogeny_partitions.tsv", sep="\t", index=False)
    shutil.copy(tree, out / "myrtaceae_cds_fasttree.nwk")

    # 8. evidence integration --------------------------------------------------
    log.info("[8] Evidence integration")
    integ = integrate(tg, hom, pep, nep, sd, tir, sl, cov, pep_cmp, sd_cmp, ident, args.target)
    integ.to_csv(tabs / "S23_integrated_evidence_target.tsv", sep="\t", index=False)
    summ = species_summary(meta, ident, pep_cmp, pep, sd_cmp, sl_cmp, blm)
    summ.to_csv(tabs / "S24_species_summary.tsv", sep="\t", index=False)

    # 9. statistics ------------------------------------------------------------
    st = {}
    if len(mf):
        sig = mf[mf.q_BH < PARAMS["fdr_alpha"]]
        for name, d in (("all_windows", mf), ("significant", sig)):
            if len(d) >= 5:
                w = stats.wilcoxon(d.motif_conservation, d.flank_conservation, alternative="greater")
                st[f"motif_vs_flank_{name}"] = dict(n=len(d), median_motif=float(d.motif_conservation.median()),
                                                    median_flank=float(d.flank_conservation.median()),
                                                    W=float(w.statistic), p=float(w.pvalue))
    if len(pep_ctrl):
        st["PEP_fraction_q005"] = dict(target_5prime=float((pep.q_BH < 0.05).mean()),
                                       controls={c: float((g.p_empirical < 0.05).mean()) for c, g in pep_ctrl.groupby("control")})
        u = stats.mannwhitneyu(pep.z, pep_ctrl.z, alternative="greater")
        st["PEP_z_target_vs_controls"] = dict(U=float(u.statistic), p=float(u.pvalue))
    if len(nep):
        u = stats.mannwhitneyu(nep.YRTA_per_100nt, nep_ctrl.YRTA_per_100nt, alternative="greater")
        st["YRTA_density_target_vs_controls"] = dict(median_target=float(nep.YRTA_per_100nt.median()),
                                                     median_controls=float(nep_ctrl.YRTA_per_100nt.median()),
                                                     U=float(u.statistic), p=float(u.pvalue))
    st["SD"] = dict(fraction_present=float(sd.SD_present.mean()), n=int(len(sd)),
                    fraction_q005=float((sd.q_BH < 0.05).mean()))
    for org, g in sd_ref.groupby("species"):
        st["SD"][f"fraction_present_{org}"] = float(g.SD_present.mean())
    st["TIR"] = dict(median_z=float(tir.z.median()), wilcoxon_p_z_gt_0=float(stats.wilcoxon(tir.z, alternative="greater").pvalue),
                     n_less_structured_q005=int((tir.q_less_structured < 0.05).sum()))
    st["stemloops"] = dict(n=int(len(sl)), n_q005=int((sl.q_BH < 0.05).sum()), median_z=float(sl.z.median()))
    (out / "statistics.json").write_text(json.dumps(st, indent=2, default=float))

    meta_run = dict(pipeline_version=__version__, target=args.target, target_species=tg.organism,
                    n_myrtaceae=len(genomes), parameters=PARAMS, software=software_versions(),
                    runtime_min=round((time.time() - t0) / 60, 2), anti_SD=antisd,
                    background_ACGT=[round(x, 4) for x in bg],
                    note="All predictions are in silico hypotheses requiring experimental validation.")
    (out / "run_metadata.json").write_text(json.dumps(meta_run, indent=2, default=str))
    write_excel(out, tabs)
    log.info("Done in %.1f min", (time.time() - t0) / 60)


def integrate(tg, hom, pep, nep, sd, tir, sl, cov, pep_cmp, sd_cmp, ident, target) -> pd.DataFrame:
    a = PARAMS["fdr_alpha"]
    c1 = PARAMS["conservation_tier1"]
    others = ident[ident.accession != target]
    up_id = others[others.region == "UP"].groupby("gene").identity.mean()
    dn_id = others[others.region == "DN"].groupby("gene").identity.mean()
    pc = pep_cmp[pep_cmp.accession != target].copy()
    pc["both"] = pc.minus35_identical & pc.minus10_identical
    pep_cons = pc.groupby("gene").both.mean()
    sdc = sd_cmp[sd_cmp.accession != target].groupby("gene").SD_present.mean()
    rows = []
    for l in tg.loci:
        g = l.gene
        r = dict(gene=g, type=l.ftype, strand="+" if l.strand == 1 else "-", start_codon=l.start_codon,
                 upstream_identity_myrtaceae=round(float(up_id.get(g, np.nan)), 4),
                 downstream_identity_myrtaceae=round(float(dn_id.get(g, np.nan)), 4))
        h5 = hom[(hom.Mdubia_ortholog_gene == g) & hom.homology_supported & ~hom.element_class.str.contains("terminator")]
        h3 = hom[(hom.Mdubia_ortholog_gene == g) & hom.homology_supported & hom.element_class.str.contains("terminator")]
        r["validated_5prime_homologs"] = "; ".join(f"{x.element_id} ({x.anchored_identity:.2f})" for x in h5.itertuples())
        r["validated_3prime_homologs"] = "; ".join(f"{x.element_id} ({x.anchored_identity:.2f})" for x in h3.itertuples())
        p = pep[pep.id == g]
        if len(p):
            p = p.iloc[0]
            cons = float(pep_cons.get(g, np.nan))
            r.update(PEP_minus35=p.minus35, PEP_minus10=p.minus10, PEP_spacer=p.spacer, PEP_minus35_pos=p.minus35_pos,
                     PEP_score=p.pep_score, PEP_q=p.q_BH, PEP_motif_conservation=round(cons, 4))
            hp = len(h5[h5.element_class.str.contains("promoter")]) > 0
            if hp or (p.q_BH < a and cons >= c1):
                r["promoter_tier"] = 1
            elif p.q_BH < a:
                r["promoter_tier"] = 2
            elif p.p_empirical < a:
                r["promoter_tier"] = 3
        n = nep[nep.id == g]
        if len(n):
            n = n.iloc[0]
            r.update(NEP_YRTA_per_100nt=n.YRTA_per_100nt, NEP_GAA_YRTA=n.GAA_YRTA_composite, NEP_q=n.q_BH)
        s = sd[sd.gene == g]
        if len(s):
            s = s.iloc[0]
            cons = float(sdc.get(g, np.nan))
            r.update(SD_core=s.sd_core, SD_dG=s.dG_kcal_mol, SD_dist=s.sd_core_to_start_nt, SD_present=bool(s.SD_present),
                     SD_q=s.q_BH, SD_conservation=round(cons, 4))
            # SD is called by the hybridisation-energy threshold (Ma et al. 2002);
            # the shuffle-based P (reported in Table S10) is conservative because
            # dinucleotide shuffling largely preserves purine runs, so SD tiers use
            # dG plus cross-species conservation rather than the empirical P.
            if s.SD_present and cons >= c1:
                r["SD_tier"] = 1
            elif s.SD_present and cons >= 0.5:
                r["SD_tier"] = 2
            elif s.SD_present:
                r["SD_tier"] = 3
        t = tir[tir.gene == g]
        if len(t):
            r.update(TIR_mfe=t.iloc[0].mfe, TIR_z=t.iloc[0].z)
        x = sl[sl.gene == g]
        if len(x):
            x = x.iloc[0]
            cv = cov[cov.gene == g]
            pm = float(cv.fraction_pairing_maintained.iloc[0]) if len(cv) else np.nan
            r.update(stemloop_mfe=x.mfe, stemloop_z=x.z, stemloop_q=x.q_BH, stemloop_structure=x.structure,
                     stemloop_pairing_maintained=pm, stemloop_compensatory=int(cv.compensatory.iloc[0]) if len(cv) else 0)
            if (x.q_BH < a and pm >= c1) or (len(h3) and x.q_BH < a):
                r["stemloop_tier"] = 1
            elif x.q_BH < a:
                r["stemloop_tier"] = 2
            elif x.p_empirical < a:
                r["stemloop_tier"] = 3
        rows.append(r)
    return pd.DataFrame(rows)


def species_summary(meta, ident, pep_cmp, pep, sd_cmp, sl_cmp, blm) -> pd.DataFrame:
    a = PARAMS["fdr_alpha"]
    sig = set(pep[pep.q_BH < a].id)
    pc = pep_cmp[pep_cmp.gene.isin(sig)].copy()
    pc["both"] = pc.minus35_identical & pc.minus10_identical
    rows = []
    for m in meta[meta.role != "reference (validated elements)"].itertuples():
        acc = m.accession
        i = ident[ident.accession == acc]
        rows.append(dict(accession=acc, species=m.species, genus=m.genus, tribe=m.tribe, length_bp=m.length_bp, GC=m.GC,
                         mean_upstream_identity=round(float(i[i.region == "UP"].identity.mean()), 4),
                         mean_downstream_identity=round(float(i[i.region == "DN"].identity.mean()), 4),
                         sig_PEP_motifs_conserved=round(float(pc[pc.accession == acc].both.mean()), 4) if len(pc) else np.nan,
                         SD_fraction=round(float(sd_cmp[sd_cmp.accession == acc].SD_present.mean()), 4),
                         stemloop_mean_mfe=round(float(sl_cmp[sl_cmp.accession == acc].mfe.mean()), 2) if len(sl_cmp[sl_cmp.accession == acc]) else np.nan,
                         validated_elements_hit_qcov50=int(((blm.sseqid == acc) & (blm.qcov >= 0.5)).sum())))
    return pd.DataFrame(rows)


TABLE_DESCRIPTIONS = {
    "S1_genomes": "Plastomes analysed: target, Myrtaceae comparison set and reference species (with SHA-256 of cached GenBank files).",
    "S2_target_spacers": "Intergenic spacers of the M. dubia plastome with orientation class (tandem, divergent, convergent) and GC content.",
    "S3_target_loci": "Unique genes of M. dubia (IR copies collapsed) with coordinates, start codons and free upstream/downstream spacer.",
    "S4_reference_library": "Reference library of experimentally validated (or candidate) plastid regulatory elements extracted from reference plastomes.",
    "S5_inventory_mapping": "Mapping of each of the 48 inventory records to reference loci or class-level models.",
    "S6_homology_validated_elements": "Homology of reference elements to M. dubia: anchored global identity vs empirical null, BLASTn, synteny and projected PEP motifs.",
    "S7_PEP_promoters": "De novo PEP (sigma70-type) promoter scan of M. dubia 5' spacers with dinucleotide-shuffle null and BH-FDR.",
    "S8_PEP_negative_controls": "PEP scan of negative-control windows (convergent spacers, coding sequence) for calibration.",
    "S9_NEP_motifs": "NEP YRTA density and GAA-box/YRTA composite (type Ib) motif enrichment in 5' spacers.",
    "S9b_NEP_controls": "NEP motif counts in negative-control windows.",
    "S10_SD_target": "Shine-Dalgarno analysis: anti-SD hybridisation free energy for every M. dubia CDS.",
    "S11_TIR_structure": "Translation-initiation-region (-60..+30) folding energy and z-scores vs shuffled sequence.",
    "S12_3prime_stemloops": "3' stem-loop candidates (best 60-nt window) with null-model z-scores and BH-FDR.",
    "S13_PEP_tobacco_positive_control": "Positive control: PEP scan of Nicotiana tabacum 5' spacers (validated PEP promoters psbA, rbcL, rrn16, atpB).",
    "S14_SD_reference_species": "SD analysis for reference species (tobacco, Arabidopsis, barley, Chlamydomonas).",
    "S15_region_identity_long": "Per-species identity of M. dubia 5' (UP, -300..-1) and 3' (DN, +1..+150) regions.",
    "S16_PEP_conservation_long": "Re-scored PEP promoters and -35/-10 identity in each Myrtaceae plastome.",
    "S17_motif_vs_flank_conservation": "Positional conservation of PEP -35/-10 hexamers versus flanking spacer positions.",
    "S18_SD_conservation_long": "Anti-SD hybridisation energy of orthologous SD windows in each Myrtaceae plastome.",
    "S19_stemloop_conservation_long": "Folding energy and base-pair classification of orthologous 3' stem-loops.",
    "S20_stemloop_covariation": "Summary of identical, compensatory, consistent and disrupted base pairs per M. dubia stem-loop.",
    "S21_validated_elements_vs_myrtaceae": "BLASTn hits of reference elements in every Myrtaceae plastome.",
    "S22_phylogeny_partitions": "Genes and alignment lengths used for the concatenated CDS phylogeny.",
    "S23_integrated_evidence_target": "Integrated per-gene evidence and confidence tiers for M. dubia regulatory elements.",
    "S24_species_summary": "Species-level summary of regulatory-element conservation across Myrtaceae.",
}


def write_excel(out: Path, tabs: Path, name: str = "Online_Resource_2_Supplementary_Tables.xlsx") -> Path:
    from openpyxl.styles import Font, PatternFill, Alignment
    path = out / name
    files = sorted(tabs.glob("S*.tsv"), key=lambda p: (int(re.match(r"S(\d+)", p.stem).group(1)), p.stem))
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        readme = pd.DataFrame([dict(sheet=p.stem[:31], description=TABLE_DESCRIPTIONS.get(p.stem, "")) for p in files])
        title = pd.DataFrame({"README": [
            "In silico identification of plastid regulatory elements in Myrciaria dubia (PX251231.1) and comparative analysis across Myrtaceae.",
            f"Generated by cpreg_pipeline.py v{__version__}. All predictions are in silico hypotheses requiring experimental validation.",
            "Coordinates are 1-based (GenBank). Positions relative to the start codon: first nt of the start codon = 0; -1 = first nt upstream.",
            "p_empirical = (k+1)/(N+1) from dinucleotide-preserving shuffles; q_BH = Benjamini-Hochberg FDR."]})
        title.to_excel(xw, sheet_name="README", index=False)
        readme.to_excel(xw, sheet_name="README", index=False, startrow=6)
        for p in files:
            df = pd.read_csv(p, sep="\t")
            df.to_excel(xw, sheet_name=p.stem[:31], index=False)
        for ws in xw.book.worksheets:
            for c in ws[1]:
                c.font = Font(bold=True, color="FFFFFF")
                c.fill = PatternFill("solid", fgColor="1F4E79")
                c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in list(col)[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)
    return path


if __name__ == "__main__":
    main()
