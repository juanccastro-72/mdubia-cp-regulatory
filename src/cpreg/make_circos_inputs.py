#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate circos_work/{repeats,gc}.tsv for the detailed Circos figures.

Self-repeats (direct/inverted, including the inverted repeat) are found by BLASTn
of the target plastome against itself; GC content is computed in 500-bp windows.
"""
from __future__ import annotations
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", default="genomes")
    ap.add_argument("--target", default="PX251231.1")
    ap.add_argument("--out", default="circos_work")
    ap.add_argument("--window", type=int, default=500)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rec = SeqIO.read(str(Path(a.genomes_dir) / f"{a.target}.gb"), "genbank")
    seq = str(rec.seq); L = len(seq)
    fa = out / "target.fa"; fa.write_text(f">{a.target}\n{seq}\n")
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl"], check=True, stdout=subprocess.DEVNULL)
    tsv = out / "self.tsv"
    subprocess.run(["blastn", "-query", str(fa), "-db", str(fa), "-word_size", "20", "-evalue", "1e-10",
                    "-outfmt", "6 qstart qend sstart send length pident sstrand", "-dust", "no", "-out", str(tsv)],
                   check=True)
    d = pd.read_csv(tsv, sep="\t", header=None,
                    names=["qs", "qe", "ss", "se", "length", "pid", "strand"])
    d = d[d.length <= L * 0.95]
    d = d[(d.length >= 100) & (d.pid >= 90)].copy()
    d["a"] = d[["qs", "ss"]].min(axis=1)
    d["b"] = d[["qe", "se"]].max(axis=1)
    d["kind"] = np.where(d.strand == "plus", "direct", "inverted")
    d = d.drop_duplicates(subset=["a", "b", "kind"])
    d.to_csv(out / "repeats.tsv", sep="\t", index=False)
    W = a.window
    gc = [(i + 1, (seq[i:i + W].count("G") + seq[i:i + W].count("C")) /
           max(1, sum(seq[i:i + W].count(x) for x in "ACGT"))) for i in range(0, L, W)]
    pd.DataFrame(gc, columns=["start", "gc"]).to_csv(out / "gc.tsv", sep="\t", index=False)
    print(f"wrote {out}/repeats.tsv ({len(d)} repeats) and {out}/gc.tsv ({len(gc)} windows)")


if __name__ == "__main__":
    main()
