#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpreg_synteny.py - whole-plastome synteny analysis of Myrciaria dubia (PX251231.1)
against the 121 other Myrtaceae plastomes and 12 reference plastomes.

Alignments are computed with nucmer (MUMmer); collinear and inverted blocks are
parsed from show-coords. Outputs:

  * Circos (pyCirclize) synteny plots of M. dubia vs one representative plastome
    per Myrtaceae tribe and vs the main reference species (tobacco)       -> Fig8
  * a family-wide synteny-conservation Circos: the M. dubia plastome with outer
    rings giving, per 500-bp window, the fraction of the 121 Myrtaceae that are
    collinear / inverted / unaligned, plus per-tribe collinearity rings  -> Fig9
  * a tribe-level synteny summary (blocks, inversions, coverage)          -> FigS6
  * machine-readable tables (synteny blocks; per-window conservation;
    per-genome summary).

Fonts: Times New Roman (rendered with the metric-compatible Liberation Serif),
10 pt base; scientific names and gene names in italics. SVG text is live/editable.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from Bio import SeqIO, Phylo

from cpreg_style import C, apply_style, save, TRIBE_ORDER, TRIBE_COL, tribe5, it

apply_style()

MM = 1 / 25.4
MIN_BLOCK = 300        # nt, minimum aligned block length kept
WINDOW = 500           # nt, window for family-wide conservation


# ---------------------------------------------------------------------------
# nucmer synteny
# ---------------------------------------------------------------------------
def write_fasta(acc: str, seq: str, path: Path):
    path.write_text(f">{acc}\n{seq}\n")


def nucmer_blocks(ref_fa: Path, qry_fa: Path, workdir: Path, tag: str) -> pd.DataFrame:
    pref = str(workdir / tag)
    subprocess.run(["nucmer", "--maxmatch", "-c", "100", "-p", pref, str(ref_fa), str(qry_fa)],
                   check=True, stderr=subprocess.DEVNULL)
    subprocess.run(["delta-filter", "-1", pref + ".delta"], check=True,
                   stdout=open(pref + ".1delta", "w"), stderr=subprocess.DEVNULL)
    r = subprocess.run(["show-coords", "-rclTH", pref + ".1delta"], capture_output=True, text=True, check=True)
    rows = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
        length = int(f[4]); ident = float(f[6])
        strand = 1 if (e2 >= s2) else -1
        if length < MIN_BLOCK:
            continue
        rows.append(dict(ref_start=min(s1, e1), ref_end=max(s1, e1),
                         qry_start=min(s2, e2), qry_end=max(s2, e2),
                         length=length, identity=ident, strand=strand))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Circos pairwise synteny (pyCirclize)
# ---------------------------------------------------------------------------
def circos_pair(ax_none, ref_acc, ref_len, qry_acc, qry_len, ref_species, qry_species,
                blocks: pd.DataFrame, regions, outpath: Path, name: str):
    from pycirclize import Circos
    sectors = {ref_acc: ref_len, qry_acc: qry_len}
    circos = Circos(sectors, space=5, start=5, end=355)
    col = {ref_acc: C["blue"], qry_acc: C["aqua"]}
    disp = {ref_acc: ref_species, qry_acc: qry_species}
    for sector in circos.sectors:
        tr = sector.add_track((94, 100))
        tr.axis(fc=col[sector.name], ec="none")
        # ticks every 20 kb
        sector.text(f"$\\it{{{disp[sector.name].replace(' ', chr(92)+' ')}}}$", r=112, size=10)
        tr.xticks_by_interval(20000, label_formatter=lambda v: f"{v/1000:.0f}", label_size=7,
                              tick_length=1.5)
    # IR/LSC/SSC ring on reference
    rsec = circos.get_sector(ref_acc)
    rc = {"LSC": C["blue"], "SSC": C["aqua"], "IRa": C["orange"], "IRb": C["orange"]}
    rt = rsec.add_track((88, 93))
    for s, e, lab in regions:
        rt.rect(s, e, fc=rc.get(lab, C["muted"]), ec="white", lw=0.3)
    # ribbons
    for b in blocks.itertuples():
        color = C["blue"] if b.strand == 1 else C["red"]
        circos.link((ref_acc, b.ref_start, b.ref_end), (qry_acc, b.qry_start, b.qry_end),
                    color=color, alpha=0.45 if b.strand == 1 else 0.6, r1=88, r2=88)
    fig = circos.plotfig(dpi=100)
    fig.set_size_inches(90 * MM, 90 * MM)
    return fig


def fig_circos_tribes(out, gdir, refdir, res, reps, work):
    """Fig8: multi-panel Circos, M. dubia vs one representative per tribe + tobacco."""
    tg = SeqIO.read(str(gdir / "PX251231.1.gb"), "genbank")
    ref_len = len(tg); ref_seq = str(tg.seq)
    regions = plastome_regions(tg)
    write_fasta("Mdubia", ref_seq, work / "ref.fa")
    from pycirclize import Circos
    panels = reps
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    fig = plt.figure(figsize=(180 * MM, nrow * 82 * MM))
    blocks_all = []
    for i, (tribe, acc, species, folder) in enumerate(panels):
        rec = SeqIO.read(str(folder / f"{acc}.gb"), "genbank")
        write_fasta(acc, str(rec.seq), work / "qry.fa")
        blocks = nucmer_blocks(work / "ref.fa", work / "qry.fa", work, f"p{i}")
        blocks["tribe"] = tribe; blocks["accession"] = acc; blocks_all.append(blocks)
        inv = int((blocks.strand == -1).sum())
        cov = 100 * blocks[blocks.strand == 1].eval("ref_end-ref_start").clip(lower=0).sum() / ref_len
        sectors = {"Mdubia": ref_len, acc: len(rec)}
        circ = Circos(sectors, space=6, start=8, end=352)
        for sector in circ.sectors:
            isref = sector.name == "Mdubia"
            tr = sector.add_track((92, 99))
            tr.axis(fc=C["blue"] if isref else C["aqua"], ec="none")
            tr.xticks_by_interval(40000, label_formatter=lambda v: f"{v/1000:.0f}", label_size=6, tick_length=1.2)
        for b in blocks.itertuples():
            color = C["blue"] if b.strand == 1 else C["red"]
            circ.link(("Mdubia", b.ref_start, b.ref_end), (acc, b.qry_start, b.qry_end),
                      color=color, alpha=0.4 if b.strand == 1 else 0.65, r1=92, r2=92)
        ax = fig.add_subplot(nrow, ncol, i + 1, polar=True)
        circ.plotfig(ax=ax)
        ax.set_title(f"{tribe}\n$\\it{{{species.replace(' ', chr(92)+' ')}}}$\n{cov:.0f}% collinear, {inv} inversion(s)",
                     size=9, pad=10, linespacing=1.4)
    fig.subplots_adjust(wspace=0.08, hspace=0.55, top=0.86, bottom=0.10)
    # legend at the bottom, outside panels
    handles = [Patch(fc=C["blue"], label="Collinear block (same strand)"),
               Patch(fc=C["red"], label="Inverted block (opposite strand)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle("Whole-plastome synteny of $\\it{Myrciaria\\ dubia}$ vs Myrtaceae tribe representatives and tobacco",
                 size=11, y=0.995)
    save(fig, out, "Fig8")
    df = pd.concat(blocks_all, ignore_index=True)
    return df


# ---------------------------------------------------------------------------
# Family-wide synteny conservation
# ---------------------------------------------------------------------------
def plastome_regions(rec):
    regions = []
    for f in rec.features:
        n = " ".join(f.qualifiers.get("note", [])).lower()
        if f.type in ("misc_feature", "repeat_region") and any(k in n for k in ("lsc", "ssc", "ira", "irb", "inverted")):
            lab = "LSC" if "lsc" in n else "SSC" if "ssc" in n else ("IRb" if ("irb" in n or "repeat b" in n) else "IRa")
            regions.append((int(f.location.start) + 1, int(f.location.end), lab))
    return regions


def family_conservation(gdir, refdir, res, work, all_genomes):
    tg = SeqIO.read(str(gdir / "PX251231.1.gb"), "genbank")
    L = len(tg); nwin = int(np.ceil(L / WINDOW))
    write_fasta("Mdubia", str(tg.seq), work / "ref.fa")
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t").set_index("accession")
    status = {}  # acc -> array nwin: 1 collinear, -1 inverted, 0 unaligned
    per_block = []
    for k, (acc, folder) in enumerate(all_genomes):
        if acc == "PX251231.1":
            continue
        rec = SeqIO.read(str(folder / f"{acc}.gb"), "genbank")
        write_fasta(acc, str(rec.seq), work / "qry.fa")
        b = nucmer_blocks(work / "ref.fa", work / "qry.fa", work, "fw")
        arr = np.zeros(nwin, dtype=np.int8)
        for r in b.itertuples():
            w0, w1 = (r.ref_start - 1) // WINDOW, (r.ref_end - 1) // WINDOW
            for w in range(w0, min(w1 + 1, nwin)):
                if arr[w] == 0:
                    arr[w] = r.strand
        status[acc] = arr
        b["accession"] = acc; per_block.append(b.assign(
            tribe=tribe5(meta.loc[acc, "tribe"]) if acc in meta.index else "reference"))
    myr = [a for a, _ in all_genomes if a != "PX251231.1" and a in meta.index and meta.loc[a, "role"] != "reference (validated elements)"]
    M = np.stack([status[a] for a in myr])
    frac_coll = (M == 1).mean(0)
    frac_inv = (M == -1).mean(0)
    frac_unal = (M == 0).mean(0)
    cons = pd.DataFrame(dict(window=np.arange(nwin), start=np.arange(nwin) * WINDOW + 1,
                             frac_collinear=frac_coll, frac_inverted=frac_inv, frac_unaligned=frac_unal))
    # per tribe collinearity
    tribe_rings = {}
    for t in TRIBE_ORDER:
        accs = [a for a in myr if tribe5(meta.loc[a, "tribe"]) == t]
        if accs:
            tribe_rings[t] = (np.stack([status[a] for a in accs]) == 1).mean(0)
    return tg, L, cons, tribe_rings, pd.concat(per_block, ignore_index=True), status, myr


def fig_family_circos(out, gdir, refdir, res, work, all_genomes):
    from pycirclize import Circos
    tg, L, cons, tribe_rings, blocks, status, myr = family_conservation(gdir, refdir, res, work, all_genomes)
    regions = plastome_regions(tg)
    circ = Circos({"Mdubia": L}, start=3, end=357)
    sec = circ.get_sector("Mdubia")
    # gene track
    gt = sec.add_track((90, 99)); gt.axis(fc="#f2f1ec", ec="none")
    fcol = gene_colormap()
    for f in tg.features:
        if f.type not in ("CDS", "tRNA", "rRNA"):
            continue
        g = f.qualifiers.get("gene", [""])[0]
        for p in f.location.parts:
            gt.rect(int(p.start) + 1, int(p.end), r_lim=(95, 99) if p.strand == 1 else (90, 94),
                    fc=fcol(g, f.type), ec="white", lw=0.1)
    # region ring
    rt = sec.add_track((86, 89))
    rc = {"LSC": C["blue"], "SSC": C["aqua"], "IRa": C["orange"], "IRb": C["orange"]}
    for s, e, lab in regions:
        rt.rect(s, e, fc=rc.get(lab, C["muted"]), ec="white", lw=0.3)
        mid = (s + e) / 2
        sec.text(lab, r=84.5, size=8, color=C["ink"], adjust_rotation=False) if False else None
    # stacked conservation ring (collinear / inverted / unaligned) drawn as per-window rects
    R0, R1 = 70, 85; span = R1 - R0
    st = sec.add_track((R0, R1)); st.axis(fc="white", ec=C["axis"], lw=0.4)
    for r in cons.itertuples():
        a0 = r.start; b0 = min(r.start + WINDOW, L)
        if a0 >= b0:
            continue
        fc = r.frac_collinear; fi = r.frac_inverted
        def seg(lo, hi, color, alpha=1.0):
            if hi - lo > 1e-4:
                st.rect(a0, b0, r_lim=(R0 + span * lo, R0 + span * hi), fc=color, ec="none", alpha=alpha)
        seg(0, fc, C["blue"], 0.8)
        seg(fc, fc + fi, "#c9c7bf")
        seg(fc + fi, 1.0, C["red"], 0.5)
    # per-tribe collinearity heat rings (outer -> inner follows TRIBE_ORDER)
    r0 = 68
    cmap = plt.get_cmap("Blues")
    ring_order = []
    for t in TRIBE_ORDER:
        if t not in tribe_rings:
            continue
        tr = sec.add_track((r0 - 4, r0)); r0 -= 5
        ring_order.append(t)
        vals = tribe_rings[t]
        for w, v in enumerate(vals):
            a0 = w * WINDOW + 1; b0 = min((w + 1) * WINDOW, L)
            if a0 >= b0:
                continue
            tr.rect(a0, b0, fc=cmap(0.15 + 0.85 * v), ec="none")
    sec.text(f"$\\it{{Myrciaria\\ dubia}}$\nPX251231.1\n{L:,} bp\nsynteny vs\n{len(myr)} Myrtaceae",
             r=14, size=11, color=C["ink"])
    fig = circ.plotfig(dpi=100)
    fig.set_size_inches(185 * MM, 225 * MM)
    ax = fig.axes[0]
    ax.set_position([0.02, 0.01, 0.96, 0.80])   # circle in lower 80 %, header band on top
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize
    fig.text(0.5, 0.985, "Family-wide plastome synteny conservation across Myrtaceae",
             ha="center", va="top", size=12, fontweight="bold")
    h1 = [Patch(fc=C["blue"], alpha=0.8, label="Collinear (same strand)"),
          Patch(fc="#c9c7bf", label="Inverted (opposite strand)"),
          Patch(fc=C["red"], alpha=0.5, label="Unaligned")]
    leg1 = fig.legend(handles=h1, loc="upper left", bbox_to_anchor=(0.04, 0.955), frameon=False,
                      fontsize=9, title="Outer stacked ring: conservation (per 500 bp)", title_fontsize=9)
    leg1._legend_box.align = "left"
    h2 = [Patch(fc=cmap(0.85), label=f"{i+1}. {t}") for i, t in enumerate(ring_order)]
    leg2 = fig.legend(handles=h2, loc="upper right", bbox_to_anchor=(0.96, 0.955), frameon=False,
                      fontsize=9, title="Inner heat rings (outer→inner): tribe collinearity", title_fontsize=9)
    leg2._legend_box.align = "left"
    sm = cm.ScalarMappable(norm=Normalize(0, 1), cmap=cmap)
    cax = fig.add_axes([0.40, 0.90, 0.22, 0.011])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("Tribe collinearity (fraction)", fontsize=8.5)
    cb.ax.tick_params(labelsize=7)
    save(fig, out, "Fig9")
    return cons, blocks


def gene_colormap():
    fcol = {"photosystem": C["green"], "ribosome": C["yellow"], "rna_pol": C["red"], "atp": C["violet"],
            "ndh": C["magenta"], "rubisco": C["aqua"], "trna": C["ink2"], "rrna": C["orange"], "other": C["blue"]}
    def cl(g, t):
        if t == "tRNA": return fcol["trna"]
        if t == "rRNA": return fcol["rrna"]
        if g.startswith(("psa", "psb", "pet", "paf", "pbf", "ycf3", "ycf4")): return fcol["photosystem"]
        if g.startswith(("rps", "rpl")): return fcol["ribosome"]
        if g.startswith("rpo"): return fcol["rna_pol"]
        if g.startswith("atp"): return fcol["atp"]
        if g.startswith("ndh"): return fcol["ndh"]
        if g.startswith("rbc"): return fcol["rubisco"]
        return fcol["other"]
    return cl


# ---------------------------------------------------------------------------
# Tribe-level summary figure (FigS6)
# ---------------------------------------------------------------------------
def fig_tribe_summary(out, res, blocks_family, status, myr):
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t").set_index("accession")
    rows = []
    for acc in myr:
        b = blocks_family[blocks_family.accession == acc]
        if b.empty:
            continue
        tribe = tribe5(meta.loc[acc, "tribe"])
        coll = b[b.strand == 1].eval("ref_end-ref_start").clip(lower=0).sum()
        inv = int((b.strand == -1).sum())
        rows.append(dict(accession=acc, tribe=tribe, n_blocks=len(b), inversions=inv,
                         mean_identity=b.identity.mean(),
                         collinear_fraction=coll / meta.loc["PX251231.1", "length_bp"]))
    df = pd.DataFrame(rows)
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 68 * MM), gridspec_kw=dict(wspace=0.42))
    order = [t for t in TRIBE_ORDER if t in df.tribe.unique()]
    rng = np.random.default_rng(7)
    for ax, col, lab in [(axs[0], "collinear_fraction", "Collinear fraction vs M. dubia"),
                         (axs[1], "mean_identity", "Mean block identity (%)"),
                         (axs[2], "inversions", "Inverted blocks (n)")]:
        for i, t in enumerate(order):
            v = df[df.tribe == t][col].values
            ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=10, color=TRIBE_COL[t], alpha=0.8, linewidth=0)
            if len(v):
                ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=C["ink"], lw=1.2)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=35, ha="right")
        ax.set_ylabel(lab)
    axs[0].set_title("a", loc="left", fontweight="bold")
    axs[1].set_title("b", loc="left", fontweight="bold")
    axs[2].set_title("c", loc="left", fontweight="bold")
    save(fig, out, "FigS6")
    return df


def build_genome_list(res, gdir, refdir):
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t")
    lst = []
    for r in meta.itertuples():
        folder = refdir if r.role == "reference (validated elements)" else gdir
        lst.append((r.accession, folder))
    return lst


REPS = [
    ("Myrteae", "KY392759.1", "Plinia aureana"),
    ("Syzygieae", "OR727280.1", "Syzygium versteegii"),
    ("Eucalypteae", "KC180806.1", "Allosyncarpia ternata"),
    ("Melaleuceae", "MT700493.1", "Melaleuca leucadendra"),
    ("Other Myrtaceae", "MW837774.1", "Xanthostemon chrysanthus"),
    ("Reference", "NC_001879.2", "Nicotiana tabacum"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--genomes-dir", default="genomes")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    res, gdir, out = Path(a.results), Path(a.genomes_dir), Path(a.out)
    refdir = gdir / "ref"
    work = Path(tempfile.mkdtemp(prefix="synteny_"))
    tabs = res / "tables"
    reps = [(t, acc, sp, refdir if t == "Reference" else gdir) for t, acc, sp in REPS]
    blocks_rep = fig_circos_tribes(out, gdir, refdir, res, reps, work)
    blocks_rep.to_csv(tabs / "S25_synteny_blocks_representatives.tsv", sep="\t", index=False)
    all_g = build_genome_list(res, gdir, refdir)
    cons, blocks_family = fig_family_circos(out, gdir, refdir, res, work, all_g)
    cons.to_csv(tabs / "S26_synteny_conservation_windows.tsv", sep="\t", index=False)
    blocks_family.to_csv(tabs / "S27_synteny_blocks_all.tsv", sep="\t", index=False)
    # rebuild status for tribe summary
    status = {}
    for acc in blocks_family.accession.unique():
        status[acc] = None
    meta = pd.read_csv(tabs / "S1_genomes.tsv", sep="\t").set_index("accession")
    myr = [a for a in blocks_family.accession.unique() if a in meta.index and meta.loc[a, "role"] != "reference (validated elements)"]
    summ = fig_tribe_summary(out, res, blocks_family, status, myr)
    summ.to_csv(tabs / "S28_synteny_per_genome_summary.tsv", sep="\t", index=False)
    print("synteny done")


if __name__ == "__main__":
    main()
