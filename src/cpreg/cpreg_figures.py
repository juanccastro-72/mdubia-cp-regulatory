#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpreg_figures.py - main and supplementary figures for the M. dubia plastid
regulatory-element study. Reads the tables written by cpreg_pipeline.py.

Each figure is written as an editable SVG (live Times New Roman text, every mark
a vector object -> CorelDRAW / Illustrator / Inkscape) and a 600 dpi JPEG.
Base font Times New Roman 10 pt; scientific and gene names in italics. Legends
are placed outside the plotting area; canvases are sized to avoid overlap.
Widths follow Springer: 174 mm (double column) / 84 mm (single column).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

from cpreg_style import (C, SEQ_BLUE, TRIBE_ORDER, TRIBE_COL, MM, apply_style, save, panel, it, tribe5, BASE)

apply_style()
W2, W1 = 174 * MM, 84 * MM
MONO = "Liberation Mono"


# ---------------------------------------------------------------------------
# Fig. 1  Workflow
# ---------------------------------------------------------------------------
def fig_workflow(out):
    fig, ax = plt.subplots(figsize=(W2, 118 * MM))
    ax.set_xlim(0, 100); ax.set_ylim(0, 62); ax.axis("off")

    def box(x, y, w, h, title, body, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc="white", ec=col, lw=1.0))
        ax.add_patch(Rectangle((x - 0.4, y + h - 3.6), w + 0.8, 4.0, fc=col, ec="none"))
        ax.text(x + w / 2, y + h - 1.6, title, ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax.text(x + 1.0, y + h - 5.2, body, ha="left", va="top", fontsize=7.5, color=C["ink"], linespacing=1.4)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=1.0, color=C["ink2"], shrinkA=0, shrinkB=0))

    box(1, 43, 30, 18, "1  Input data",
        f"{it('M. dubia')} plastome (PX251231.1)\n121 Myrtaceae plastomes (NCBI)\n12 reference plastomes\nInventory: 48 validated elements\n(cached GenBank, SHA-256)", C["blue"])
    box(35, 43, 30, 18, "2  Validated-element library",
        "45 loci from reference plastomes\n(promoters, 5'UTRs, SD,\nIEE, PPR site, 3'UTRs)\n+ class-level motif models", C["blue"])
    box(69, 43, 30, 18, f"3  Homology to {it('M. dubia')}",
        "BLASTn + synteny check\nAnchored global alignment vs\n500 random spacer windows\n(empirical P, BH-FDR)", C["blue"])
    box(1, 21, 30, 19, "4  De novo prediction",
        "PEP -35/-10 PSSM (spacer 15-19)\nNEP YRTA / GAA-box composite\nSD: anti-SD duplex dG\nTIR and 3' stem-loop MFE\nNull: dinucleotide shuffles", C["orange"])
    box(35, 21, 30, 19, "5  Calibration",
        "Negative controls: convergent\nspacers, coding windows\nPositive control: tobacco\nvalidated PEP promoters\nModel-species SD frequency", C["orange"])
    box(69, 21, 30, 19, "6  Comparative genomics",
        "Ortholog mapping in 121 Myrtaceae\nMotif vs flank conservation\nCross-species re-scoring\nStem-loop covariation\nSynteny + CDS phylogeny", C["aqua"])
    box(6, 1, 88, 16, "7  Evidence integration and outputs",
        "Tier 1: syntenic homology to a validated element, or significant prediction conserved in ≥ 95 % of Myrtaceae\n"
        "Tier 2: significant de novo prediction (FDR < 0.05)      Tier 3: nominal prediction (exploratory)\n"
        "Outputs: Online Resource tables (xlsx), figures (SVG/JPEG), synteny blocks, Newick tree, run metadata", C["violet"])
    arrow(31.5, 52, 34.5, 52); arrow(65.5, 52, 68.5, 52)
    arrow(16, 42.5, 16, 40.5); arrow(84, 42.5, 84, 40.5)
    arrow(31.5, 30.5, 34.5, 30.5); arrow(65.5, 30.5, 68.5, 30.5)
    arrow(16, 20.5, 16, 17.5); arrow(50, 20.5, 50, 17.5); arrow(84, 20.5, 84, 17.5)
    save(fig, out, "Fig1")


# ---------------------------------------------------------------------------
# Fig. 2  Circular plastome map
# ---------------------------------------------------------------------------
def _gene_class(g, t):
    if t == "tRNA": return "trna"
    if t == "rRNA": return "rrna"
    if g.startswith(("psa", "psb", "pet", "paf", "pbf", "ycf3", "ycf4")): return "photosystem"
    if g.startswith(("rps", "rpl")): return "ribosome"
    if g.startswith("rpo"): return "rna_pol"
    if g.startswith("atp"): return "atp"
    if g.startswith("ndh"): return "ndh"
    if g.startswith("rbc"): return "rubisco"
    return "other"


def fig_map(out, res, gdir, target):
    rec = SeqIO.read(str(gdir / f"{target}.gb"), "genbank")
    L = len(rec)
    integ = pd.read_csv(res / "tables" / "S23_integrated_evidence_target.tsv", sep="\t")
    loci = pd.read_csv(res / "tables" / "S3_target_loci.tsv", sep="\t")
    fig = plt.figure(figsize=(W2, 165 * MM))
    ax = fig.add_axes([0.0, 0.0, 0.70, 1.0], projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.axis("off"); ax.set_ylim(0, 10)

    def th(p):
        return 2 * np.pi * p / L

    fcol = {"photosystem": C["green"], "ribosome": C["yellow"], "rna_pol": C["red"], "atp": C["violet"],
            "ndh": C["magenta"], "rubisco": C["aqua"], "trna": C["ink2"], "rrna": C["orange"], "other": C["blue"]}
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA"):
            continue
        g = f.qualifiers.get("gene", [""])[0]
        for p in f.location.parts:
            r0 = 7.7 if p.strand == 1 else 6.7
            ax.bar(th((p.start + p.end) / 2), 0.9, width=max(th(p.end - p.start), 0.004), bottom=r0,
                   color=fcol[_gene_class(g, f.type)], edgecolor="white", linewidth=0.15)
    # gene labels: protein-coding + rRNA only (skip tRNA to declutter), radial, de-collided
    done = set(); placed = []
    feats = [f for f in rec.features if f.type in ("CDS", "rRNA")]
    for f in feats:
        g = f.qualifiers.get("gene", [""])[0]
        if not g or g in done:
            continue
        done.add(g)
        mid = (int(f.location.start) + int(f.location.end)) / 2
        if int(f.location.end) - int(f.location.start) > 20000:
            continue
        a = th(mid); deg = np.degrees(a)
        # avoid collisions by nudging radius when two labels are within ~3.3 deg
        r = 8.85
        for pd_deg, pr in placed:
            if abs(((deg - pd_deg + 180) % 360) - 180) < 3.3 and abs(pr - r) < 0.55:
                r += 0.6
        placed.append((deg, r))
        rot = 90 - deg if deg <= 180 else 270 - deg
        ax.text(a, r, g, rotation=rot, rotation_mode="anchor", ha="left" if deg <= 180 else "right",
                va="center", fontsize=5.4, style="italic", color=C["ink"])
    # LSC / IR / SSC ring
    regions = []
    for f in rec.features:
        n = " ".join(f.qualifiers.get("note", [])).lower()
        if f.type in ("misc_feature", "repeat_region") and any(k in n for k in ("lsc", "ssc", "ira", "irb", "inverted")):
            lab = "LSC" if "lsc" in n else "SSC" if "ssc" in n else ("IRb" if ("irb" in n or "inverted repeat b" in n) else "IRa")
            regions.append((int(f.location.start), int(f.location.end), lab))
    rc = {"LSC": C["blue"], "SSC": C["aqua"], "IRa": C["orange"], "IRb": C["orange"]}
    for s, e, lab in regions:
        ax.bar(th((s + e) / 2), 0.45, width=th(e - s), bottom=5.95, color=rc[lab], edgecolor="white", linewidth=0.5)
        ax.text(th((s + e) / 2), 5.3, f"{lab}\n{(e - s) / 1000:.1f} kb", ha="center", va="center", fontsize=7, color=C["ink"])
    loci = loci.set_index("gene")
    tracks = [("promoter_tier", 4.3, C["red"]), ("SD_tier", 3.5, C["violet"]), ("stemloop_tier", 2.7, C["green"])]
    for col, r, colr in tracks:
        ax.bar(np.linspace(0, 2 * np.pi, 361), 0.55, width=2 * np.pi / 360, bottom=r, color=C["neutral"], edgecolor="none")
        for x in integ.itertuples():
            t = getattr(x, col, np.nan)
            if pd.isna(t) or x.gene not in loci.index:
                continue
            pos = loci.loc[x.gene]
            pos = pos.iloc[0] if isinstance(pos, pd.DataFrame) else pos
            p = pos.five_prime_coord if col != "stemloop_tier" else pos.three_prime_coord
            size = {1: 0.55, 2: 0.38, 3: 0.2}[int(t)]
            ax.bar(th(p), size, width=0.012, bottom=r, color=colr, alpha={1: 1, 2: 0.75, 3: 0.4}[int(t)], edgecolor="none")
    hom = pd.read_csv(res / "tables" / "S6_homology_validated_elements.tsv", sep="\t")
    ax.bar(np.linspace(0, 2 * np.pi, 361), 0.55, width=2 * np.pi / 360, bottom=1.9, color=C["neutral"], edgecolor="none", zorder=0)
    for x in hom[hom.homology_supported & hom.Mdubia_ortholog_coord.notna()].itertuples():
        a, b = [int(v) for v in x.Mdubia_ortholog_coord.split()[0].split("-")]
        ax.bar(th((a + b) / 2), 0.55, width=max(th(b - a), 0.01), bottom=1.9, color=C["blue"], edgecolor="none")
    ax.text(0, 0, f"{it('Myrciaria dubia')}\nplastome\n{L:,} bp\nPX251231.1", ha="center", va="center", fontsize=10,
            fontweight="bold")
    # legend (own axes, right)
    lax = fig.add_axes([0.71, 0.08, 0.28, 0.84]); lax.axis("off")
    y = 0.99
    lax.text(0, y, "Genes (outer: + strand; inner: − strand)", fontsize=9, fontweight="bold", va="top"); y -= 0.05
    items = [("Photosynthesis (psa/psb/pet)", fcol["photosystem"]), ("Ribosomal proteins", fcol["ribosome"]),
             ("RNA polymerase (rpo)", fcol["rna_pol"]), ("ATP synthase (atp)", fcol["atp"]), ("NDH complex (ndh)", fcol["ndh"]),
             ("Rubisco (rbcL)", fcol["rubisco"]), ("Other genes", fcol["other"]), ("rRNA", fcol["rrna"]), ("tRNA", fcol["trna"])]
    for lab, c in items:
        lax.add_patch(Rectangle((0, y - 0.02), 0.08, 0.03, color=c, transform=lax.transAxes)); lax.text(0.12, y - 0.005, lab, fontsize=8, va="center"); y -= 0.042
    y -= 0.03
    lax.text(0, y, "Regulatory tracks (outer→inner)", fontsize=9, fontweight="bold", va="top"); y -= 0.05
    for lab, c in [("PEP promoter (5' spacer)", C["red"]), ("Shine–Dalgarno motif", C["violet"]),
                   ("3' stem-loop", C["green"]), ("Homolog of validated element", C["blue"])]:
        lax.add_patch(Rectangle((0, y - 0.02), 0.08, 0.03, color=c, transform=lax.transAxes)); lax.text(0.12, y - 0.005, lab, fontsize=8, va="center"); y -= 0.042
    y -= 0.02
    lax.text(0, y, "Bar height/opacity: Tier 1 > 2 > 3", fontsize=8, va="top", color=C["ink2"])
    save(fig, out, "Fig2")


# ---------------------------------------------------------------------------
# Fig. 3  Homology to validated elements
# ---------------------------------------------------------------------------
def fig_homology(out, res, gdir, refdir, target):
    hom = pd.read_csv(res / "tables" / "S6_homology_validated_elements.tsv", sep="\t")
    hom = hom[hom.anchored_identity.notna()].copy()
    order = hom.sort_values(["element_class", "anchored_identity"]).element_id.tolist()
    hom = hom.set_index("element_id").loc[order].reset_index()
    fig = plt.figure(figsize=(W2, 200 * MM))
    ax = fig.add_axes([0.26, 0.34, 0.34, 0.63])
    supp_col = {"strong (BLAST + anchored alignment)": C["blue"], "BLAST only": C["aqua"],
                "anchored alignment only": C["yellow"], "none": C["muted"]}
    for i, r in enumerate(hom.itertuples()):
        ax.plot([r.null_identity_mean, r.null_identity_p99], [i, i], color=C["axis"], lw=3, solid_capstyle="round", zorder=1)
        ax.plot([r.null_identity_p99, r.anchored_identity], [i, i], color=C["grid"], lw=0.8, zorder=1)
        ax.scatter(r.anchored_identity, i, s=28, color=supp_col[r.homology_support], edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_yticks(range(len(hom))); ax.set_yticklabels([r.element_id for r in hom.itertuples()], fontsize=7)
    ax.set_xlabel(f"Identity with {it('M. dubia')} orthologous region")
    ax.set_xlim(0.5, 1.0); ax.grid(axis="x", color=C["grid"], lw=0.4); ax.set_axisbelow(True)
    ax.set_ylim(-0.7, len(hom) - 0.3)
    panel(ax, "a", x=-0.52, y=1.0)
    handles = [Line2D([], [], marker="o", ls="", color=c, label=k, markersize=6) for k, c in supp_col.items()]
    handles.append(Line2D([], [], color=C["axis"], lw=3, label="Null: mean – 99th percentile"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=1, frameon=False, fontsize=8,
              title="Homology support", title_fontsize=8.5)
    # b: class summary
    ax2 = fig.add_axes([0.70, 0.73, 0.28, 0.24])
    cl = hom.groupby("element_class").agg(n=("element_id", "size"), sup=("homology_supported", "sum")).sort_values("n")
    ax2.barh(cl.index, cl.n, color=C["neutral"], edgecolor=C["axis"], lw=0.5, height=0.65, label="Tested")
    ax2.barh(cl.index, cl.sup, color=C["blue"], height=0.65, label="Supported")
    for i, (n, s) in enumerate(zip(cl.n, cl.sup)):
        ax2.text(n + 0.2, i, f"{int(s)}/{int(n)}", va="center", fontsize=7.5)
    ax2.set_xlabel("Reference elements"); ax2.set_xlim(0, cl.n.max() + 3)
    ax2.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax2.tick_params(labelsize=7.5)
    panel(ax2, "b", x=-0.62)
    # c: psbA promoter alignment
    ax3 = fig.add_axes([0.66, 0.06, 0.33, 0.58]); ax3.axis("off"); panel(ax3, "c", x=-0.06, y=1.0)
    from cpreg_pipeline import Plastome, make_aligner
    seqs = {}
    for lab, acc, path in [("N. tabacum", "NC_001879.2", refdir), ("A. thaliana", "NC_000932.1", refdir),
                           ("H. vulgare", "NC_008590.1", refdir), ("L. sativa", "NC_007578.1", refdir),
                           ("M. truncatula", "NC_003119.8", refdir), ("M. dubia", target, gdir)]:
        pg = Plastome(path / f"{acc}.gb")
        l = pg.locus("psbA", "CDS")
        seqs[lab] = pg.upstream(l.five, l.strand, 160)
    ref = seqs["N. tabacum"]; al = make_aligner()
    i35 = ref.find("TTGACA"); i10 = ref.find("TATACT", i35)
    a0, a1 = i35 - 8, i10 + 22
    ax3.text(0, 1.0, f"{it('psbA')} promoter (−35 / −10) and 5' leader", fontsize=9, fontweight="bold", va="top", transform=ax3.transAxes)
    yy = 0.9; dx = 0.0150
    for lab, s in seqs.items():
        aln = al.align(ref, s)[0]
        qmap = {}
        for (ts, te), (qs, qe) in zip(*aln.aligned):
            for k in range(te - ts):
                qmap[ts + k] = qs + k
        row = "".join(s[qmap[k]] if k in qmap else "-" for k in range(a0, a1))
        ax3.text(0, yy, it(lab), fontsize=8, va="center", transform=ax3.transAxes,
                 fontweight="bold" if lab == "M. dubia" else "normal")
        for k, ch in enumerate(row):
            pos = a0 + k
            inbox = (i35 <= pos < i35 + 6) or (i10 <= pos < i10 + 6)
            same = ch == ref[pos]
            xx = 0.30 + k * dx
            if inbox:
                ax3.add_patch(Rectangle((xx - 0.007, yy - 0.035), dx, 0.07, transform=ax3.transAxes,
                                        color="#fde3d6" if same else "#fbd0d0", lw=0))
            ax3.text(xx, yy, ch, fontsize=6.5, family=MONO, ha="center", va="center", transform=ax3.transAxes,
                     color=C["ink"] if same else C["red"])
        yy -= 0.11
    ax3.text(0.30 + (i35 - a0 + 2.5) * dx, yy + 0.03, "−35", ha="center", fontsize=8, fontweight="bold", color=C["orange"], transform=ax3.transAxes)
    ax3.text(0.30 + (i10 - a0 + 2.5) * dx, yy + 0.03, "−10", ha="center", fontsize=8, fontweight="bold", color=C["orange"], transform=ax3.transAxes)
    ax3.text(0, yy - 0.08, f"Red: differences from {it('N. tabacum')}. Shaded: hexamers of the\nvalidated tobacco PpsbA (Staub and Maliga 1993).",
             fontsize=7.5, color=C["ink2"], transform=ax3.transAxes, va="top")
    save(fig, out, "Fig3")


# ---------------------------------------------------------------------------
# Fig. 4  PEP / NEP promoter predictions
# ---------------------------------------------------------------------------
def draw_logo(ax, seqs, title):
    seqs = [s for s in seqs if len(s) == 6 and set(s) <= set("ACGT")]
    n = len(seqs)
    cols = {"A": C["green"], "C": C["blue"], "G": C["yellow"], "T": C["red"]}
    for i in range(6):
        cnt = {b: sum(s[i] == b for s in seqs) / n for b in "ACGT"}
        H = -sum(p * np.log2(p) for p in cnt.values() if p > 0)
        R = max(2 - (H + 3 / (2 * np.log(2) * n)), 0)
        y = 0
        for b, p in sorted(cnt.items(), key=lambda x: x[1]):
            h = p * R
            if h < 0.02:
                continue
            t = ax.text(i + 0.5, y, b, ha="center", va="bottom", fontweight="bold", color=cols[b], family=MONO)
            t.set_fontsize(max(2, 26 * h)); y += h
    ax.set_xlim(0, 6); ax.set_ylim(0, 2); ax.set_xticks(np.arange(6) + 0.5); ax.set_xticklabels(range(1, 7))
    ax.set_ylabel("Bits"); ax.set_title(f"{title} (n = {n})", fontsize=9)


def _box_strip(ax, groups, rng, ylab):
    for i, (lab, v, c) in enumerate(groups):
        v = np.asarray(v)
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=C["ink"], lw=1),
                   boxprops=dict(color=C["ink2"], lw=0.6), whiskerprops=dict(color=C["ink2"], lw=0.6), capprops=dict(lw=0))
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=c, alpha=0.7, linewidth=0, zorder=3)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups], fontsize=8)
    ax.set_ylabel(ylab)


def fig_promoters(out, res):
    pep = pd.read_csv(res / "tables" / "S7_PEP_promoters.tsv", sep="\t")
    ctrl = pd.read_csv(res / "tables" / "S8_PEP_negative_controls.tsv", sep="\t")
    nep = pd.read_csv(res / "tables" / "S9_NEP_motifs.tsv", sep="\t")
    nctrl = pd.read_csv(res / "tables" / "S9b_NEP_controls.tsv", sep="\t")
    mf = pd.read_csv(res / "tables" / "S17_motif_vs_flank_conservation.tsv", sep="\t")
    fig, axs = plt.subplots(2, 3, figsize=(W2, 125 * MM), gridspec_kw=dict(wspace=0.5, hspace=0.62))
    rng = np.random.default_rng(1)
    # a
    ax = axs[0, 0]
    _box_strip(ax, [("5' spacers", pep.z, C["blue"]),
                    ("Convergent", ctrl[ctrl.control == "convergent spacer"].z, C["muted"]),
                    ("Coding", ctrl[ctrl.control == "coding sequence"].z, C["muted"])], rng, "PEP score z (vs shuffles)")
    top = pep.nsmallest(4, "p_empirical").sort_values("z")
    for r, dy in zip(top.itertuples(), np.linspace(-1.0, 1.0, len(top))):
        ax.annotate(it(r.id), (0.14, r.z), xytext=(0.42, r.z + dy), fontsize=7, va="center",
                    color=C["ink"], arrowprops=dict(arrowstyle="-", lw=0.3, color=C["muted"]))
    ax.axhline(0, color=C["axis"], lw=0.5); panel(ax, "a", x=-0.34)
    # b
    ax = axs[0, 1]
    for lab, v, c in [("5' spacers", pep.p_empirical, C["blue"]), ("Controls", ctrl.p_empirical, C["muted"])]:
        v = np.sort(np.asarray(v)); n = len(v); e = (np.arange(1, n + 1) - 0.5) / n
        ax.scatter(-np.log10(e), -np.log10(v), s=7, color=c, label=lab, linewidth=0)
    ax.plot([0, 3.2], [0, 3.2], color=C["ink2"], lw=0.6, ls="--")
    ax.set_xlabel("Expected −log₁₀ P"); ax.set_ylabel("Observed −log₁₀ P")
    ax.legend(frameon=False, loc="upper left"); panel(ax, "b", x=-0.34)
    # c
    ax = axs[0, 2]
    for r in mf.itertuples():
        sig = r.q_BH < 0.10
        ax.plot([0, 1], [r.flank_conservation, r.motif_conservation], color=C["blue"] if sig else C["grid"],
                lw=0.8 if sig else 0.5, zorder=3 if sig else 1)
    ax.scatter(np.zeros(len(mf)), mf.flank_conservation, s=6, color=C["muted"], zorder=4)
    ax.scatter(np.ones(len(mf)), mf.motif_conservation, s=6, color=C["blue"], zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Flanking\nspacer", "−35/−10\nhexamers"]); ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel("Positional identity")
    st = json.loads((res / "statistics.json").read_text()).get("motif_vs_flank_all_windows", {})
    if st:
        ax.set_title(f"Wilcoxon P = {st['p']:.1e}", fontsize=8, loc="center", fontweight="normal")
    panel(ax, "c", x=-0.34)
    # d, e
    topn = pep.nsmallest(15, "p_empirical")
    draw_logo(axs[1, 0], topn.minus35.tolist(), "−35 box, top 15"); panel(axs[1, 0], "d", x=-0.34)
    draw_logo(axs[1, 1], topn.minus10.tolist(), "−10 box, top 15"); panel(axs[1, 1], "e", x=-0.34)
    # f
    _box_strip(axs[1, 2], [("5' spacers", nep.YRTA_per_100nt, C["blue"]),
                           ("Convergent", nctrl[nctrl.control == "convergent spacer"].YRTA_per_100nt, C["muted"]),
                           ("Coding", nctrl[nctrl.control == "coding sequence"].YRTA_per_100nt, C["muted"])],
               rng, "YRTA per 100 nt (NEP core)")
    panel(axs[1, 2], "f", x=-0.34)
    save(fig, out, "Fig4")


# ---------------------------------------------------------------------------
# Fig. 5  SD and TIR
# ---------------------------------------------------------------------------
def fig_sd(out, res):
    prof = pd.read_csv(res / "work" / "sd_profile_target.tsv", sep="\t")
    sd = pd.read_csv(res / "tables" / "S10_SD_target.tsv", sep="\t")
    sdr = pd.read_csv(res / "tables" / "S14_SD_reference_species.tsv", sep="\t")
    summ = pd.read_csv(res / "tables" / "S24_species_summary.tsv", sep="\t")
    tir = pd.read_csv(res / "tables" / "S11_TIR_structure.tsv", sep="\t")
    fig, axs = plt.subplots(1, 3, figsize=(W2, 66 * MM), gridspec_kw=dict(wspace=0.5, width_ratios=[1.25, 1.1, 1]))
    ax = axs[0]
    g = prof.groupby("offset").dG; mu, se = g.mean(), g.sem()
    ax.fill_between(mu.index, mu - 1.96 * se, mu + 1.96 * se, color=SEQ_BLUE[2], lw=0)
    ax.plot(mu.index, mu.values, color=C["blue"], lw=1.4)
    ax.axvspan(-20, -1, color=C["neutral"], zorder=0)
    ax.axvline(0, color=C["ink2"], lw=0.5, ls="--")
    ax.text(1, mu.min() - 0.05, "AUG", fontsize=8, color=C["ink2"])
    ax.text(-10.5, mu.max() + 0.03, "SD window", fontsize=8, color=C["ink2"], va="bottom", ha="center")
    ax.set_xlabel("Window start vs start codon (nt)"); ax.set_ylabel("Mean anti-SD ΔG (kcal/mol)")
    panel(ax, "a", x=-0.3)
    ax = axs[1]
    labels = ["M. dubia", "Myrtaceae", "N. tabacum", "A. thaliana", "H. vulgare", "C. reinhardtii"]
    myr = summ[summ.accession != "PX251231.1"].SD_fraction.values
    vals = [sd.SD_present.mean(), None]
    for org in ["Nicotiana tabacum", "Arabidopsis thaliana", "Hordeum vulgare subsp. vulgare", "Chlamydomonas reinhardtii"]:
        d = sdr[sdr.species == org]; vals.append(d.SD_present.mean() if len(d) else np.nan)
    for i, v in enumerate(vals):
        if i == 1:
            ax.boxplot(myr, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=C["ink"]),
                       boxprops=dict(color=C["ink2"], lw=0.6), whiskerprops=dict(color=C["ink2"], lw=0.6), capprops=dict(lw=0))
            ax.scatter(i + np.random.default_rng(2).uniform(-0.15, 0.15, len(myr)), myr, s=5, color=C["aqua"], linewidth=0, alpha=0.7)
        else:
            ax.bar(i, v, width=0.6, color=C["blue"] if i == 0 else C["muted"])
            ax.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=7.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([it(x) if i != 1 else f"{x}\n(n=121)" for i, x in enumerate(labels)], rotation=40, ha="right", fontsize=7.5)
    ax.set_ylabel("Fraction of CDS with SD"); ax.set_ylim(0, 1); panel(ax, "b", x=-0.32)
    ax = axs[2]
    ax.hist(tir.z, bins=np.arange(-4, 4.01, 0.5), color=C["blue"], edgecolor="white", lw=0.5)
    ax.axvline(0, color=C["ink2"], lw=0.6, ls="--")
    ax.set_xlabel("TIR folding z-score"); ax.set_ylabel("Genes")
    ax.annotate("more structured", xy=(-3.6, 0.5), fontsize=7, color=C["ink2"], rotation=90, va="center")
    panel(ax, "c", x=-0.3)
    save(fig, out, "Fig5")


# ---------------------------------------------------------------------------
# Fig. 6  3' stem-loops
# ---------------------------------------------------------------------------
def fig_stemloops(out, res):
    sl = pd.read_csv(res / "tables" / "S12_3prime_stemloops.tsv", sep="\t").sort_values("z")
    cov = pd.read_csv(res / "tables" / "S20_stemloop_covariation.tsv", sep="\t")
    fig, axs = plt.subplots(1, 2, figsize=(W2, 100 * MM), gridspec_kw=dict(wspace=0.32, width_ratios=[1.35, 1]))
    ax = axs[0]
    cols = [C["green"] if q < 0.05 else (SEQ_BLUE[4] if p < 0.05 else C["axis"]) for q, p in zip(sl.q_BH, sl.p_empirical)]
    ax.bar(range(len(sl)), sl.z, color=cols, width=0.8)
    ax.set_xticks([]); ax.set_xlabel("Protein-coding genes ranked by z-score")
    for i, r in enumerate(sl.itertuples()):
        if r.p_empirical < 0.05:
            ax.text(i, r.z - 0.3, r.gene, rotation=90, ha="center", va="top", fontsize=6.5, style="italic")
    ax.set_ylabel("3' stem-loop MFE z-score"); ax.axhline(0, color=C["ink2"], lw=0.5)
    ax.axhline(-2, color=C["ink2"], lw=0.5, ls=":"); ax.set_xlim(-1, len(sl))
    ax.legend(handles=[Patch(color=C["green"], label="FDR < 0.05"), Patch(color=SEQ_BLUE[4], label="P < 0.05"),
                       Patch(color=C["axis"], label="n.s.")], frameon=False, loc="lower right", fontsize=8)
    panel(ax, "a", x=-0.13)
    ax = axs[1]
    top = cov.merge(sl[["gene", "p_empirical"]], on="gene").sort_values("z_target").head(15)
    cats = [("identical", C["blue"], "Identical pair"), ("compensatory", C["green"], "Compensatory (double)"),
            ("consistent_single", C["aqua"], "Consistent (single)"), ("disrupted", C["red"], "Disrupted"),
            ("gapped", C["axis"], "Gap")]
    left = np.zeros(len(top)); tot = top[[c for c, _, _ in cats]].sum(1).values
    for c, col, lab in cats:
        v = top[c].values / np.maximum(tot, 1)
        ax.barh(range(len(top)), v, left=left, color=col, label=lab, height=0.75, edgecolor="white", lw=0.3)
        left += v
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, style="italic", fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("Fraction of base pairs across 121 Myrtaceae"); ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    panel(ax, "b", x=-0.26)
    save(fig, out, "Fig6")


# ---------------------------------------------------------------------------
# Phylogeny helpers
# ---------------------------------------------------------------------------
def tree_layout(tree):
    tips = tree.get_terminals()
    ypos = {t: i for i, t in enumerate(tips)}; xpos = {}

    def setx(cl, x):
        xpos[cl] = x
        for c in cl.clades:
            setx(c, x + (c.branch_length or 0))
    setx(tree.root, 0)

    def sety(cl):
        if cl.is_terminal():
            return ypos[cl]
        ys = [sety(c) for c in cl.clades]; ypos[cl] = (min(ys) + max(ys)) / 2; return ypos[cl]
    sety(tree.root)
    return xpos, ypos, tips


def draw_tree(ax, tree, xpos, ypos, lw=0.5, support=False):
    for cl in tree.find_clades():
        for c in cl.clades:
            ax.plot([xpos[cl], xpos[cl]], [ypos[cl], ypos[c]], color=C["ink2"], lw=lw)
            ax.plot([xpos[cl], xpos[c]], [ypos[c], ypos[c]], color=C["ink2"], lw=lw)
        if support and not cl.is_terminal() and cl.confidence is not None and len(cl.clades) and cl.confidence < 0.9:
            ax.text(xpos[cl], ypos[cl], f"{cl.confidence:.2f}", fontsize=4.5, ha="right", va="bottom", color=C["red"])


def load_tree(res, meta):
    tree = Phylo.read(str(res / "myrtaceae_cds_fasttree.nwk"), "newick")
    outg = [t for t in tree.get_terminals() if t.name and t.name in set(meta[meta.genus == "Heteropyxis"].accession)]
    if outg:
        tree.root_with_outgroup(outg[0])
    tree.ladderize(reverse=True)
    return tree


# ---------------------------------------------------------------------------
# Fig. 7  Phylogeny + conservation heatmap
# ---------------------------------------------------------------------------
def fig_conservation(out, res, target):
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t")
    meta = meta[meta.role != "reference (validated elements)"].set_index("accession")
    ident = pd.read_csv(res / "tables" / "S15_region_identity_long.tsv", sep="\t")
    integ = pd.read_csv(res / "tables" / "S23_integrated_evidence_target.tsv", sep="\t")
    tree = load_tree(res, meta.reset_index())
    xpos, ypos, tips = tree_layout(tree)
    sel = integ[(integ.promoter_tier <= 2) | (integ.SD_tier <= 1) | (integ.stemloop_tier <= 1) |
                integ.validated_5prime_homologs.notna() | integ.validated_3prime_homologs.notna()]
    cols5 = [g for g in sel.gene if g in set(ident[ident.region == "UP"].gene)][:38]
    cols3 = [g for g in sel[sel.type == "CDS"].gene if g in set(ident[ident.region == "DN"].gene)][:28]
    up = ident[ident.region == "UP"].pivot(index="accession", columns="gene", values="identity")
    dn = ident[ident.region == "DN"].pivot(index="accession", columns="gene", values="identity")
    names = [t.name for t in tips]
    M = np.hstack([up.reindex(index=names, columns=cols5).values, np.full((len(names), 1), np.nan),
                   dn.reindex(index=names, columns=cols3).values])
    fig = plt.figure(figsize=(W2, 250 * MM))
    axt = fig.add_axes([0.01, 0.07, 0.13, 0.86]); draw_tree(axt, tree, xpos, ypos, lw=0.4)
    axt.set_ylim(len(tips) - 0.5, -0.5); axt.axis("off")
    axl = fig.add_axes([0.14, 0.07, 0.20, 0.86]); axl.set_ylim(len(tips) - 0.5, -0.5); axl.set_xlim(0, 1); axl.axis("off")
    for i, a in enumerate(names):
        tr = tribe5(meta.loc[a, "tribe"]) if a in meta.index else "Other Myrtaceae"
        axl.add_patch(Rectangle((0, i - 0.45), 0.045, 0.9, color=TRIBE_COL[tr], lw=0))
        sp = meta.loc[a, "species"] if a in meta.index else a
        axl.text(0.07, i, sp, fontsize=4.3, style="italic", va="center",
                 fontweight="bold" if a == target else "normal", color=C["red"] if a == target else C["ink"])
    ax = fig.add_axes([0.35, 0.07, 0.60, 0.86])
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(SEQ_BLUE); cmap.set_bad("white")
    bounds = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 1.0001]
    im = ax.imshow(M, aspect="auto", cmap=cmap, norm=BoundaryNorm(bounds, cmap.N), interpolation="none")
    labs = [f"{g} 5'" for g in cols5] + [""] + [f"{g} 3'" for g in cols3]
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=90, fontsize=5, style="italic")
    ax.xaxis.tick_top(); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    cax = fig.add_axes([0.45, 0.038, 0.32, 0.009])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label(f"Identity with {it('M. dubia')} region", fontsize=8)
    cb.ax.tick_params(labelsize=6.5)
    axg = fig.add_axes([0.02, 0.015, 0.42, 0.03]); axg.axis("off")
    xpos_leg = [0.0, 0.17, 0.37, 0.60, 0.80]
    for xp, t in zip(xpos_leg, TRIBE_ORDER):
        axg.add_patch(Rectangle((xp, 0.3), 0.022, 0.4, color=TRIBE_COL[t], transform=axg.transAxes))
        axg.text(xp + 0.03, 0.5, t, fontsize=6.5, va="center", transform=axg.transAxes)
    save(fig, out, "Fig7")


# ---------------------------------------------------------------------------
# Supplementary figures
# ---------------------------------------------------------------------------
def figS1_spacers(out, res):
    igs = pd.read_csv(res / "tables" / "S2_target_spacers.tsv", sep="\t")
    igs = igs[~igs.IR_duplicate & (igs.orientation != "undetermined")]
    order = ["divergent", "tandem(+)", "tandem(-)", "convergent"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 68 * MM), gridspec_kw=dict(wspace=0.3))
    rng = np.random.default_rng(3)
    for ax, col, lab in [(axs[0], "length", "Spacer length (nt)"), (axs[1], "GC", "Spacer GC fraction")]:
        for i, o in enumerate(order):
            v = igs[igs.orientation == o][col].values
            ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=C["ink"]),
                       boxprops=dict(color=C["ink2"], lw=0.6), whiskerprops=dict(color=C["ink2"], lw=0.6), capprops=dict(lw=0))
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=7, color=C["blue"], alpha=0.7, linewidth=0)
            ax.text(i, np.nanmax(igs[col]) * 1.04, f"n={len(v)}", ha="center", fontsize=7)
        ax.set_xticks(range(4)); ax.set_xticklabels(order, fontsize=8); ax.set_ylabel(lab)
    panel(axs[0], "a", x=-0.2); panel(axs[1], "b", x=-0.2)
    save(fig, out, "FigS1")


def figS2_tobacco(out, res):
    t = pd.read_csv(res / "tables" / "S13_PEP_tobacco_positive_control.tsv", sep="\t").sort_values("p_empirical").reset_index(drop=True)
    known = {"psbA", "rbcL", "rrn16", "atpB"}
    fig, ax = plt.subplots(figsize=(W1 * 1.2, 72 * MM))
    ax.scatter(t.index + 1, -np.log10(t.p_empirical), s=8, color=C["muted"], linewidth=0)
    k = t[t.id.isin(known)]
    ax.scatter(k.index + 1, -np.log10(k.p_empirical), s=24, color=C["red"], zorder=3)
    for r in k.itertuples():
        ax.annotate(f"{r.id} ({r.minus35}/{r.minus10})", (r.Index + 1, -np.log10(r.p_empirical)), xytext=(7, 2),
                    textcoords="offset points", fontsize=8, style="italic")
    ax.set_xlabel("Rank of tobacco 5' spacer"); ax.set_ylabel("−log₁₀ empirical P (PEP model)")
    save(fig, out, "FigS2")


def figS3_sd_conservation(out, res, target):
    sdc = pd.read_csv(res / "tables" / "S18_SD_conservation_long.tsv", sep="\t")
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t")
    meta = meta[meta.role != "reference (validated elements)"]
    tree = load_tree(res, meta)
    _, _, tips = tree_layout(tree)
    names = [t.name for t in tips]
    M = sdc.pivot_table(index="accession", columns="gene", values="dG").reindex(index=names)
    genes = M.columns[np.argsort(M.mean().values)]; M = M[genes]
    fig, ax = plt.subplots(figsize=(W2, 165 * MM))
    im = ax.imshow(M.values, aspect="auto", cmap="Blues_r", vmin=-12, vmax=0, interpolation="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=5, style="italic")
    sp = meta.set_index("accession").species
    ax.set_yticks(range(len(names))); ax.set_yticklabels([sp.get(a, a) for a in names], fontsize=3.4, style="italic")
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label("Anti-SD duplex ΔG (kcal/mol)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    save(fig, out, "FigS3")


def figS4_tree(out, res, target):
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t")
    meta = meta[meta.role != "reference (validated elements)"]
    tree = load_tree(res, meta)
    xpos, ypos, tips = tree_layout(tree)
    fig, ax = plt.subplots(figsize=(W2, 250 * MM))
    draw_tree(ax, tree, xpos, ypos, lw=0.6, support=True)
    mx = max(xpos.values()); sp = meta.set_index("accession")
    for t in tips:
        r = sp.loc[t.name] if t.name in sp.index else None
        lab = f"{r.species} ({t.name})" if r is not None else t.name
        tr = tribe5(r.tribe) if r is not None else "Other Myrtaceae"
        ax.text(xpos[t] + mx * 0.01, ypos[t], lab, fontsize=4.6, va="center", style="italic",
                color=C["red"] if t.name == target else C["ink"], fontweight="bold" if t.name == target else "normal")
        ax.add_patch(Rectangle((mx * 1.46, ypos[t] - 0.45), mx * 0.02, 0.9, color=TRIBE_COL[tr], lw=0))
    ax.set_ylim(len(tips), -1); ax.set_xlim(-mx * 0.02, mx * 1.55)
    ax.plot([0, 0.005], [len(tips) - 0.3] * 2, color=C["ink"], lw=1.0)
    ax.text(0.0025, len(tips) - 0.9, "0.005 subst./site", ha="center", fontsize=7); ax.axis("off")
    ax.legend(handles=[Patch(color=TRIBE_COL[t], label=t) for t in TRIBE_ORDER], loc="lower left", frameon=False, fontsize=8)
    save(fig, out, "FigS4")


def arc_diagram(ax, seq, struct, title):
    n = len(seq); stack = []
    for i, c in enumerate(struct):
        if c == "(":
            stack.append(i)
        elif c == ")" and stack:
            j = stack.pop(); r = (i - j) / 2
            th = np.linspace(0, np.pi, 30)
            ax.plot(j + r + r * np.cos(th), r * np.sin(th) * 0.6, color=C["blue"], lw=0.6)
    for i, b in enumerate(seq):
        ax.text(i, -1.2, b, fontsize=4.5, ha="center", va="top", family=MONO)
    ax.set_xlim(-1, n); ax.set_ylim(-4, n * 0.32); ax.axis("off"); ax.set_title(title, fontsize=8, style="italic", loc="left")


def figS5_structures(out, res):
    sl = pd.read_csv(res / "tables" / "S12_3prime_stemloops.tsv", sep="\t").sort_values("z").head(8)
    fig, axs = plt.subplots(4, 2, figsize=(W2, 160 * MM), gridspec_kw=dict(hspace=0.55, wspace=0.1))
    for ax, r in zip(axs.flat, sl.itertuples()):
        arc_diagram(ax, r.stemloop_seq, r.structure, f"{r.gene}  (MFE {r.mfe} kcal/mol; z = {r.z})")
    save(fig, out, "FigS5")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--genomes-dir", default="genomes")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--target", default="PX251231.1")
    a = ap.parse_args()
    res, gdir, out = Path(a.results), Path(a.genomes_dir), Path(a.out)
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    fig_workflow(out)
    fig_map(out, res, gdir, a.target)
    fig_homology(out, res, gdir, gdir / "ref", a.target)
    fig_promoters(out, res)
    fig_sd(out, res)
    fig_stemloops(out, res)
    fig_conservation(out, res, a.target)
    figS1_spacers(out, res)
    figS2_tobacco(out, res)
    figS3_sd_conservation(out, res, a.target)
    figS4_tree(out, res, a.target)
    figS5_structures(out, res)


if __name__ == "__main__":
    main()
