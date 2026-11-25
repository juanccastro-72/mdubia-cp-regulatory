#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cpreg_circos_detail.py - detailed focal-genome Circos of the Myrciaria dubia
chloroplast genome (Fig. 10), in the style of a comparative-genomics Circos,
foregrounding the identified regulatory elements.

The M. dubia plastome occupies ~45 % of the circle (drawn to scale); 22
comparison plastomes (representatives of the five Myrtaceae tribes plus tobacco
and Arabidopsis) share the remaining arc on a common compressed scale. Tracks on
the M. dubia arc, outer -> inner:
  1. genes, split by strand, coloured by functional category, with gene labels;
  2. regulatory-element tracks (PEP, NEP, SD, 3' stem-loops, validated-element
     homologues) with the confidence tier encoded by bar height/opacity;
  3. five per-tribe conservation heat rings (fraction of each tribe's plastomes
     syntenic per 1-kb window);
  4. GC content (vs mean) and dispersed/IR repeats.
Synteny ribbons link M. dubia to each comparison plastome, coloured by tribe,
with opacity scaled to BLASTn identity.

Reuses tables written by cpreg_pipeline.py / cpreg_synteny.py (S1, S3, S6, S23,
S27) plus circos_work/{repeats,gc}.tsv. Fonts: Times New Roman, italic names.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from Bio import SeqIO
from cpreg_style import C, MM, apply_style, it, tribe5

apply_style()

GENE_CAT = {  # functional category -> colour (chloroplast)
    "psii_psi_pet": ("Photosystems / cyt b6f (psa, psb, pet)", C["green"]),
    "atp": ("ATP synthase (atp)", C["violet"]),
    "ndh": ("NADH dehydrogenase (ndh)", C["magenta"]),
    "rbcl": ("Rubisco (rbcL)", C["aqua"]),
    "rpo": ("RNA polymerase (rpo)", C["red"]),
    "ribosomal": ("Ribosomal proteins (rps, rpl)", C["yellow"]),
    "other": ("Other genes (matK, ycf, …)", C["ink2"]),
    "rrna": ("rRNA", "#8c4a9c"),
    "trna": ("tRNA", C["muted"]),
}
TRIBE_RAMP = {  # per-tribe heat ramps (light -> saturated)
    "Myrteae": ["#eef6ff", C["blue"]],
    "Syzygieae": ["#fdeee6", C["orange"]],
    "Eucalypteae": ["#e7f6f0", C["aqua"]],
    "Melaleuceae": ["#fdf3dd", C["yellow"]],
    "Other Myrtaceae": ["#fdeef4", C["magenta"]],
}
TRIBE_COL = {"Myrteae": C["blue"], "Syzygieae": C["orange"], "Eucalypteae": C["aqua"],
             "Melaleuceae": C["yellow"], "Other Myrtaceae": C["magenta"], "reference": C["ink2"]}
WIN = 1000


def gclass(g, t):
    if t == "tRNA": return "trna"
    if t == "rRNA": return "rrna"
    g = g.lower()
    if g.startswith(("psa", "psb", "pet")): return "psii_psi_pet"
    if g.startswith("rbc"): return "rbcl"
    if g.startswith("atp"): return "atp"
    if g.startswith("ndh"): return "ndh"
    if g.startswith("rpo"): return "rpo"
    if g.startswith(("rps", "rpl")): return "ribosomal"
    return "other"


def ramp(name):
    return LinearSegmentedColormap.from_list(name, TRIBE_RAMP[name])


def build(res: Path, gdir: Path, out: Path, comp, group_colors, ribbon_title, subtitle, name,
          md_frac=0.45):
    """comp: list of (accession, group_label); group_colors: {group_label: hex}."""
    from pycirclize import Circos
    refdir = gdir / "ref"
    meta = pd.read_csv(res / "tables" / "S1_genomes.tsv", sep="\t").set_index("accession")
    rec = SeqIO.read(str(gdir / "PX251231.1.gb"), "genbank"); L = len(rec); seq = str(rec.seq)
    integ = pd.read_csv(res / "tables" / "S23_integrated_evidence_target.tsv", sep="\t")
    loci = pd.read_csv(res / "tables" / "S3_target_loci.tsv", sep="\t").set_index("gene")
    hom = pd.read_csv(res / "tables" / "S6_homology_validated_elements.tsv", sep="\t")
    blocks = pd.read_csv(res / "tables" / "S27_synteny_blocks_all.tsv", sep="\t")
    reps = pd.read_csv("circos_work/repeats.tsv", sep="\t")
    gc = pd.read_csv("circos_work/gc.tsv", sep="\t")
    comp_accs = [a for a, _ in comp]
    GCOL = group_colors
    comp_len = {a: int(meta.loc[a, "length_bp"]) for a in comp_accs}

    # ---- sector sizes: M. dubia focal, others compressed ----
    tot_comp = sum(comp_len.values())
    unit = L / md_frac * (1 - md_frac)      # total "span" for all comp genomes
    comp_span = {a: comp_len[a] / tot_comp * unit for a in comp_accs}
    sectors = {"Mdubia": L}
    sectors.update({a: comp_span[a] for a in comp_accs})
    circos = Circos(sectors, space=1.3, start=2, end=358)

    md = circos.get_sector("Mdubia")

    # ---------- M. dubia: gene ring (split by strand) ----------
    gtrack = md.add_track((93, 100))
    gtrack.axis(fc="#faf9f5", ec=C["axis"], lw=0.3)
    label_pos, label_txt = [], []
    seen = set()
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA"):
            continue
        g = f.qualifiers.get("gene", [""])[0]
        cat = gclass(g, f.type)
        for p in f.location.parts:
            rlim = (96.6, 100) if p.strand == 1 else (93, 96.4)
            gtrack.rect(int(p.start) + 1, int(p.end), r_lim=rlim, fc=GENE_CAT[cat][1], ec="white", lw=0.1)
        if f.type in ("CDS", "rRNA") and g and g not in seen:
            seen.add(g)
            label_pos.append((int(f.location.start) + int(f.location.end)) / 2)
            label_txt.append(g)
    # gene labels are drawn after plotfig() with de-collision (see _gene_labels)
    # ---------- regulatory-element tracks ----------
    reg_defs = [("promoter_tier", "PEP promoter", C["red"], "five_prime_coord"),
                ("NEP", "NEP promoter", C["orange"], "five_prime_coord"),
                ("SD_tier", "Shine–Dalgarno", C["violet"], "five_prime_coord"),
                ("stemloop_tier", "3' stem-loop", C["green"], "three_prime_coord")]
    r_hi = 91.5
    for col, lab, colr, coord in reg_defs:
        tr = md.add_track((r_hi - 2.6, r_hi)); r_hi -= 3.1
        tr.axis(fc="#f7f6f1", ec="none")
        if col == "NEP":
            nep = pd.read_csv(res / "tables" / "S9_NEP_motifs.tsv", sep="\t")
            nep = nep.merge(loci.reset_index()[["gene", "five_prime_coord"]], left_on="id", right_on="gene", how="left")
            vmax = nep.YRTA_per_100nt.quantile(0.95)
            for x in nep.itertuples():
                if pd.isna(x.five_prime_coord):
                    continue
                v = min(x.YRTA_per_100nt / vmax, 1)
                if v <= 0.02:
                    continue
                a0 = max(int(x.five_prime_coord) - 300, 1); b0 = min(int(x.five_prime_coord), L)
                if a0 < b0:
                    tr.rect(a0, b0, r_lim=(r_hi + 3.1 - 2.6, r_hi + 3.1 - 2.6 + 2.6 * v), fc=colr, ec="none", alpha=0.85)
            continue
        for x in integ.itertuples():
            t = getattr(x, col, np.nan)
            if pd.isna(t) or x.gene not in loci.index:
                continue
            pos = loci.loc[x.gene, coord]
            pos = int(pos.iloc[0]) if hasattr(pos, "iloc") else int(pos)
            h = {1: 2.6, 2: 1.7, 3: 0.9}[int(t)]; al = {1: 1.0, 2: 0.8, 3: 0.45}[int(t)]
            base = r_hi + 3.1 - 2.6
            tr.rect(max(pos - 400, 1), min(pos + 400, L), r_lim=(base, base + h), fc=colr, ec="none", alpha=al)
    # validated-element homologues ring
    trh = md.add_track((r_hi - 2.6, r_hi)); r_hi -= 3.4
    trh.axis(fc="#f7f6f1", ec="none")
    for x in hom[hom.homology_supported & hom.Mdubia_ortholog_coord.notna()].itertuples():
        a, b = [int(v) for v in x.Mdubia_ortholog_coord.split()[0].split("-")]
        trh.rect(a, b, r_lim=(r_hi + 3.4 - 2.6, r_hi + 3.4), fc=C["blue"], ec="none", alpha=0.9)

    # ---------- per-tribe conservation heat rings ----------
    nwin = int(np.ceil(L / WIN))
    tribes = ["Myrteae", "Syzygieae", "Eucalypteae", "Melaleuceae", "Other Myrtaceae"]
    status = {}  # acc -> set of covered windows
    for a, g in blocks.groupby("accession"):
        wins = set()
        for r in g.itertuples():
            for w in range((r.ref_start - 1) // WIN, min((r.ref_end - 1) // WIN + 1, nwin)):
                wins.add(w)
        status[a] = wins
    myr = [a for a in status if a in meta.index and meta.loc[a, "role"] != "reference (validated elements)"]
    r_hi -= 0.5
    for trb in tribes:
        accs = [a for a in myr if tribe5(meta.loc[a, "tribe"]) == trb]
        if not accs:
            continue
        frac = np.array([sum(w in status[a] for a in accs) / len(accs) for w in range(nwin)])
        tr = md.add_track((r_hi - 2.3, r_hi)); r_hi -= 2.6
        cmap = ramp(trb)
        for w, v in enumerate(frac):
            a0 = w * WIN + 1; b0 = min((w + 1) * WIN, L)
            if a0 < b0:
                tr.rect(a0, b0, fc=cmap(0.1 + 0.9 * v), ec="none")

    # ---------- GC content ----------
    r_hi -= 0.5
    gtr = md.add_track((r_hi - 6, r_hi)); r_hi -= 7
    meang = float((gc.gc).mean())
    x = np.minimum(gc.start.values + WIN / 2, L - 1.0)
    gtr.axis(fc="white", ec=C["axis"], lw=0.3)
    gtr.fill_between(x, gc.gc.values, meang, color=C["muted"], alpha=0.8, lw=0)
    # ---------- repeats ring (inverted repeats as links at the very inside) ----------
    rtr = md.add_track((r_hi - 1.6, r_hi)); rtr.axis(fc="#f0efe9", ec="none"); r_inner = r_hi - 2
    for rp in reps.itertuples():
        col = C["red"] if rp.kind == "inverted" else C["ink2"]
        rtr.rect(int(rp.a), int(rp.b), fc=col, ec="none", alpha=0.6)

    # ---------- comparison genomes: arcs + thin gene ticks ----------
    for a, trb in comp:
        sec = circos.get_sector(a)
        col = GCOL[trb]
        band = sec.add_track((97, 100)); band.axis(fc=col, ec="white", lw=0.2)
        # thin gene ticks from annotation
        folder = refdir if trb == "reference" else gdir
        try:
            crec = SeqIO.read(str(folder / f"{a}.gb"), "genbank"); cl = len(crec)
            gt = sec.add_track((93.5, 96.5)); gt.axis(fc="#faf9f5", ec="none")
            scale = comp_span[a] / cl
            for f in crec.features:
                if f.type not in ("CDS", "rRNA", "tRNA"):
                    continue
                g = f.qualifiers.get("gene", [""])[0]
                s0 = int(f.location.start) * scale; e0 = int(f.location.end) * scale
                if e0 - s0 < 1:
                    e0 = s0 + comp_span[a] * 0.001
                gt.rect(min(s0, comp_span[a]), min(e0, comp_span[a]), fc=GENE_CAT[gclass(g, f.type)][1], ec="none")
        except Exception:
            pass

    # ---------- synteny ribbons ----------
    ident_norm = Normalize(85, 100)
    for a, trb in comp:
        col = GCOL[trb]
        sub = blocks[(blocks.accession == a) & (blocks.length >= 1000)]
        scale = comp_span[a] / comp_len[a]
        for b in sub.itertuples():
            alpha = 0.15 + 0.55 * float(ident_norm(min(max(b.identity, 85), 100)))
            q0, q1 = b.qry_start * scale, b.qry_end * scale
            if b.strand == -1:
                q0, q1 = q1, q0
            circos.link(("Mdubia", b.ref_start, b.ref_end), (a, q0, q1),
                        color=col, alpha=alpha, r1=r_inner, r2=96.8, lw=0.05)

    fig = circos.plotfig(dpi=100)
    fig.set_size_inches(205 * MM, 250 * MM)
    ax = fig.axes[0]
    ax.set_position([0.06, 0.27, 0.88, 0.66])
    off = ax.get_theta_offset(); direction = ax.get_theta_direction()

    def screen_deg(rad):
        return np.degrees(off + direction * rad)

    def radial_label(rad, r, text, color, size, bold=False, italic=True):
        sd = screen_deg(rad) % 360
        if 90 < sd < 270:
            rot = sd + 180; ha = "right"
        else:
            rot = sd; ha = "left"
        ax.text(rad, r, text, rotation=rot, rotation_mode="anchor", ha=ha, va="center",
                fontsize=size, color=color, style="italic" if italic else "normal",
                fontweight="bold" if bold else "normal")

    # ---- M. dubia gene labels, greedy radial de-collision ----
    # each label needs ~ANG_MIN degrees of tangential room at its radius; when the
    # local angular slot is taken, push the label outward by RSTEP and draw a leader.
    items = sorted(zip(label_pos, label_txt), key=lambda z: md.x_to_rad(z[0]))
    ANG_MIN, RSTEP, R0 = 2.6, 3.0, 101.5
    placed = []  # (screen_deg, radius)
    rmax_used = R0
    for pos, g in items:
        rad = md.x_to_rad(pos); sd = screen_deg(rad) % 360
        r = R0
        while any(abs(((sd - psd + 180) % 360) - 180) < ANG_MIN and abs(pr - r) < RSTEP - 0.3
                  for psd, pr in placed):
            r += RSTEP
        placed.append((sd, r)); rmax_used = max(rmax_used, r)
        ax.plot([rad, rad], [100.2, r - 0.5], color=C["muted"], lw=0.3, zorder=1)
        radial_label(rad, r, g, C["ink"], 5.6)
    ax.set_rmax(max(118, rmax_used + 10))

    # ---- comparison genome labels (species + size + % covered) ----
    su = pd.read_csv(res / "tables" / "S28_synteny_per_genome_summary.tsv", sep="\t")
    su = su.set_index("accession")
    for a, trb in comp:
        sec = circos.get_sector(a)
        mid = sec.x_to_rad(comp_span[a] / 2)
        sp = meta.loc[a, "species"]
        parts = sp.split()
        abbr = f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else sp
        cov = int(round(su.loc[a, "collinear_fraction"] * 100)) if a in su.index else None
        size_kb = comp_len[a] / 1000
        stat = f"  ·  {size_kb:.0f} kb · {cov}%" if cov is not None else f"  ·  {size_kb:.0f} kb"
        # one italic species label + a roman stat segment, on a single radial line
        sd = screen_deg(mid) % 360
        flip = 90 < sd < 270
        ha = "right" if flip else "left"
        rot = sd + 180 if flip else sd
        ax.text(mid, 101.5, (stat + "  " + abbr) if flip else (abbr + stat), rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=6.2, color=GCOL[trb], style="italic")

    fig.text(0.5, 0.99, it("Myrciaria dubia") + " chloroplast genome — " + ribbon_title, ha="center", va="top", size=14, color=C["red"], fontweight="bold")
    fig.text(0.5, 0.972, f"{L/1000:.1f} kb  ·  GC 37.1 %  ·  78 protein-coding · 4 rRNA · 29 tRNA genes (IR copies collapsed)  ·  {len(comp)} comparison plastomes",
             ha="center", va="top", size=9.5, color=C["ink2"])
    _legend(fig, GCOL, ribbon_title)
    _save(fig, out, name)


def _legend(fig, group_colors, ribbon_title):
    ax = fig.add_axes([0, 0, 1, 0.20]); ax.axis("off")
    rw, dy = 0.013, 0.085
    # col 1: gene categories
    x0 = 0.015
    ax.text(x0, 0.97, "Gene categories", fontsize=9, fontweight="bold", va="top")
    ax.text(x0, 0.90, "outer ring: + strand / − strand", fontsize=6.5, color=C["ink2"], va="top")
    y = 0.80
    for key, (lab, col) in GENE_CAT.items():
        ax.add_patch(plt.Rectangle((x0, y - 0.028), rw, 0.045, color=col, transform=ax.transAxes))
        ax.text(x0 + 0.018, y, lab, fontsize=7, va="center"); y -= dy
    # col 2: regulatory elements
    x1 = 0.30
    ax.text(x1, 0.97, "Regulatory-element rings", fontsize=9, fontweight="bold", va="top")
    ax.text(x1, 0.90, "outer → inner; bar height/opacity = Tier 1 > 2 > 3", fontsize=6.5, color=C["ink2"], va="top")
    reg = [("PEP promoter", C["red"]), ("NEP promoter (YRTA density)", C["orange"]),
           ("Shine–Dalgarno motif", C["violet"]), ("3' stem-loop", C["green"]),
           ("Homolog of validated element", C["blue"])]
    y = 0.80
    for lab, col in reg:
        ax.add_patch(plt.Rectangle((x1, y - 0.028), rw, 0.045, color=col, transform=ax.transAxes))
        ax.text(x1 + 0.018, y, lab, fontsize=7, va="center"); y -= dy
    # col 3: tribe conservation
    x2 = 0.56
    ax.text(x2, 0.97, "Per-tribe conservation rings", fontsize=9, fontweight="bold", va="top")
    ax.text(x2, 0.90, "fraction of tribe syntenic per 1 kb (0→1)", fontsize=6.5, color=C["ink2"], va="top")
    y = 0.80
    for trb in ["Myrteae", "Syzygieae", "Eucalypteae", "Melaleuceae", "Other Myrtaceae"]:
        cmap = ramp(trb)
        for i in range(12):
            ax.add_patch(plt.Rectangle((x2 + i * 0.0055, y - 0.028), 0.0055, 0.045, color=cmap(i / 11), transform=ax.transAxes))
        ax.text(x2 + 0.075, y, trb, fontsize=7, va="center"); y -= dy
    # col 4: ribbons (group colours supplied by caller)
    x3 = 0.80
    ax.text(x3, 0.97, f"Synteny ribbons ({ribbon_title})", fontsize=9, fontweight="bold", va="top")
    ax.text(x3, 0.90, "opacity ∝ BLASTn identity (≥ 85 %)", fontsize=6.5, color=C["ink2"], va="top")
    y = 0.80
    step = 0.072 if len(group_colors) <= 6 else 0.052
    for grp, col in group_colors.items():
        ax.add_patch(plt.Rectangle((x3, y - 0.024), rw, 0.04, color=col, transform=ax.transAxes))
        ax.text(x3 + 0.018, y, grp, fontsize=6.8, va="center"); y -= step
    ax.text(x3, y + 0.005, "Innermost M. dubia rings: GC (grey, vs mean)\nand repeats (red = inverted / IR).",
            fontsize=6.3, color=C["ink2"], va="top", linespacing=1.4)


def _save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    svg = out / f"{name}.svg"; fig.savefig(svg, format="svg", bbox_inches="tight")
    svg.write_text(svg.read_text().replace("Liberation Serif", "Times New Roman"))
    fig.savefig(out / f"{name}.jpg", dpi=600, bbox_inches="tight", pil_kwargs={"quality": 95})
    plt.close(fig)


COMP22 = [
    ("KU318111.1", "Myrteae"), ("KX364403.1", "Myrteae"), ("ON641286.1", "Myrteae"),
    ("MW353255.1", "Myrteae"), ("OR663915.1", "Myrteae"), ("MK044696.1", "Myrteae"),
    ("MT975437.1", "Syzygieae"), ("MZ929420.1", "Syzygieae"), ("OP778218.1", "Syzygieae"), ("OR712057.1", "Syzygieae"),
    ("KC180807.1", "Eucalypteae"), ("MT948196.1", "Eucalypteae"), ("MW446439.1", "Eucalypteae"),
    ("KC180794.1", "Eucalypteae"), ("MW446421.1", "Eucalypteae"),
    ("MT700493.1", "Melaleuceae"), ("MT731621.1", "Melaleuceae"),
    ("MW837774.1", "Other Myrtaceae"), ("ON422330.1", "Other Myrtaceae"), ("MK726014.1", "Other Myrtaceae"),
    ("NC_001879.2", "reference"), ("NC_000932.1", "reference"),
]


# tribe colours (real tribes) and reference-group colours
TRIBE9 = {"Myrteae": C["blue"], "Syzygieae": C["orange"], "Eucalypteae": C["aqua"],
          "Melaleuceae": C["yellow"], "Backhousieae": C["green"], "Chamelaucieae": C["violet"],
          "Lophostemoneae": "#b5651d", "Xanthostemoneae": C["magenta"],
          "Heteropyxideae (Psiloxyloideae)": C["ink2"]}
REFGRP = {"Eudicot": C["blue"], "Monocot": C["orange"], "Bryophyte": C["green"], "Green alga": C["magenta"]}

TRIBE_REPS = [  # (accession, tribe) - >=2 per tribe where available
    ("KU318111.1", "Myrteae"), ("KX364403.1", "Myrteae"),
    ("MT975437.1", "Syzygieae"), ("MN746306.1", "Syzygieae"),
    ("HM347959.1", "Eucalypteae"), ("KC180801.1", "Eucalypteae"),
    ("MT700493.1", "Melaleuceae"), ("MT731621.1", "Melaleuceae"),
    ("ON422330.1", "Backhousieae"), ("ON470204.1", "Chamelaucieae"),
    ("OM640421.1", "Lophostemoneae"), ("MW837774.1", "Xanthostemoneae"),
    ("MK726014.1", "Heteropyxideae (Psiloxyloideae)")]
REF12 = [("NC_001879.2", "Eudicot"), ("NC_008096.2", "Eudicot"), ("NC_007578.1", "Eudicot"),
         ("NC_043882.1", "Eudicot"), ("NC_000932.1", "Eudicot"), ("NC_003119.8", "Eudicot"),
         ("NC_009259.1", "Eudicot"), ("NC_008590.1", "Monocot"), ("NC_001666.2", "Monocot"),
         ("NC_001320.1", "Monocot"), ("NC_037507.1", "Bryophyte"), ("NC_005353.1", "Green alga")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--genomes-dir", default="genomes")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--which", default="all", choices=["all", "fig10", "fig11", "fig12"])
    a = ap.parse_args()
    res, gdir, out = Path(a.results), Path(a.genomes_dir), Path(a.out)
    comp10 = COMP22
    gcol10 = {"Myrteae": C["blue"], "Syzygieae": C["orange"], "Eucalypteae": C["aqua"],
              "Melaleuceae": C["yellow"], "Other Myrtaceae": C["magenta"], "reference": C["ink2"]}
    if a.which in ("all", "fig10"):
        build(res, gdir, out, comp10, gcol10, "22 Myrtaceae + references",
              f"{158.664:.1f} kb  ·  GC 37.1 %  ·  78 protein-coding · 4 rRNA · 29 tRNA genes (IR copies collapsed)  ·  22 comparison plastomes",
              "Fig10", md_frac=0.48)
        print("Fig10 done")
    if a.which in ("all", "fig11"):
        build(res, gdir, out, TRIBE_REPS, TRIBE9, "Myrtaceae tribes (≥2 per tribe)",
              "Synteny of the M. dubia plastome with two representatives of each Myrtaceae tribe (one where only one is available)",
              "Fig11", md_frac=0.5)
        print("Fig11 done")
    if a.which in ("all", "fig12"):
        build(res, gdir, out, REF12, REFGRP, "12 reference plastomes",
              "Synteny of the M. dubia plastome with the 12 reference plastomes carrying experimentally validated regulatory elements",
              "Fig12", md_frac=0.5)
        print("Fig12 done")


if __name__ == "__main__":
    main()
