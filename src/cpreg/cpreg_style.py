#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared plotting style for the camu-camu plastid-regulatory-element figures.

Times New Roman (rendered with the metric-compatible Liberation Serif and
relabelled to 'Times New Roman' in the SVG so it is correct in CorelDRAW),
10 pt base font, italics for scientific and gene names. SVG text stays live.
"""
from __future__ import annotations
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4

C = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100", magenta="#e87ba4",
         green="#008300", violet="#6250d6", red="#e34948", ink="#0b0b0b", ink2="#52514e",
         muted="#898781", grid="#e1e0d9", axis="#c3c2b7", neutral="#f0efec")
SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
TRIBE_ORDER = ["Myrteae", "Syzygieae", "Eucalypteae", "Melaleuceae", "Other Myrtaceae"]
TRIBE_COL = {"Myrteae": C["blue"], "Syzygieae": C["orange"], "Eucalypteae": C["aqua"],
             "Melaleuceae": C["yellow"], "Other Myrtaceae": C["muted"]}

# base size 10; a few sub-elements step down slightly to prevent overlap in dense panels
BASE = 10


def apply_style():
    plt.rcParams.update({
        "font.family": "Liberation Serif",          # metric-compatible with Times New Roman
        "mathtext.fontset": "stix",                  # serif math, italics for $\it{...}$
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "xtick.labelsize": BASE - 1, "ytick.labelsize": BASE - 1, "legend.fontsize": BASE - 1,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.edgecolor": C["ink2"], "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "figure.dpi": 100,
    })


def it(name: str) -> str:
    """Return a mathtext italic label for a scientific/gene name (spaces preserved)."""
    return "$\\it{" + name.replace(" ", "\\ ") + "}$"


def tribe5(t: str) -> str:
    return t if t in TRIBE_ORDER[:4] else "Other Myrtaceae"


def panel(ax, letter, x=-0.14, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=BASE + 2, fontweight="bold",
            va="bottom", ha="left")


def save(fig, out: Path, name: str, jpg=True):
    out.mkdir(parents=True, exist_ok=True)
    svg = out / f"{name}.svg"
    fig.savefig(svg, format="svg", bbox_inches="tight")
    txt = svg.read_text().replace("Liberation Serif", "Times New Roman")
    svg.write_text(txt)
    if jpg:
        fig.savefig(out / f"{name}.jpg", dpi=600, bbox_inches="tight", pil_kwargs={"quality": 95})
    plt.close(fig)
