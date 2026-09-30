"""Shared paths and helpers for the yeast membrane-protein structure pipeline."""
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("MP_DATA", ROOT / "data"))
IDMAPPING = DATA / "idmapping_2026_09_30.txt"
STRUCT_DIR = DATA / "structures"        # AlphaFold models, <ACC>.pdb
UNIPROT_DIR = DATA / "uniprot"          # UniProt JSON entries, <ACC>.json
RESULTS = Path(os.environ.get("MP_RESULTS", ROOT / "results"))
FIGURES = RESULTS / "figures"


def load_proteins(path=IDMAPPING):
    """Return one row per unique UniProt entry.

    The ID-mapping table can contain the same entry more than once (e.g. the
    query 'ctr1' also hits HNM1, whose gene synonyms include CTR1), so rows are
    de-duplicated on the UniProt accession, keeping the first occurrence.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["Length"] = df["Length"].astype(int)
    df = df.drop_duplicates("Entry").reset_index(drop=True)
    df["gene"] = df["Gene Names"].str.split().str[0]
    df["pfam"] = df["Pfam"].fillna("").str.strip(";").str.split(";").str[0]
    return df


# Readable names for the Pfam families in this set (first Pfam per entry).
PFAM_NAMES = {
    "PF00083": "Sugar_tr (MFS)",
    "PF07690": "MFS_1",
    "PF00324": "AA_permease (APC)",
    "PF06687": "SUR7",
    "PF01184": "Gpr1/Fun34/YaaH",
    "PF13520": "AA_permease_2 (APC)",
    "PF16944": "KCH",
    "PF04479": "RTA1",
}

# Categorical palette (fixed order, colour-blind validated) + neutral "Other".
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER_COLOR = "#a3a29c"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"


def family_label(pfam):
    return PFAM_NAMES.get(pfam, "Other")


def family_colors():
    """Map family label -> colour, in PFAM_NAMES order (largest family first)."""
    colors = {name: CATEGORICAL[i] for i, name in enumerate(PFAM_NAMES.values())}
    colors["Other"] = OTHER_COLOR
    return colors


def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(TEXT_2)
    ax.tick_params(colors=TEXT_2, labelcolor=TEXT_2)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
