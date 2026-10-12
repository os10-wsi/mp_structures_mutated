"""Structural summary: transmembrane topologies and length distribution.

Inputs : data/idmapping_*.txt, data/uniprot/<ACC>.json (from 01_fetch_data.py)
Outputs: results/topology_summary.tsv
         results/topology_classes.tsv
         results/figures/length_histogram.{png,pdf}
         results/figures/tm_count.{png,pdf}
         results/figures/topology_classes.{png,pdf}
         results/figures/topology_map.{png,pdf}

Topology is taken from the UniProt 'Transmembrane', 'Intramembrane' and
'Topological domain' features. A protein's topology class is its number of TM
helices plus the side of its N- and C-termini, e.g. "12 TM | N-in / C-in",
where "in" = cytoplasmic and "out" = extracellular/lumenal/vacuolar etc.
If UniProt only annotates one terminus, the other is inferred from the parity
of the TM count; if neither is annotated the side is shown as "?".
Length histograms only need the ID-mapping table and are always produced.
"""
import json
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from common import (FAMILY_BLUE_GREY, FIGURES, GREY, INK, RESULTS, UNIPROT_DIR,
                    family_label, load_proteins, paper_style)

paper_style()


def side(description):
    d = (description or "").lower()
    if "cytoplasm" in d:
        return "in"
    if d:
        return "out"
    return None


def parse_uniprot(acc):
    path = UNIPROT_DIR / f"{acc}.json"
    if not path.exists():
        return None
    entry = json.loads(path.read_text())
    tms, intra, domains = [], [], []
    for f in entry.get("features", []):
        loc = f["location"]
        start, end = loc["start"].get("value"), loc["end"].get("value")
        if start is None or end is None:
            continue
        if f["type"] == "Transmembrane":
            tms.append((start, end))
        elif f["type"] == "Intramembrane":
            intra.append((start, end))
        elif f["type"] == "Topological domain":
            domains.append((start, end, side(f.get("description"))))
    tms.sort()
    domains.sort()
    n_side = domains[0][2] if domains and domains[0][0] == 1 else None
    c_side = None
    if domains and tms and domains[-1][0] > tms[-1][1]:
        c_side = domains[-1][2]
    flip = {"in": "out", "out": "in"}
    n_tm = len(tms)
    if n_side and not c_side and n_tm:
        c_side = n_side if n_tm % 2 == 0 else flip[n_side]
    if c_side and not n_side and n_tm:
        n_side = c_side if n_tm % 2 == 0 else flip[c_side]
    return {"n_tm": n_tm, "n_intramembrane": len(intra),
            "n_term": n_side or "?", "c_term": c_side or "?",
            "tm_segments": ";".join(f"{s}-{e}" for s, e in tms),
            "intramembrane_segments": ";".join(f"{s}-{e}" for s, e in intra),
            "topo_domains": ";".join(f"{s}-{e}:{d}" for s, e, d in domains)}


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def length_histogram(df):
    bins = np.arange(250, 901, 50)
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.hist(df.Length, bins=bins, color=GREY, edgecolor="white", linewidth=0.6)
    med = df.Length.median()
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.15)
    ax.axvline(med, ymax=1 / 1.15, color=INK, lw=0.6, ls="--")
    ax.text(med, top * 1.02, f"median {med:.0f} aa", ha="center",
            va="bottom", fontsize=6.5)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Number of proteins")
    ax.set_xticks(bins[::2])
    save(fig, "length_histogram")


def tm_count_plot(df):
    counts = df.n_tm.value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.bar(counts.index, counts.values, color=GREY, width=0.75)
    ax.set_xticks(range(1, int(counts.index.max()) + 1))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Number of transmembrane helices")
    ax.set_ylabel("Number of proteins")
    save(fig, "tm_count")


def topology_class_plot(classes):
    fig, ax = plt.subplots(figsize=(3.6, 0.17 * len(classes) + 0.6))
    y = np.arange(len(classes))[::-1]
    ax.barh(y, classes["n_proteins"], color=GREY, height=0.7)
    for yi, (n, genes) in zip(y, zip(classes["n_proteins"], classes["genes"])):
        ax.text(n + 0.4, yi, textwrap.fill(genes, 60), va="center",
                fontsize=4.5 if len(genes) > 60 else 5.5, color="#555555",
                linespacing=1.0)
    ax.set_yticks(y, classes["topology"].str.replace(" | ", ", ",
                                                     regex=False))
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Number of proteins")
    ax.set_xlim(0, classes.n_proteins.max() * 1.05)
    save(fig, "topology_classes")


def topology_map(df):
    """One track per protein: sequence line with TM helices as boxes,
    coloured by Pfam family in shades of blue and grey."""
    order = {f: i for i, f in enumerate(FAMILY_BLUE_GREY)}
    d = (df.assign(fo=df.family.map(order))
         .sort_values(["fo", "n_tm", "Length"]).reset_index(drop=True))
    fig, ax = plt.subplots(figsize=(4.5, 0.11 * len(d) + 0.5))
    for i, r in d.iterrows():
        y = len(d) - i
        ax.plot([1, r.Length], [y, y], color="#ececec", lw=0.6,
                solid_capstyle="butt", zorder=1)
        for seg in filter(None, r.tm_segments.split(";")):
            s, e = map(int, seg.split("-"))
            ax.add_patch(plt.Rectangle((s, y - 0.32), e - s + 1, 0.64,
                                       color=FAMILY_BLUE_GREY[r.family],
                                       lw=0, zorder=2))
    ax.set_yticks(range(len(d), 0, -1), d.gene, fontsize=5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, d.Length.max() + 10)
    ax.set_ylim(0.3, len(d) + 0.7)
    ax.set_xlabel("Residue")
    ax.spines["left"].set_visible(False)
    present = [f for f in FAMILY_BLUE_GREY if f in set(d.family)]
    handles = [plt.Rectangle((0, 0), 1, 1, color=FAMILY_BLUE_GREY[f])
               for f in present]
    ax.legend(handles, present, fontsize=5.5, frameon=False,
              handlelength=1.2, handleheight=0.8, loc="upper left",
              bbox_to_anchor=(1.0, 1.0))
    save(fig, "topology_map")


def main():
    RESULTS.mkdir(exist_ok=True)
    df = load_proteins()
    df["family"] = df.pfam.map(family_label)
    length_histogram(df)

    topo = {acc: parse_uniprot(acc) for acc in df.Entry}
    have = [a for a, t in topo.items() if t is not None]
    print(f"UniProt topology available for {len(have)}/{len(df)} proteins")
    if not have:
        print("Run 01_fetch_data.py first to get topology; "
              "only the length histogram was made.")
        return
    df = df[df.Entry.isin(have)].copy()
    df = df.join(pd.DataFrame([topo[a] for a in df.Entry], index=df.index))
    df["topology"] = (df.n_tm.astype(str) + " TM | N-" + df.n_term
                      + " / C-" + df.c_term)
    cols = ["Entry", "gene", "Entry Name", "Protein names", "Length", "pfam",
            "family", "n_tm", "n_intramembrane", "n_term", "c_term",
            "topology", "tm_segments", "intramembrane_segments",
            "topo_domains"]
    df[cols].to_csv(RESULTS / "topology_summary.tsv", sep="\t", index=False)

    classes = (df.sort_values("gene").groupby("topology")
               .agg(n_proteins=("gene", "size"),
                    n_tm=("n_tm", "first"),
                    genes=("gene", lambda g: ", ".join(g)))
               .reset_index()
               .sort_values(["n_proteins", "n_tm"], ascending=[False, False]))
    classes.to_csv(RESULTS / "topology_classes.tsv", sep="\t", index=False)
    print(f"{len(classes)} distinct topology classes "
          f"({df.n_tm.nunique()} distinct TM counts)")
    print(classes[["topology", "n_proteins"]].to_string(index=False))

    tm_count_plot(df)
    topology_class_plot(classes)
    topology_map(df)


if __name__ == "__main__":
    main()
