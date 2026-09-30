"""Structural similarity network from Foldseek all-vs-all hit probabilities.

Nodes are proteins; edges are Foldseek hit probabilities (symmetrised as the
mean of the two search directions; self-hits removed). The script

  1. writes results/network/similarity_network.gexf for Gephi, with node
     attributes (gene, family, Pfam, length, TM count, topology), edge weight =
     probability, and pre-computed positions/colours so it opens ready-laid-out;
  2. reproduces the Gephi layout in Python — Fruchterman-Reingold to
     equilibrium, then ForceAtlas2 with node sizes to prevent overlap — and
     draws results/figures/similarity_network.{png,pdf};
  3. writes the probability matrix (results/network/prob_matrix.tsv) and the
     node / edge tables.

For the Gephi route from the paper: File > Open the .gexf, Layout >
Fruchterman Reingold (run until stable), then Layout > ForceAtlas 2 with
"Prevent Overlap" ticked.

Usage: python 04_similarity_network.py [--min-prob 0.5] [--seed 7]
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from common import (FIGURES, RESULTS, TEXT, TEXT_2, family_colors,
                    family_label, load_proteins)

FS_OUT = RESULTS / "foldseek" / "foldseek_easy_allvsall"
NET_DIR = RESULTS / "network"
COLS = ["query", "target", "alntmscore", "qtmscore", "ttmscore", "prob"]


def load_hits():
    hits = pd.read_csv(FS_OUT, sep="\t", names=COLS)
    for c in ("query", "target"):   # Foldseek may keep the file extension
        hits[c] = hits[c].str.replace(r"\.(pdb|cif)(\.gz)?$", "", regex=True)
    hits = hits[hits["query"] != hits["target"]]
    # best hit per ordered pair (multiple alignments can be reported)
    return hits.sort_values("prob").drop_duplicates(["query", "target"],
                                                     keep="last")


def symmetric(hits, col, ids):
    m = (hits.pivot(index="query", columns="target", values=col)
         .reindex(index=ids, columns=ids).fillna(0.0))
    a = (m.values + m.values.T) / 2
    np.fill_diagonal(a, 1.0)
    return pd.DataFrame(a, index=ids, columns=ids)


def node_table(ids):
    df = load_proteins().set_index("Entry")
    topo_path = RESULTS / "topology_summary.tsv"
    if topo_path.exists():
        topo = pd.read_csv(topo_path, sep="\t").set_index("Entry")
        df = df.join(topo[["n_tm", "topology"]])
    else:
        df["n_tm"], df["topology"] = -1, "unknown"
    df["family"] = df.pfam.map(family_label)
    return df.loc[ids]


def layout(G, sizes, seed):
    # Fruchterman-Reingold to (near) equilibrium
    pos = nx.spring_layout(G, weight="weight", iterations=2000,
                           threshold=1e-6, seed=seed)
    # ForceAtlas2, preventing overlap via node sizes
    pos = nx.forceatlas2_layout(G, pos=pos, weight="weight", max_iter=2000,
                                scaling_ratio=2.0, gravity=1.0,
                                node_size=sizes, seed=seed)
    return pos


def draw(G, pos, nodes, colors, min_prob):
    fig, ax = plt.subplots(figsize=(9, 8))
    edges = sorted(G.edges(data="weight"), key=lambda e: e[2])
    for u, v, w in edges:
        (x0, y0), (x1, y1) = pos[u], pos[v]
        ax.plot([x0, x1], [y0, y1], color="#8a8983",
                lw=0.2 + 1.2 * w, alpha=0.08 + 0.4 * w, zorder=1)
    sizes = 30 + 0.25 * nodes.Length
    xy = np.array([pos[n] for n in nodes.index])
    ax.scatter(xy[:, 0], xy[:, 1], s=sizes,
               c=[colors[f] for f in nodes.family],
               edgecolors="white", linewidths=1.2, zorder=2)
    span = np.ptp(xy, axis=0).max()
    for (x, y), gene, s in zip(xy, nodes.gene, sizes):
        ax.text(x, y + span * 0.012 + np.sqrt(s) * span / 900, gene,
                ha="center", va="bottom", fontsize=6.5, color=TEXT, zorder=3)
    ax.set_aspect("equal")
    ax.axis("off")
    present = [f for f in colors if f in set(nodes.family)]
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=colors[f])
               for f in present]
    leg = ax.legend(handles, present, title="Pfam family", frameon=False,
                    fontsize=7.5, title_fontsize=8,
                    bbox_to_anchor=(1.0, 1.0), loc="upper left")
    leg.get_title().set_color(TEXT)
    ax.set_title("Structural similarity network (Foldseek)", loc="left",
                 fontsize=12, fontweight="bold", color=TEXT, pad=16)
    ax.text(0, 1.0, f"edges: hit probability ≥ {min_prob:g} "
            f"(width ∝ probability); node size ∝ length",
            transform=ax.transAxes, fontsize=7.5, color=TEXT_2, va="bottom")
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"similarity_network.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-prob", type=float, default=0.5,
                    help="minimum symmetrised Foldseek probability for an edge")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    hits = load_hits()
    ids = sorted(set(hits["query"]) | set(hits["target"]))
    nodes = node_table(ids)
    prob = symmetric(hits, "prob", ids)
    tm = symmetric(hits, "alntmscore", ids)

    NET_DIR.mkdir(parents=True, exist_ok=True)
    labelled = prob.rename(index=nodes.gene, columns=nodes.gene)
    labelled.to_csv(NET_DIR / "prob_matrix.tsv", sep="\t", float_format="%.3f")

    colors = family_colors()
    G = nx.Graph()
    for acc, r in nodes.iterrows():
        G.add_node(acc, label=r.gene, family=r.family, pfam=r.pfam,
                   length=int(r.Length), n_tm=int(r.n_tm),
                   topology=str(r.topology))
    iu = np.triu_indices(len(ids), k=1)
    for i, j in zip(*iu):
        p = prob.values[i, j]
        if p >= args.min_prob:
            G.add_edge(ids[i], ids[j], weight=float(p),
                       alntmscore=float(tm.values[i, j]))
    print(f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
          f"(prob >= {args.min_prob}), "
          f"{nx.number_connected_components(G)} connected components")

    sizes = {n: 0.02 + 0.00003 * nodes.loc[n, "Length"] for n in G}
    pos = layout(G, sizes, args.seed)
    draw(G, pos, nodes, colors, args.min_prob)

    # Gephi export with positions and colours baked in
    for n in G:
        h = colors[G.nodes[n]["family"]].lstrip("#")
        r, g, b = (int(h[k:k + 2], 16) for k in (0, 2, 4))
        G.nodes[n]["viz"] = {
            "position": {"x": float(pos[n][0] * 1000),
                         "y": float(pos[n][1] * 1000), "z": 0.0},
            "color": {"r": r, "g": g, "b": b, "a": 1.0},
            "size": float(5 + nodes.loc[n, "Length"] / 50)}
    nx.write_gexf(G, NET_DIR / "similarity_network.gexf")
    nx.to_pandas_edgelist(G).rename(columns={"source": "Source",
                                             "target": "Target",
                                             "weight": "Weight"}) \
        .to_csv(NET_DIR / "edges.csv", index=False)
    nodes.assign(Id=nodes.index, Label=nodes.gene)[
        ["Id", "Label", "family", "pfam", "Length", "n_tm", "topology"]] \
        .to_csv(NET_DIR / "nodes.csv", index=False)
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    with open(NET_DIR / "clusters.txt", "w") as fh:
        for k, c in enumerate(comps, 1):
            fh.write(f"cluster {k} ({len(c)}): "
                     + ", ".join(sorted(nodes.loc[list(c), "gene"])) + "\n")


if __name__ == "__main__":
    main()
