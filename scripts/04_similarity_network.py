"""Structural similarity network from Foldseek all-vs-all hit probabilities.

Nodes are proteins; edges are Foldseek hit probabilities (symmetrised as the
mean of the two search directions; self-hits removed). The script

  1. writes results/network/similarity_network.gexf for Gephi, with node
     attributes (gene, family, Pfam, length, TM count, topology), edge weight =
     probability, and pre-computed positions/colours so it opens ready-laid-out;
  2. reproduces the Gephi layout in Python — Fruchterman-Reingold to
     equilibrium, then ForceAtlas2, then an overlap-removal pass — and
     draws results/figures/similarity_network{,_nolabels}.{png,pdf};
  3. writes the probability matrix (results/network/prob_matrix.tsv) and the
     node / edge tables.

For the Gephi route from the paper: File > Open the .gexf, Layout >
Fruchterman Reingold (run until stable), then Layout > ForceAtlas 2 with
"Prevent Overlap" ticked.

Usage: python 04_similarity_network.py [--min-prob 0.01] [--cluster-prob 0.5]
                                       [--seed 7]
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from common import FIGURES, RESULTS, family_label, load_proteins

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


def prevent_overlap(pos, radii, iters=500):
    """Push apart nodes closer than the sum of their radii (like Gephi's
    'Prevent Overlap'). Positions are first scaled to a unit extent."""
    keys = list(pos)
    xy = np.array([pos[k] for k in keys], dtype=float)
    xy = (xy - xy.mean(0)) / np.ptp(xy, axis=0).max()
    r = np.array([radii[k] for k in keys])
    rng = np.random.default_rng(0)
    for _ in range(iters):
        diff = xy[:, None] - xy[None]
        d = np.sqrt((diff ** 2).sum(-1))
        np.fill_diagonal(d, np.inf)
        overlap = (r[:, None] + r[None]) - d
        if (overlap <= 0).all():
            break
        d = np.where(d == 0, 1e-9, d)
        diff += (d == 1e-9)[..., None] * rng.normal(0, 1e-4, diff.shape)
        push = np.clip(overlap, 0, None)[..., None] * diff / d[..., None] / 2
        xy += push.sum(1) * 0.5
    return dict(zip(keys, xy))


def layout(G, radii, seed):
    # Layout weight = probability cubed, so the many weak (<0.1) hits only
    # nudge nodes while confident hits pull families into tight clusters.
    H = nx.Graph()
    H.add_nodes_from(G)
    H.add_weighted_edges_from((u, v, w ** 3)
                              for u, v, w in G.edges(data="weight"))
    # Fruchterman-Reingold to (near) equilibrium
    pos = nx.spring_layout(H, weight="weight", iterations=2000,
                           threshold=1e-6, seed=seed)
    # ForceAtlas2 from the FR positions; gravity keeps unconnected nodes near
    pos = nx.forceatlas2_layout(H, pos=pos, weight="weight", max_iter=3000,
                                scaling_ratio=2.0, gravity=10.0, seed=seed)
    # networkx's node_size option is numerically unstable here, so
    # overlap is removed as a separate step
    return prevent_overlap(pos, radii)


# Colour classes, in the spirit of the SCOP-class colouring in the Domainome
# network: two transporter folds by Pfam, the rest split by TM count.
FOLD_PFAMS = {
    "MFS fold": {"PF00083", "PF07690"},
    "LeuT fold (APC)": {"PF00324", "PF13520", "PF00474", "PF02133"},
}
CLASS_COLORS = {
    "MFS fold": "#d7262b",
    "LeuT fold (APC)": "#2c3a96",
    "Other, ≥8 TM": "#a3a6dc",
    "≤7 TM": "#efa63a",
}
SHORT_PFAM = {
    "PF00083": "Sugar_tr", "PF07690": "MFS_1", "PF00324": "AA_permease",
    "PF13520": "AA_permease_2", "PF06687": "SUR7", "PF01184": "Gpr1/Fun34",
    "PF16944": "KCH", "PF04479": "RTA1",
}


def fold_class(r):
    for name, pfams in FOLD_PFAMS.items():
        if r.pfam in pfams:
            return name
    return "Other, ≥8 TM" if r.n_tm >= 8 else "≤7 TM"


def draw(G, strong, pos, nodes, labels=True, suffix=""):
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica",
                                             "Liberation Sans", "DejaVu Sans"],
                         "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6, 6))
    col = {n: CLASS_COLORS[c] for n, c in nodes["fold_class"].items()}
    for u, v, w in sorted(G.edges(data="weight"), key=lambda e: e[2]):
        (x0, y0), (x1, y1) = pos[u], pos[v]
        ax.plot([x0, x1], [y0, y1], color=col[u], lw=0.15 + 0.6 * w,
                alpha=0.05 + 0.35 * w, zorder=1, solid_capstyle="round")
    xy = np.array([pos[n] for n in nodes.index])
    ax.scatter(xy[:, 0], xy[:, 1], s=26, c=[col[n] for n in nodes.index],
               linewidths=0, zorder=2)
    ax.set_aspect("equal")
    ax.axis("off")
    pad = 0.08
    ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    if labels:
        # one label per cluster (connected component of confident hits)
        centre = xy.mean(0)
        for comp in nx.connected_components(strong):
            if len(comp) < 2:
                continue
            grp = nodes.loc[sorted(comp)]
            top = grp.fold_class.value_counts().index[0]
            if top in FOLD_PFAMS:
                name = top
            else:
                name = "/".join(SHORT_PFAM[p] for p in
                                grp.pfam.value_counts().index
                                if p in SHORT_PFAM)
            # place by the cluster core so a weakly attached member
            # (e.g. BOR1) does not drag the label away
            pts = np.array([pos[n] for n in comp])
            c = np.median(pts, axis=0)
            out = c - centre
            out = out / (np.linalg.norm(out) or 1)
            radius = np.percentile(np.linalg.norm(pts - c, axis=1), 80)
            x, y = c + out * (radius + 0.06)
            ax.text(x, y, f"{name}\n({len(comp)})", ha="center",
                    va="center", fontsize=7, linespacing=1.1, zorder=3)
        counts = nodes.fold_class.value_counts()
        lines = [f"{len(nodes)} proteins"] + [
            f"{c} ({counts.get(c, 0)})" for c in CLASS_COLORS]
        cols = ["black"] + list(CLASS_COLORS.values())
        for k, (t, c) in enumerate(zip(lines, cols)):
            ax.text(0.0, 1.0 - k * 0.035, t, color=c, fontsize=7.5,
                    transform=ax.transAxes, va="top")
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"similarity_network{suffix}.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-prob", type=float, default=0.01,
                    help="minimum symmetrised Foldseek probability for an "
                         "edge in the network (the paper used all of them)")
    ap.add_argument("--cluster-prob", type=float, default=0.5,
                    help="probability threshold for reporting clusters")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    hits = load_hits()
    ids = sorted(set(hits["query"]) | set(hits["target"]))
    nodes = node_table(ids)
    nodes["fold_class"] = nodes.apply(fold_class, axis=1)
    prob = symmetric(hits, "prob", ids)
    tm = symmetric(hits, "alntmscore", ids)

    NET_DIR.mkdir(parents=True, exist_ok=True)
    labelled = prob.rename(index=nodes.gene, columns=nodes.gene)
    labelled.to_csv(NET_DIR / "prob_matrix.tsv", sep="\t", float_format="%.3f")

    G = nx.Graph()
    for acc, r in nodes.iterrows():
        G.add_node(acc, label=r.gene, family=r.family, pfam=r.pfam,
                   fold_class=r.fold_class, length=int(r.Length),
                   n_tm=int(r.n_tm), topology=str(r.topology))
    iu = np.triu_indices(len(ids), k=1)
    for i, j in zip(*iu):
        p = prob.values[i, j]
        if p >= args.min_prob:
            G.add_edge(ids[i], ids[j], weight=float(p),
                       alntmscore=float(tm.values[i, j]))
    strong = nx.Graph([(u, v) for u, v, w in G.edges(data="weight")
                       if w >= args.cluster_prob])
    strong.add_nodes_from(G)
    print(f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
          f"(prob >= {args.min_prob}); "
          f"{nx.number_connected_components(strong)} clusters at prob >= "
          f"{args.cluster_prob}")

    radii = {n: 0.011 for n in G}
    pos = layout(G, radii, args.seed)
    draw(G, strong, pos, nodes, labels=True)
    draw(G, strong, pos, nodes, labels=False, suffix="_nolabels")

    # Gephi export with positions and colours baked in
    for n in G:
        h = CLASS_COLORS[G.nodes[n]["fold_class"]].lstrip("#")
        r, g, b = (int(h[k:k + 2], 16) for k in (0, 2, 4))
        G.nodes[n]["viz"] = {
            "position": {"x": float(pos[n][0] * 1000),
                         "y": float(pos[n][1] * 1000), "z": 0.0},
            "color": {"r": r, "g": g, "b": b, "a": 1.0},
            "size": 10.0}
    nx.write_gexf(G, NET_DIR / "similarity_network.gexf")
    nx.to_pandas_edgelist(G).rename(columns={"source": "Source",
                                             "target": "Target",
                                             "weight": "Weight"}) \
        .to_csv(NET_DIR / "edges.csv", index=False)
    nodes.assign(Id=nodes.index, Label=nodes.gene)[
        ["Id", "Label", "family", "pfam", "fold_class", "Length", "n_tm",
         "topology"]].to_csv(NET_DIR / "nodes.csv", index=False)
    comps = sorted(nx.connected_components(strong), key=len, reverse=True)
    with open(NET_DIR / "clusters.txt", "w") as fh:
        fh.write(f"# connected components at prob >= {args.cluster_prob}\n")
        for k, c in enumerate(comps, 1):
            fh.write(f"cluster {k} ({len(c)}): "
                     + ", ".join(sorted(nodes.loc[list(c), "gene"])) + "\n")


if __name__ == "__main__":
    main()
