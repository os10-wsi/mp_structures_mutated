"""Gallery of all structures, rainbow-coloured N (blue) -> C (red), 8 x 7 grid.

Each AlphaFold model is rendered with PyMOL as a cartoon coloured by residue
position. Where UniProt TM helices are known (results/topology_summary.tsv),
the model is rotated so the membrane normal is vertical, cytoplasm at the
bottom, which makes the transporters directly comparable; otherwise PyMOL's
default orientation (principal axes) is used.

Panels are ordered by Pfam family, then gene name, and labelled with the gene
name and UniProt accession.

Outputs: results/renders/<ACC>.png, results/figures/structure_grid.{png,pdf}
         (with --min-plddt X: renders_plddtX/, structure_grid_plddtX.*)

Usage: python 05_structure_grid.py [--cols 8] [--rows 7] [--min-plddt 0]
  --min-plddt X  hides terminal residues with pLDDT < X (long disordered tails
                 can dominate a small panel). Default 0 = show everything.
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol2
from matplotlib.colors import LinearSegmentedColormap

from common import (FIGURES, RESULTS, STRUCT_DIR, TEXT, TEXT_2, family_colors,
                    family_label, load_proteins)

PX = 900


def membrane_normal(ca, tm_segments):
    """Unit membrane normal from TM helix directions (first TM points +).
    `ca` maps resi -> CA xyz."""
    vecs = []
    for k, (s, e) in enumerate(tm_segments):
        if s in ca and e in ca:
            v = ca[e] - ca[s]
            vecs.append(v if k % 2 == 0 else -v)   # consecutive TMs antiparallel
    if not vecs:
        return None
    normal = np.sum(vecs, axis=0)
    return normal / np.linalg.norm(normal)


def orient_by_membrane(cmd, obj, tm_segments, n_term):
    ca = {}
    cmd.iterate_state(1, f"{obj} and name CA",
                      "ca[int(resv)] = (x, y, z)", space={"ca": ca})
    ca = {k: np.array(v) for k, v in ca.items()}
    normal = membrane_normal(ca, tm_segments)
    if normal is None:
        return False
    # first TM runs from the N-terminal side to the other side
    if n_term == "out":
        normal = -normal          # make +y point out of the cell
    # plain iterate/alter keeps this independent of PyMOL's NumPy ABI
    coords = []
    cmd.iterate_state(1, obj, "coords.append((x, y, z))",
                      space={"coords": coords})
    xyz = np.array(coords)
    centre = xyz.mean(axis=0)
    xyz -= centre
    inplane = xyz - np.outer(xyz @ normal, normal)
    _, _, vt = np.linalg.svd(inplane, full_matrices=False)
    xaxis = vt[0] - (vt[0] @ normal) * normal
    xaxis /= np.linalg.norm(xaxis)
    zaxis = np.cross(xaxis, normal)
    R = np.vstack([xaxis, normal, zaxis])        # rows = new axes
    new = iter((xyz @ R.T).tolist())
    cmd.alter_state(1, obj, "(x, y, z) = next(new)",
                    space={"new": new, "next": next})
    # identity view: camera looks down -z, y up
    view = list(cmd.get_view())
    view[:9] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    cmd.set_view(view)
    return True


def trim_termini(cmd, obj, min_plddt):
    plddt = {}
    cmd.iterate(f"{obj} and name CA", "plddt[int(resv)] = b",
                space={"plddt": plddt})
    res = sorted(plddt)
    good = [r for r in res if plddt[r] >= min_plddt]
    if not good:
        return
    if good[0] > res[0]:
        cmd.remove(f"{obj} and resi {res[0]}-{good[0] - 1}")
    if good[-1] < res[-1]:
        cmd.remove(f"{obj} and resi {good[-1] + 1}-{res[-1]}")


def render(acc, topo_row, min_plddt, render_dir, force=False):
    out = render_dir / f"{acc}.png"
    if out.exists() and not force:
        return out
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(str(STRUCT_DIR / f"{acc}.pdb"), "m")
        cmd.remove("hydro or solvent or hetatm")
        if min_plddt > 0:
            trim_termini(cmd, "m", min_plddt)
        cmd.hide("everything")
        cmd.show("cartoon")
        cmd.spectrum("count", "rainbow", "m and name CA")
        cmd.set("cartoon_transparency", 0)
        cmd.set("cartoon_fancy_helices", 1)
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)
        cmd.set("ray_shadows", 0)
        cmd.set("antialias", 2)
        cmd.set("ambient", 0.4)
        cmd.set("specular", 0.2)
        cmd.set("orthoscopic", 1)
        oriented = False
        if topo_row is not None and isinstance(topo_row.tm_segments, str):
            segs = [tuple(map(int, s.split("-")))
                    for s in topo_row.tm_segments.split(";") if s]
            oriented = orient_by_membrane(cmd, "m", segs, topo_row.n_term)
        if not oriented:
            cmd.orient("m")
        cmd.zoom("m", buffer=2, complete=1)
        cmd.png(str(out), width=PX, height=PX, dpi=300, ray=1)
    return out


def rainbow_cmap():
    """Sample PyMOL's own 'rainbow' spectrum so the colour bar matches."""
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.pseudoatom("cbar", pos=[0, 0, 0])
        for i in range(1, 100):
            cmd.pseudoatom("cbar", pos=[i, 0, 0], resi=i + 1)
        cmd.spectrum("count", "rainbow", "cbar")
        cols = []
        cmd.iterate("cbar", "cols.append(color)", space={"cols": cols})
        rgb = [cmd.get_color_tuple(c) for c in cols]
    return LinearSegmentedColormap.from_list("pymol_rainbow", rgb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=8)
    ap.add_argument("--rows", type=int, default=7)
    ap.add_argument("--min-plddt", type=float, default=0)
    ap.add_argument("--force", action="store_true", help="re-render all")
    args = ap.parse_args()

    df = load_proteins()
    df["family"] = df.pfam.map(family_label)
    fam_order = {f: i for i, f in enumerate(family_colors())}
    df = df[[(STRUCT_DIR / f"{a}.pdb").exists() for a in df.Entry]]
    df = (df.assign(fo=df.family.map(fam_order))
          .sort_values(["fo", "gene"]).reset_index(drop=True))
    topo_path = RESULTS / "topology_summary.tsv"
    topo = (pd.read_csv(topo_path, sep="\t").set_index("Entry")
            if topo_path.exists() else None)
    if len(df) > args.cols * args.rows:
        print(f"warning: {len(df)} structures but only "
              f"{args.cols * args.rows} panels; extra ones are dropped")

    tag = f"_plddt{args.min_plddt:g}" if args.min_plddt > 0 else ""
    render_dir = RESULTS / f"renders{tag}"
    render_dir.mkdir(parents=True, exist_ok=True)
    for acc in df.Entry:
        row = topo.loc[acc] if topo is not None and acc in topo.index else None
        render(acc, row, args.min_plddt, render_dir, args.force)
        print("rendered", acc)

    colors = family_colors()
    fig, axes = plt.subplots(args.rows, args.cols,
                             figsize=(args.cols * 1.9, args.rows * 2.15 + 0.6))
    for ax in axes.flat:
        ax.axis("off")
    for ax, r in zip(axes.flat, df.itertuples()):
        ax.imshow(mpimg.imread(render_dir / f"{r.Entry}.png"))
        ax.set_title(r.gene, fontsize=9, fontweight="bold", color=TEXT, pad=7)
        ax.text(0.5, -0.02, f"{r.Entry} · {r.Length} aa", ha="center",
                va="top", transform=ax.transAxes, fontsize=6, color=TEXT_2)
        ax.add_patch(plt.Rectangle((0.0, 1.0), 1.0, 0.025,
                                   transform=ax.transAxes, clip_on=False,
                                   color=colors[r.family], lw=0))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.07,
                        wspace=0.05, hspace=0.35)
    fig.suptitle("AlphaFold models of S. cerevisiae membrane proteins",
                 x=0.01, ha="left", fontsize=13, fontweight="bold", color=TEXT)

    # N -> C colour bar
    cax = fig.add_axes([0.68, 0.03, 0.28, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=rainbow_cmap()), cax=cax,
                      orientation="horizontal", ticks=[0, 1])
    cb.ax.set_xticklabels(["N-terminus", "C-terminus"], fontsize=8)
    cb.outline.set_visible(False)
    # family colour key (strip above each panel)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[f])
               for f in colors if f in set(df.family)]
    fig.legend(handles, [f for f in colors if f in set(df.family)],
               ncol=3, fontsize=7, frameon=False, loc="lower left",
               bbox_to_anchor=(0.01, 0.0), title="Pfam family (bar above panel)",
               title_fontsize=7)
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"structure_grid{tag}.{ext}", dpi=300)
    plt.close(fig)
    print(f"wrote {FIGURES / f'structure_grid{tag}.png'} ({len(df)} structures)")


if __name__ == "__main__":
    main()
