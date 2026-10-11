"""AQR1: GFP surface localisation of point mutants vs their DMS fitness.

Inputs : data/aqr1_gfp_localisation.tsv   membrane/internal GFP ratio per mutant
                                           (mean, ± error bar, lower, upper, significance)
         the AQR1 DMS table, given on the command line. Two layouts are recognised:
           * the mpdms canonical table, data/processed/AQR1_fitness.parquet
             (columns pos, wt, mut, vclass, pass_filter, score_z, score_raw, se, rep1..N)
           * the lab's fitness_estimation.tsv / fitness_singles_aqr1*.txt layout
             (columns 'wt aa', 'pos', 'mut aa', 'mean fitness', 'fitness sd', 'wt')
Outputs: results/aqr1_gfp_vs_fitness.tsv          merged table, one row per mutant
         results/figures/aqr1_gfp_vs_fitness.{png,pdf}

Usage:
    python 06_aqr1_gfp_vs_fitness.py /path/to/AQR1_fitness.parquet
    python 06_aqr1_gfp_vs_fitness.py /path/to/fitness_singles_aqr1_003.txt --score "mean fitness"

Each GFP row names a mutant (Q149A, Q238A, R504A, or WT). The script looks the
mutant up in the DMS table by position, wild-type and mutant residue, and plots
DMS fitness (x, ± its error: `se` or `fitness sd`) against the membrane/internal GFP ratio
(y, with the error bar from the GFP table). The WT point takes the WT row of the
DMS table if there is one, otherwise the median of the synonymous variants (which
is 0 by construction for score_z). Pearson and Spearman correlations are printed
and written on the figure; with four points they are descriptive only.
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import CATEGORICAL, DATA, FIGURES, GRID, RESULTS, TEXT, TEXT_2, style_axes

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.titlesize": 11, "axes.titleweight": "bold",
                     "axes.labelcolor": TEXT_2, "text.color": TEXT})

GFP_TABLE = DATA / "aqr1_gfp_localisation.tsv"
MUTANT_RE = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY\*])\b")


# ------------------------------------------------------------------ GFP table
def load_gfp(path=GFP_TABLE):
    """Return the GFP table with tidy column names and a parsed `variant` column."""
    sep = "\t" if Path(path).suffix.lower() in (".tsv", ".txt") else ","
    raw = pd.read_csv(path, sep=sep, keep_default_na=False)
    cols = {c: c.strip().lower() for c in raw.columns}
    pick = {}
    for c, lc in cols.items():
        if lc.startswith("condition"):
            pick["condition"] = c
        elif lc.startswith("mean"):
            pick["gfp_mean"] = c
        elif lc.startswith("error"):
            pick["gfp_err"] = c
        elif lc == "lower":
            pick["gfp_lower"] = c
        elif lc == "upper":
            pick["gfp_upper"] = c
        elif lc.startswith("significance"):
            pick["significance"] = c
    missing = {"condition", "gfp_mean"} - set(pick)
    if missing:
        sys.exit(f"GFP table {path}: cannot find column(s) {sorted(missing)}")
    df = raw.rename(columns={v: k for k, v in pick.items()})
    for c in ("gfp_mean", "gfp_err", "gfp_lower", "gfp_upper"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "gfp_lower" not in df:
        df["gfp_lower"] = df["gfp_mean"] - df["gfp_err"]
    if "gfp_upper" not in df:
        df["gfp_upper"] = df["gfp_mean"] + df["gfp_err"]
    if "significance" not in df:
        df["significance"] = ""
    df["significance"] = df["significance"].fillna("").astype(str).str.strip()

    def variant(cond):
        m = MUTANT_RE.search(str(cond))
        if m:
            return f"{m.group(1)}{int(m.group(2))}{m.group(3)}"
        if re.search(r"\bWT\b|wild[- ]?type", str(cond), flags=re.I):
            return "WT"
        return None

    df["variant"] = df["condition"].map(variant)
    bad = df[df["variant"].isna()]
    if len(bad):
        sys.exit("GFP table: cannot read a mutant (e.g. Q149A) or 'WT' from: "
                 + ", ".join(bad["condition"].astype(str)))
    return df


# ------------------------------------------------------------------ DMS table
def read_table(path):
    path = Path(path)
    ext = path.suffix.lower()
    if ext in (".parquet", ".pq"):
        return pd.read_parquet(path)
    if ext == ".csv":
        return pd.read_csv(path)
    # .tsv / .txt / anything else: let pandas sniff the delimiter
    return pd.read_csv(path, sep=None, engine="python")


def _is_true(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t", "y"])


def canonical_dms(df, score):
    """mpdms canonical table -> (variants, wt_fitness, score_label).

    variants has columns pos, wt, mut, fitness, fitness_err.
    """
    score = score or ("score_z" if "score_z" in df else "score_raw")
    if score not in df:
        sys.exit(f"DMS table has no column {score!r}; columns: {list(df.columns)}")
    d = df.copy()
    if "pass_filter" in d:
        d = d[d["pass_filter"].astype(bool)]
    reps = sorted([c for c in d.columns if re.fullmatch(r"rep\d+", c)], key=lambda c: int(c[3:]))
    if score == "score_raw":
        reps = [c + "_raw" for c in reps if c + "_raw" in d]
    if "se" in d and score == "score_z" and d["se"].notna().any():
        se = d["se"]
    elif "se_raw" in d and score == "score_raw" and d["se_raw"].notna().any():
        se = d["se_raw"]
    elif reps:
        se = d[reps].std(axis=1, ddof=1) / np.sqrt(d[reps].notna().sum(axis=1))
    else:
        se = pd.Series(np.nan, index=d.index)
    out = pd.DataFrame({"pos": pd.to_numeric(d["pos"], errors="coerce"),
                        "wt": d["wt"].astype(str).str.upper(),
                        "mut": d["mut"].astype(str).str.upper(),
                        "fitness": pd.to_numeric(d[score], errors="coerce"),
                        "fitness_err": se})
    if "vclass" in d:
        syn = out.loc[d["vclass"] == "synonymous", "fitness"]
    else:
        syn = out.loc[out["mut"] == "=", "fitness"]
    wt = float(syn.median()) if syn.notna().any() else (0.0 if score == "score_z" else np.nan)
    wt_src = "median of synonymous variants" if syn.notna().any() else "0 by definition"
    out = out[out["mut"].str.fullmatch(r"[A-Z]|\*")]
    return out, wt, wt_src, score


def lab_dms(df, score):
    """fitness_estimation.tsv layout -> (variants, wt_fitness, score_label)."""
    lc = {c: c.strip().lower() for c in df.columns}
    inv = {v: k for k, v in lc.items()}
    need = ["wt aa", "pos", "mut aa"]
    if not all(k in inv for k in need):
        sys.exit("DMS table: expected either mpdms columns (pos, wt, mut, score_z) or lab columns "
                 f"('wt aa', 'pos', 'mut aa', 'mean fitness'); got {list(df.columns)}")
    score = score or ("mean fitness" if "mean fitness" in inv else None)
    if score is None or score.lower() not in inv:
        sys.exit(f"DMS table has no column {score!r}; columns: {list(df.columns)}")
    score_col = inv[score.lower()]
    fitness = pd.to_numeric(df[score_col], errors="coerce")
    if "fitness sd" in inv and score.lower() == "mean fitness":
        se = pd.to_numeric(df[inv["fitness sd"]], errors="coerce")
    else:
        reps = [c for c in df.columns if re.fullmatch(r"rescaled_fitness_rep\d+", c.strip())]
        rep = df[reps].apply(pd.to_numeric, errors="coerce") if reps else None
        se = rep.std(axis=1, ddof=1) / np.sqrt(rep.notna().sum(axis=1)) if reps else np.nan
    pos = pd.to_numeric(df[inv["pos"]], errors="coerce")
    wt_aa = df[inv["wt aa"]].fillna("").astype(str).str.strip().str.upper()
    mut_aa = df[inv["mut aa"]].fillna("").astype(str).str.strip().str.upper()
    if "stop" in inv:
        mut_aa = mut_aa.where(~_is_true(df[inv["stop"]]), "*")
    # WT fitness: the explicit WT row, else the synonymous rows (aa_ham == 0, no position)
    wt_rows = _is_true(df[inv["wt"]]) if "wt" in inv else pd.Series(False, index=df.index)
    if wt_rows.any() and fitness[wt_rows].notna().any():
        wt, wt_src = float(fitness[wt_rows].mean()), "WT row of the DMS table"
    else:
        syn = (pd.to_numeric(df[inv["aa_ham"]], errors="coerce") == 0) if "aa_ham" in inv else pos.isna()
        syn &= ~wt_rows
        wt = float(fitness[syn].median()) if fitness[syn].notna().any() else np.nan
        wt_src = "median of synonymous variants"
    out = pd.DataFrame({"pos": pos, "wt": wt_aa, "mut": mut_aa, "fitness": fitness, "fitness_err": se})
    out = out[out["pos"].notna() & out["mut"].str.fullmatch(r"[A-Z]|\*")]
    return out, wt, wt_src, score


def load_dms(path, score=None):
    df = read_table(path)
    if {"pos", "wt", "mut"} <= set(df.columns) and ("score_z" in df or "score_raw" in df):
        return canonical_dms(df, score)
    return lab_dms(df, score)


# ------------------------------------------------------------------ merge + plot
def merge(gfp, variants, wt_fitness, wt_src):
    rows = []
    for r in gfp.itertuples(index=False):
        rec = r._asdict()
        if r.variant == "WT":
            rec.update(pos=np.nan, wt="", mut="", fitness=wt_fitness, fitness_err=np.nan,
                       fitness_source=wt_src)
        else:
            m = MUTANT_RE.fullmatch(r.variant)
            pos, wt, mut = int(m.group(2)), m.group(1), m.group(3)
            hit = variants[(variants["pos"] == pos) & (variants["mut"] == mut)]
            wt_mismatch = hit[hit["wt"] != wt]
            hit = hit[hit["wt"] == wt]
            if len(hit) == 0:
                src = "not in DMS table"
                if len(wt_mismatch):
                    src += f" (position {pos} is {wt_mismatch['wt'].iloc[0]} in the DMS table, not {wt})"
                print(f"warning: {r.variant}: {src}", file=sys.stderr)
                rec.update(pos=pos, wt=wt, mut=mut, fitness=np.nan, fitness_err=np.nan,
                           fitness_source=src)
            else:
                # several rows (e.g. several codons): average them
                rec.update(pos=pos, wt=wt, mut=mut,
                           fitness=float(hit["fitness"].mean()),
                           fitness_err=float(hit["fitness_err"].mean()) if hit["fitness_err"].notna().any() else np.nan,
                           fitness_source=f"DMS table ({len(hit)} row{'s' if len(hit) > 1 else ''})")
        rows.append(rec)
    cols = ["condition", "variant", "pos", "wt", "mut", "fitness", "fitness_err", "fitness_source",
            "gfp_mean", "gfp_err", "gfp_lower", "gfp_upper", "significance"]
    return pd.DataFrame(rows)[cols]


def correlations(m):
    d = m[m["fitness"].notna() & m["gfp_mean"].notna()]
    if len(d) < 3:
        return {"n": len(d), "pearson_r": np.nan, "pearson_p": np.nan,
                "spearman_rho": np.nan, "spearman_p": np.nan}
    pr = stats.pearsonr(d["fitness"], d["gfp_mean"])
    sr = stats.spearmanr(d["fitness"], d["gfp_mean"])
    return {"n": len(d), "pearson_r": float(pr[0]), "pearson_p": float(pr[1]),
            "spearman_rho": float(sr[0]), "spearman_p": float(sr[1])}


def plot(m, corr, score_label, out_stem):
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    style_axes(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    d = m[m["fitness"].notna()]
    for i, r in enumerate(d.itertuples(index=False)):
        color = CATEGORICAL[i % len(CATEGORICAL)]
        yerr = [[r.gfp_mean - r.gfp_lower], [r.gfp_upper - r.gfp_mean]] if np.isfinite(r.gfp_lower) else None
        xerr = r.fitness_err if np.isfinite(r.fitness_err) else None
        ax.errorbar(r.fitness, r.gfp_mean, xerr=xerr, yerr=yerr, fmt="o", ms=7, color=color,
                    ecolor=color, elinewidth=1, capsize=3, mec="white", mew=0.8, zorder=3)
        label = r.variant if r.significance in ("", "n/a", "na", "nan") else f"{r.variant}  {r.significance}"
        ax.annotate(label, (r.fitness, r.gfp_mean), xytext=(7, 6), textcoords="offset points",
                    fontsize=8.5, color=TEXT, zorder=4)
    if corr["n"] >= 3:
        x = d["fitness"].to_numpy(float)
        slope, intercept = np.polyfit(x, d["gfp_mean"].to_numpy(float), 1)
        xs = np.array([x.min(), x.max()])
        ax.plot(xs, slope * xs + intercept, color=TEXT_2, lw=1, ls="--", zorder=2)
        txt = (f"Pearson r = {corr['pearson_r']:.2f} (p = {corr['pearson_p']:.2f})\n"
               f"Spearman ρ = {corr['spearman_rho']:.2f} (p = {corr['spearman_p']:.2f})\n"
               f"n = {corr['n']}")
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=8, color=TEXT_2)
    missing = m[m["fitness"].isna() & (m["variant"] != "WT")]
    if len(missing):
        ax.text(0.98, 0.02, "not in DMS table: " + ", ".join(missing["variant"]), transform=ax.transAxes,
                va="bottom", ha="right", fontsize=7.5, color=TEXT_2)
    ax.set_xlabel(f"DMS fitness ({score_label})")
    ax.set_ylabel("Membrane / internal GFP")
    ax.set_title("AQR1: GFP surface localisation vs DMS fitness")
    ax.margins(x=0.25, y=0.2)
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dms", help="AQR1 DMS table: mpdms parquet, or the lab fitness tsv/txt/csv")
    ap.add_argument("--gfp", default=GFP_TABLE, help=f"GFP localisation table (default {GFP_TABLE})")
    ap.add_argument("--score", default=None,
                    help="fitness column to use (default: score_z for an mpdms table, 'mean fitness' "
                         "for the lab layout)")
    ap.add_argument("--wt-fitness", type=float, default=None,
                    help="override the fitness value used for the WT point")
    ap.add_argument("--name", default="aqr1_gfp_vs_fitness", help="output file stem")
    args = ap.parse_args()

    gfp = load_gfp(args.gfp)
    variants, wt_fitness, wt_src, score_label = load_dms(args.dms, args.score)
    if args.wt_fitness is not None:
        wt_fitness, wt_src = args.wt_fitness, "--wt-fitness"
    m = merge(gfp, variants, wt_fitness, wt_src)
    corr = correlations(m)

    RESULTS.mkdir(parents=True, exist_ok=True)
    table = RESULTS / f"{args.name}.tsv"
    m.to_csv(table, sep="\t", index=False, float_format="%.4g")
    plot(m, corr, score_label, FIGURES / args.name)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(m[["variant", "fitness", "fitness_err", "fitness_source", "gfp_mean", "gfp_err", "significance"]]
              .to_string(index=False))
    if corr["n"] >= 3:
        print(f"\nPearson r = {corr['pearson_r']:.3f} (p = {corr['pearson_p']:.3f}), "
              f"Spearman rho = {corr['spearman_rho']:.3f} (p = {corr['spearman_p']:.3f}), n = {corr['n']}")
    else:
        print(f"\nonly {corr['n']} mutants with both values: no correlation computed")
    print(f"wrote {table} and {FIGURES / args.name}.png/.pdf")


if __name__ == "__main__":
    main()
