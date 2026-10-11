# Structural summary of S. cerevisiae membrane proteins

Input: `data/idmapping_2026_09_30.txt`, a UniProt ID-mapping table with 57 rows
and **56 unique proteins**. The query `ctr1` matches both CTR1 (P49573) and
HNM1 (P19807), because CTR1 is also a gene synonym of HNM1. Rows are
de-duplicated on accession, so the 56 proteins fill an 8 × 7 grid exactly.

## Pipeline

```bash
conda env create -f environment.yml && conda activate mp_structures
./run_all.sh
```

| Step | Script | Output |
|---|---|---|
| 1 | `01_fetch_data.py` | `data/structures/<ACC>.pdb` (AlphaFold DB), `data/uniprot/<ACC>.json` |
| 2 | `02_topology_summary.py` | `results/topology_summary.tsv`, `results/topology_classes.tsv`, figures `length_histogram`, `tm_count`, `topology_classes`, `topology_map` |
| 3 | `03_foldseek_allvsall.sh` | `results/foldseek/foldseek_easy_allvsall` |
| 4 | `04_similarity_network.py` | `results/network/similarity_network.gexf` (Gephi), `prob_matrix.tsv`, `nodes.csv`/`edges.csv`, `clusters.txt`, figure `similarity_network` |
| 5 | `05_structure_grid.py` | `results/renders/<ACC>.png`, figure `structure_grid` (8 × 7, rainbow N→C) |
| 6 | `06_aqr1_gfp_vs_fitness.py` | `results/aqr1_gfp_vs_fitness.tsv`, figure `aqr1_gfp_vs_fitness` (run on its own, see below) |

Figures go in `results/figures/`, as both PNG (300 dpi) and PDF.

### Transmembrane topology
The topology comes from the UniProt *Transmembrane*, *Intramembrane* and
*Topological domain* features. Each protein gets a class of the form
`<n> TM | N-in/out / C-in/out`, where *in* means cytoplasmic. If only one
terminus is annotated, the other is inferred from whether the TM count is even
or odd. If neither is annotated, the side is shown as `?`.

### Similarity network (Domainome protocol)
```
foldseek createdb data/structures foldseek_domainome_db
foldseek easy-search data/structures/*.pdb foldseek_domainome_db foldseek_easy_allvsall tmp \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,prob" --exhaustive-search 1 -e inf
```
Edges are Foldseek hit probabilities. Each edge takes the mean of the two
search directions, self-hits are removed, and the default cutoff is
`--min-prob 0.5`. Nodes are coloured by Pfam family.

- **Gephi, as in the paper:** open `results/network/similarity_network.gexf`.
  Run Layout → *Fruchterman Reingold* until the layout stops moving, then
  *ForceAtlas 2* with *Prevent Overlap* ticked. The file already stores
  colours, sizes and a starting layout.
- **Python:** the script runs the same two steps with networkx
  (`spring_layout` = Fruchterman-Reingold, then `forceatlas2_layout`, then
  a pass that pushes overlapping nodes apart) and saves the figure.

### Structure gallery
PyMOL cartoon coloured with `spectrum count, rainbow`: blue at the N-terminus,
red at the C-terminus. When UniProt TM helices are known, each model is rotated
so the membrane normal is vertical, with the cytoplasm at the bottom. Panels
are grouped by Pfam family, and each has a coloured bar above it. Long
disordered tails can make a model look small in its panel; use
`--min-plddt 50` to hide low-confidence termini.

### AQR1: GFP surface localisation vs DMS fitness
`data/aqr1_gfp_localisation.tsv` holds the membrane/internal GFP ratio of the
Aqr1 point mutants (mean, ± error bar, lower/upper bound, significance vs WT).
The script looks each mutant up in the AQR1 DMS table by position, wild-type
and mutant residue (Q149A, Q238A, R504A, all checked against UniProt P53943)
and plots DMS fitness (x, ± its error) against the GFP ratio (y, with the
error bar from the table), with Pearson and Spearman correlations. With four
points the correlations are descriptive only. The DMS table is not part of
this repository; pass it on the command line:

```bash
cd scripts
python 06_aqr1_gfp_vs_fitness.py /path/to/AQR1_fitness.parquet             # mpdms table, uses score_z
python 06_aqr1_gfp_vs_fitness.py /path/to/fitness_singles_aqr1_003.txt     # lab layout, uses "mean fitness"
python 06_aqr1_gfp_vs_fitness.py ... --score score_raw --wt-fitness 0      # other column / WT value
```

Both the mpdms canonical table (`pos`, `wt`, `mut`, `score_z`, `se`, `rep1..N`)
and the lab `fitness_estimation.tsv` layout (`wt aa`, `pos`, `mut aa`,
`mean fitness`, `fitness sd`, `wt` row) are recognised. The WT point takes the
WT row of the DMS table when there is one, otherwise the median of the
synonymous variants. A mutant missing from the DMS table is reported and left
out of the correlation. To add more conditions, append rows to the TSV in the
same `Aqr1 (X123Y)` format.

Set `MP_DATA` / `MP_RESULTS` to use other input and output directories.

## Results (run 2026-09-30, AlphaFold DB + UniProt)

- **Lengths:** 263–884 aa, median 568 aa (`length_histogram`).
- **Topologies:** there are 18 distinct classes and 13 distinct TM counts
  (`topology_classes`, `topology_map`). The largest class is 12 TM, N-in/C-in,
  with 24 proteins (most MFS and APC transporters). Next come 11 TM,
  N-out/C-in (6), then five classes of 3 proteins each. These are UniProt's
  sequence-based annotations: an MFS protein listed with 11 or 13 TM probably
  has one helix missed or split.
- **Structural similarity network** (prob ≥ 0.5, 310 edges, 15 components):
  - MFS: Sugar_tr + MFS_1, 22 proteins.
  - LeuT fold: APC permeases + DUR3 + FUR4, with BOR1 joined by a single
    weak edge, 13 proteins.
  - Smaller groups: SUR7 (4), Gpr1/Fun34/YaaH (3), RTA1 + IZH2 (3), KCH (2).
  - 9 singletons.
  See `results/network/clusters.txt`.
- **Structure gallery:** `structure_grid` shows every residue;
  `structure_grid_plddt50` hides terminal residues with pLDDT < 50.
