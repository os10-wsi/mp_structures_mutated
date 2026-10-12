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
search directions, and self-hits are removed. As in the paper, all hits are
kept as edges (`--min-prob 0.01`). Clusters are reported at `--cluster-prob
0.5` in `results/network/clusters.txt`. Nodes are coloured by class, in the
style of the Domainome figure:
- MFS fold (Sugar_tr, MFS_1)
- LeuT fold (APC: AA_permease, AA_permease_2, SSS/DUR3, NCS1/FUR4)
- other proteins with ≥ 8 TM
- proteins with ≤ 7 TM

`similarity_network` has cluster labels and class counts;
`similarity_network_nolabels` has no text.

- **Gephi, as in the paper:** open `results/network/similarity_network.gexf`.
  Run Layout → *Fruchterman Reingold* until the layout stops moving, then
  *ForceAtlas 2* with *Prevent Overlap* ticked. The file already stores
  colours and a starting layout.
- **Python:** the script runs the same two steps with networkx
  (`spring_layout` = Fruchterman-Reingold, then `forceatlas2_layout`, then
  a pass that pushes overlapping nodes apart). For the layout, edge weights
  are probability cubed, so the many weak hits don't smear the clusters
  together.

### Structure gallery
PyMOL cartoons in an 8 × 7 grid with no labels. By default each structure is
coloured by Pfam family, using the same blue/grey shades as the topology map;
the key is `structure_family_key`. With `--color rainbow` it is coloured N→C
instead (blue → red; files `structure_grid_rainbow*`, key
`structure_colour_key`). Panels are in row-major order, grouped by family, and
listed in `results/figures/structure_grid_order.tsv`. When UniProt TM helices
are known, each model is rotated so the membrane normal is vertical, with the
cytoplasm at the bottom. Long disordered tails can make a model look small in
its panel; use `--min-plddt 50` to hide low-confidence termini.

### Charts
The length, TM-count and topology-class charts are plain grey. In the
topology map, TM helices are coloured by Pfam family in shades of blue (the
transporter families) and grey (the rest).

Set `MP_DATA` / `MP_RESULTS` to use other input and output directories.

## Results (run 2026-09-30, AlphaFold DB + UniProt)

- **Lengths:** 263–884 aa, median 568 aa (`length_histogram`).
- **Topologies:** there are 18 distinct classes and 13 distinct TM counts
  (`topology_classes`, `topology_map`). The largest class is 12 TM, N-in/C-in,
  with 24 proteins (most MFS and APC transporters). Next come 11 TM,
  N-out/C-in (6), then five classes of 3 proteins each. These are UniProt's
  sequence-based annotations: an MFS protein listed with 11 or 13 TM probably
  has one helix missed or split.
- **Structural similarity network** (789 edges at prob ≥ 0.01; 15 clusters at prob ≥ 0.5):
  - MFS: Sugar_tr + MFS_1, 22 proteins.
  - LeuT fold: APC permeases + DUR3 + FUR4, with BOR1 joined by a single
    weak edge, 13 proteins.
  - Smaller groups: SUR7 (4), Gpr1/Fun34/YaaH (3), RTA1 + IZH2 (3), KCH (2).
  - 9 singletons.
  See `results/network/clusters.txt`.
- **Structure gallery:** `structure_grid` (family colours) and
  `structure_grid_rainbow` show every residue; the `_plddt50` versions hide
  terminal residues with pLDDT < 50.
