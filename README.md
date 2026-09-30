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
search directions, self-hits are removed, and the default cutoff is
`--min-prob 0.5`. Nodes are coloured by Pfam family.

- **Gephi, as in the paper:** open `results/network/similarity_network.gexf`.
  Run Layout → *Fruchterman Reingold* until the layout stops moving, then
  *ForceAtlas 2* with *Prevent Overlap* ticked. The file already stores
  colours, sizes and a starting layout.
- **Python:** the script runs the same two steps with networkx
  (`spring_layout` = Fruchterman-Reingold, then `forceatlas2_layout` with
  node sizes to prevent overlap) and saves the figure.

### Structure gallery
PyMOL cartoon coloured with `spectrum count, rainbow`: blue at the N-terminus,
red at the C-terminus. When UniProt TM helices are known, each model is rotated
so the membrane normal is vertical, with the cytoplasm at the bottom. Panels
are grouped by Pfam family, and each has a coloured bar above it. Long
disordered tails can make a model look small in its panel; use
`--min-plddt 50` to hide low-confidence termini.

Set `MP_DATA` / `MP_RESULTS` to use other input and output directories.
