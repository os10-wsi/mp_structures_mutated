#!/usr/bin/env bash
# Full pipeline. Needs network access for step 1 only.
set -euo pipefail
cd "$(dirname "$0")/scripts"
python 01_fetch_data.py          # AlphaFold models + UniProt topology
python 02_topology_summary.py    # TM topologies, length histogram
./03_foldseek_allvsall.sh        # Foldseek all-vs-all
python 04_similarity_network.py  # network figure + Gephi GEXF
python 05_structure_grid.py      # 8 x 7 rainbow structure gallery
