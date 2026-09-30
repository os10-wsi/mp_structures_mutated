"""Download AlphaFold models and UniProt entries for every protein in the ID mapping.

Needs network access to alphafold.ebi.ac.uk and rest.uniprot.org.
Already-downloaded files are skipped, so the script can be re-run safely.
"""
import json
import sys
import time

import requests

from common import STRUCT_DIR, UNIPROT_DIR, load_proteins

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{acc}.json"


def get(url, retries=4):
    for i in range(retries):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if i == retries - 1:
                raise
            print(f"  retry {url}: {e}", file=sys.stderr)
            time.sleep(2 ** (i + 1))


def fetch_structure(acc):
    out = STRUCT_DIR / f"{acc}.pdb"
    if out.exists():
        return True
    r = get(AFDB_API.format(acc=acc))
    if r is None:
        print(f"  {acc}: no AlphaFold model", file=sys.stderr)
        return False
    model = r.json()[0]
    pdb = get(model["pdbUrl"])
    out.write_bytes(pdb.content)
    return True


def fetch_uniprot(acc):
    out = UNIPROT_DIR / f"{acc}.json"
    if out.exists():
        return True
    r = get(UNIPROT_API.format(acc=acc))
    if r is None:
        print(f"  {acc}: not found in UniProt", file=sys.stderr)
        return False
    out.write_text(json.dumps(r.json()))
    return True


def main():
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    UNIPROT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_proteins()
    missing = []
    for acc, gene in zip(df["Entry"], df["gene"]):
        print(f"{gene:10s} {acc}")
        ok_s = fetch_structure(acc)
        ok_u = fetch_uniprot(acc)
        if not (ok_s and ok_u):
            missing.append(acc)
    print(f"\n{len(df) - len(missing)}/{len(df)} proteins complete")
    if missing:
        print("missing:", " ".join(missing))


if __name__ == "__main__":
    main()
