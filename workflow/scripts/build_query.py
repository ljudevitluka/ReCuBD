"""Build the Entrez query from the search-term file and resolve it to UIDs.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import refdb_entrez
import refdb_io


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True, help="effective config written by the workflow")
    ap.add_argument("--terms", required=True)
    ap.add_argument("--out-query", required=True)
    ap.add_argument("--out-uids", required=True)
    ap.add_argument("--out-meta", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.config_json, encoding="utf-8"))
    entrez = cfg["entrez"]
    terms = refdb_io.load_search_terms(args.terms)
    query = refdb_io.build_entrez_query(terms, entrez.get("virus_taxid_root", 10239))

    client = refdb_entrez.client_from_params(entrez)
    uids, total = client.esearch_uids("nuccore", query)

    os.makedirs(os.path.dirname(args.out_query) or ".", exist_ok=True)
    with open(args.out_query, "w", encoding="utf-8") as fh:
        fh.write(query + "\n")
    with open(args.out_uids, "w", encoding="utf-8") as fh:
        fh.write("\n".join(uids) + ("\n" if uids else ""))
    with open(args.out_meta, "w", encoding="utf-8") as fh:
        json.dump({"query": query, "hit_count": total, "uids_retrieved": len(uids),
                   "query_terms": len(terms["query"]), "name_terms": len(terms["name"]),
                   "virus_taxid_root": entrez.get("virus_taxid_root", 10239)}, fh, indent=2)

    print("Entrez hits: %d" % total)
    if total == 0:
        sys.exit("No records matched the query. Check " + args.terms + ".")


if __name__ == "__main__":
    main()
