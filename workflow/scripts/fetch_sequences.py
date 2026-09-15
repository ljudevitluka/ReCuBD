"""Download the sequence of every representative record.

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
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    entrez = json.load(open(args.config_json, encoding="utf-8"))["entrez"]
    reps = refdb_io.read_table(args.reps)
    ids = [v or a for v, a in zip(reps.version, reps.accession)]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    client = refdb_entrez.client_from_params(entrez)
    chunk = max(1, min(int(entrez.get("batch_size", 200)), 100))
    with open(args.out, "wb") as out:
        for start in range(0, len(ids), chunk):
            out.write(client.efetch("nuccore", ids[start:start + chunk], rettype="fasta", retmode="text"))
            print("sequences %d/%d" % (min(start + chunk, len(ids)), len(ids)), flush=True)

    got = sum(1 for line in open(args.out, encoding="utf-8") if line.startswith(">"))
    if got != len(ids):
        sys.exit("expected %d sequences, retrieved %d" % (len(ids), got))
    print("sequences retrieved: %d" % got)


if __name__ == "__main__":
    main()
