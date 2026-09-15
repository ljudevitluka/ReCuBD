"""Fetch complete genomes of tier-1 taxa deposited from hosts outside the target taxon.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import refdb_entrez
import refdb_gb
import refdb_io


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.config_json, encoding="utf-8"))
    entrez, aug = cfg["entrez"], cfg["augmentation"]
    os.makedirs(args.outdir, exist_ok=True)
    df = refdb_io.read_table(args.annotated)

    if not aug.get("enabled", True):
        refdb_io.write_table(pd.DataFrame(columns=refdb_gb.COLUMNS), args.out)
        print("augmentation disabled")
        return

    keep = df.included & df.tier.str.startswith("tier1")
    if aug.get("skip_placeholder_taxa", True):
        keep &= ~df.placeholder_taxon
    taxids = sorted({t for t in df.loc[keep, "taxid"] if t})

    clause = " OR ".join('"' + t + '"[Title]' for t in aug.get("title_terms", ["complete genome"]))
    if aug.get("include_refseq", True):
        clause += " OR srcdb_refseq[PROP]"

    client = refdb_entrez.client_from_params(entrez)
    step = int(aug.get("taxids_per_query", 60))
    batch = int(entrez.get("batch_size", 200))
    for i in range(0, len(taxids), step):
        done = os.path.join(args.outdir, "chunk_%06d.done" % i)
        if os.path.exists(done):
            continue
        chunk = taxids[i:i + step]
        query = "(" + " OR ".join("txid" + t + "[Organism:exp]" for t in chunk) + ") AND (" + clause + ")"
        uids, _ = client.esearch_uids("nuccore", query)
        uids = uids[: int(aug.get("max_per_chunk", 400))]
        for j, start in enumerate(range(0, len(uids), batch)):
            with open(os.path.join(args.outdir, "aug_%06d_%03d.xml" % (i, j)), "wb") as fh:
                fh.write(client.efetch("nuccore", uids[start:start + batch]))
        open(done, "w", encoding="utf-8").write(str(len(uids)) + "\n")
        print("taxid chunk %d: %d candidate genomes" % (i, len(uids)), flush=True)

    records = refdb_gb.parse_dir(args.outdir)
    out = pd.DataFrame(records, columns=refdb_gb.COLUMNS) if records else pd.DataFrame(columns=refdb_gb.COLUMNS)
    out = out.drop_duplicates("accession")
    out = out[~out.accession.isin(set(df.accession))]
    refdb_io.write_table(out, args.out)
    print("augmentation candidates: %d for %d tier-1 taxa" % (len(out), len(taxids)))


if __name__ == "__main__":
    main()
