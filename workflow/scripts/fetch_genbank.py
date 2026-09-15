"""Download full GenBank records as XML, in resumable batches.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import refdb_entrez


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--uids", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    entrez = json.load(open(args.config_json, encoding="utf-8"))["entrez"]
    batch = int(entrez.get("batch_size", 200))
    uids = [u.strip() for u in open(args.uids, encoding="utf-8") if u.strip()]
    os.makedirs(args.outdir, exist_ok=True)

    client = refdb_entrez.client_from_params(entrez)
    for start in range(0, len(uids), batch):
        path = os.path.join(args.outdir, "gb_%06d.xml" % start)
        if os.path.exists(path) and os.path.getsize(path) > 500:
            continue
        with open(path, "wb") as fh:
            fh.write(client.efetch("nuccore", uids[start:start + batch]))
        print("fetched %d/%d" % (min(start + batch, len(uids)), len(uids)), flush=True)


if __name__ == "__main__":
    main()
