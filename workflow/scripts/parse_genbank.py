"""Parse GenBank XML into a flat record table.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import refdb_gb
import refdb_io


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xmldir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    records = refdb_gb.parse_dir(args.xmldir)
    if not records:
        sys.exit("No GBSeq records parsed from " + args.xmldir)
    df = pd.DataFrame(records, columns=refdb_gb.COLUMNS).drop_duplicates("accession")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    refdb_io.write_table(df, args.out)
    print("parsed records: %d" % len(df))


if __name__ == "__main__":
    main()
