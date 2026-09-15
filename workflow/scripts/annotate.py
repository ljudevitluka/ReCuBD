"""Assign tier, probable host class and phage rejection from the tier-2 group file.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import refdb_io


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--tier2", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    curation = json.load(open(args.config_json, encoding="utf-8"))["curation"]
    groups = refdb_io.load_tier2_groups(args.tier2)
    name_re = refdb_io.term_regex(refdb_io.load_search_terms(args.terms)["name"])

    df = refdb_io.read_table(args.table)
    tokens = df.lineage.apply(lambda s: set(str(s).split("; ")))

    df["is_phage"] = tokens.apply(lambda t: bool(t & groups["phage"])) | \
        df.organism.str.contains(r"\bphage\b", case=False, regex=True)
    df["in_tier2_group"] = tokens.apply(lambda t: bool(t & groups["environmental"]))
    df["host_named_virus"] = df.organism.str.contains(name_re)

    df["tier"] = "tier1_host_and_related"
    df.loc[df.in_tier2_group & ~df.host_named_virus, "tier"] = "tier2_environmental_diet_microbiome"

    df["probable_host_class"] = "uncertain (unclassified or unlisted group)"
    df.loc[df.in_tier2_group & ~df.host_named_virus, "probable_host_class"] = "environmental, diet or microbiome"
    df.loc[df.host_named_virus, "probable_host_class"] = "named after the host taxon or a relative"
    df.loc[df.is_phage, "probable_host_class"] = "prokaryotic (microbiome)"

    df["exclusion_reason"] = df.host_reject_reason
    if curation.get("reject_phages", True):
        df.loc[df.is_phage & df.exclusion_reason.eq(""), "exclusion_reason"] = \
            "bacteriophage (true host is a bacterium in the host-taxon microbiome)"
    df["included"] = df.exclusion_reason.eq("") & df.from_host_taxon

    refdb_io.write_table(df, args.out)
    print("included: %d | excluded: %d" % (int(df.included.sum()), int((~df.included).sum())))
    print(df[df.included].tier.value_counts().to_string())


if __name__ == "__main__":
    main()
