"""Score completeness and pick one representative per taxon per segment.

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

RANK = {"complete genome": 0, "coding-complete": 1, "near-complete": 2,
        "substantial fragment": 3, "gene fragment": 4, "unscored": 5}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--augment", required=True)
    ap.add_argument("--tier2", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--out-pooled", required=True)
    ap.add_argument("--out-reps", required=True)
    args = ap.parse_args()

    thresholds = json.load(open(args.config_json, encoding="utf-8"))["curation"]["completeness"]
    near_min = float(thresholds.get("near_complete_min_pct", 85))
    frag_min = float(thresholds.get("substantial_fragment_min_pct", 50))

    groups = refdb_io.load_tier2_groups(args.tier2)
    name_re = refdb_io.term_regex(refdb_io.load_search_terms(args.terms)["name"])

    df = refdb_io.read_table(args.annotated)
    aug = refdb_io.read_table(args.augment)

    if len(aug):
        aug = refdb_gb.add_viral_taxonomy(aug)
        tokens = aug.lineage.apply(lambda s: set(str(s).split("; ")))
        aug["is_phage"] = tokens.apply(lambda t: bool(t & groups["phage"]))
        aug["in_tier2_group"] = tokens.apply(lambda t: bool(t & groups["environmental"]))
        aug["host_named_virus"] = aug.organism.str.contains(name_re)
        aug["tier"] = "tier1_host_and_related"
        aug.loc[aug.in_tier2_group & ~aug.host_named_virus, "tier"] = "tier2_environmental_diet_microbiome"
        aug["probable_host_class"] = "named after the host taxon or a relative"
        aug["host_evidence"] = "species-level reference (isolate from non-target host)"
        aug["host_scientific_name"] = aug.host.str.replace(r"\s*\(.*?\)", "", regex=True).str.strip()
        aug["host_taxid"] = ""
        aug["host_reject_reason"] = ""
        aug["exclusion_reason"] = ""
        aug["included"] = ~aug.is_phage
        aug["from_host_taxon"] = False
        aug["refseq"] = aug.accession.str.match(r"^(NC_|AC_|NG_)") | aug.db_source.str.upper().str.contains("REFSEQ")
        aug["placeholder_taxon"] = aug.organism.str.contains(r"\bsp\.$|^unclassified", case=False, regex=True)
        aug = aug[~aug.placeholder_taxon]   # a shared '... sp.' taxid pools unrelated viruses
        pool = pd.concat([df, aug], ignore_index=True).drop_duplicates("accession")
    else:
        pool = df.copy()

    pool["tax_group"] = pool.virus_family.where(
        pool.virus_family != "unclassified",
        pool.virus_order.where(pool.virus_order != "", pool.lineage_terminal))
    defs = pool.definition.str.lower()
    pool["def_complete"] = defs.str.contains(r"complete genome|complete sequence|genome assembly")
    pool["def_cds_complete"] = defs.str.contains(r"complete cds|polyprotein gene, complete")

    anchored = pool[pool.def_complete | pool.def_cds_complete].groupby("tax_group").length.max()
    fallback = pool.groupby("tax_group").length.quantile(0.95)
    ref_len = anchored.reindex(fallback.index).fillna(fallback).replace(0, pd.NA)
    pool["group_reference_length"] = pool.tax_group.map(ref_len)
    pool["pct_of_group_reference"] = (100 * pool.length / pool.group_reference_length).round(1)

    def status(row):
        if row.def_complete:
            return "complete genome"
        if row.def_cds_complete:
            return "coding-complete"
        pct = row.pct_of_group_reference
        if pd.isna(pct):
            return "unscored"
        if pct >= near_min:
            return "near-complete"
        if pct >= frag_min:
            return "substantial fragment"
        return "gene fragment"

    pool["genome_status"] = pool.apply(status, axis=1)
    refdb_io.write_table(pool, args.out_pooled)

    sel = pool[pool.included].copy()
    sel["segment_key"] = sel.segment.replace("", "-")
    sel["group_key"] = sel.taxid.where(~sel.placeholder_taxon, sel.accession)
    sel["status_rank"] = sel.genome_status.map(RANK)
    reps = (sel.sort_values(["status_rank", "refseq", "from_host_taxon", "length"],
                            ascending=[True, False, False, False])
               .groupby(["group_key", "segment_key"], as_index=False).first())
    refdb_io.write_table(reps, args.out_reps)

    print("pooled records: %d | representatives: %d" % (len(pool), len(reps)))
    print(reps.genome_status.value_counts().to_string())
    print("representatives from a non-target host: %d rows, %d taxa"
          % (int((~reps.from_host_taxon).sum()), reps.loc[~reps.from_host_taxon, "taxid"].nunique()))


if __name__ == "__main__":
    main()
