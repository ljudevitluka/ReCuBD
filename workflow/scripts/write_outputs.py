"""Write the FASTA files, metadata tables and build report.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import datetime
import re
import textwrap
from string import Template

import pandas as pd

import refdb_io

TIER1 = "tier1_host_and_related"
TIER2 = "tier2_environmental_diet_microbiome"
COLS = ["accession", "version", "organism", "taxid", "virus_family", "virus_order", "virus_class",
        "lineage", "tier", "probable_host_class", "genome_status", "length",
        "pct_of_group_reference", "segment", "moltype", "topology", "n_cds", "products", "host",
        "host_scientific_name", "host_evidence", "from_host_taxon", "country", "collection_date",
        "strain", "refseq", "pubmed", "kmer_containment_to_host_record",
        "est_ani_to_host_record_pct", "compared_with_host_accession",
        "create_date", "update_date", "definition"]
RENAME = {"organism": "virus_name", "length": "length_bp", "host": "host_qualifier_raw",
          "host_evidence": "host_taxon_evidence"}

REPORT = """
# $dataset reference database

**Build date:** $date
**Source:** NCBI Nucleotide (GenBank + RefSeq) via E-utilities
**Built with:** host-virus-refdb (Snakemake). The scope is defined entirely by
`$terms_file`, `$accept_file` and `$tier2_file`.

## Entrez query

```
$query
```

$hits records matched and $parsed were parsed. Host curation kept $curated records from the target
host taxon and rejected $excluded. Representative selection returned $reps sequences ($n_t1 tier 1,
$n_t2 tier 2), of which $n_nontarget rows ($n_nontarget_taxa taxa) are species-level references whose
isolate came from a host outside the target taxon (`from_host_taxon = False`); the other
$n_target representatives are host-taxon-derived. $t1_usable of $n_t1 tier-1 representatives are
complete, coding-complete or near-complete.

$screen_txt

## Contents

| File | Records | Description |
|---|---|---|
| `${dataset}_refdb_tier1.fasta` | $n_t1 | Representatives plausibly infecting the host taxon or its relatives |
| `${dataset}_refdb_tier2.fasta` | $n_t2 | Groups listed as environmental in the tier-2 input file (diet, water, microbiome) |
| `${dataset}_refdb_all.fasta` | $reps | Both tiers |
| `${dataset}_reference_genomes.tsv` | $reps | Metadata for every sequence in the FASTA files |
| `${dataset}_all_records.tsv` | $curated | Every curated host-taxon record, before representative selection |
| `${dataset}_species_summary.tsv` | $n_species | One row per virus taxon: hosts, countries, years, representative |
| `${dataset}_excluded_records.tsv` | $excluded | Retrieved but rejected, with reasons (audit trail) |
| `tables/augment_similarity_audit.tsv` | $n_audit | Every screened candidate reference genome, with containment score and verdict |

FASTA headers: `accession.version | virus name | segment= | family= | genome status | host= | country= | length`

## Host coverage

$host_table

## Representative genome status

$status_table

## Why records were rejected

$exc_table

## Method

1. **Retrieval** — every `use=query`/`use=both` term in `$terms_file` is OR-joined into one Entrez
   query restricted to the taxid subtree $taxroot. Full GenBank XML is fetched so that host,
   isolation source, country and collection date come from source feature qualifiers rather than
   definition lines.
2. **Host curation** — every distinct `/host` string is resolved against NCBI Taxonomy and kept only
   if its lineage contains a clade from `$accept_file`, or the raw string matches one of that file's
   vernacular patterns. Host-less records are kept when the virus name carries a host term or the
   isolation source names the host taxon; symbiont sources are rejected.
3. **Tiering** — groups listed in `$tier2_file` as `environmental` go to tier 2; groups listed as
   `phage` are rejected outright. A virus whose name carries a host term stays in tier 1 even when
   its family is listed. Tier is a triage annotation inferred from taxonomy, not an experimental
   host determination.
4. **Completeness** — `complete genome`/`coding-complete` from the GenBank definition; otherwise a
   percentage of the longest definition-complete member of the same family (or order, for
   unclassified viruses): near-complete >= $near_min%, substantial fragment >= $frag_min%, else gene
   fragment.
5. **Representatives** — one sequence per taxon per segment, preferring complete > coding-complete >
   near-complete > fragment, then RefSeq, then host-taxon-derived, then longest. Placeholder taxa
   (`... sp.`) are kept per accession because that taxid is shared by unrelated viruses.
6. **Species-level references** — where a virus reported from the host taxon has only fragments from
   it, complete genomes deposited from other hosts are retrieved as candidates; those that win a
   representative slot are flagged `from_host_taxon = False`.

## Use

```bash
makeblastdb -in ${dataset}_refdb_tier1.fasta -dbtype nucl -parse_seqids -out ${dataset}_t1
blastn -query contigs.fasta -db ${dataset}_t1 -outfmt 6 -evalue 1e-10
```

Rebuild rather than edit: the input files plus the query above are the reproducible definition of
this set, and GenBank grows.
"""


def read_fasta(path):
    seqs, key = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            key = line[1:].split()[0]
            seqs[key] = []
        elif key and line:
            seqs[key].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def tidy(frame):
    return (frame.reindex(columns=COLS).rename(columns=RENAME)
                 .sort_values(["tier", "virus_family", "virus_name", "accession"]))


def joined(values):
    return "; ".join(sorted({str(v).strip() for v in values if str(v).strip()}))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--pooled", required=True)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--host-accept", required=True)
    ap.add_argument("--tier2", required=True)
    ap.add_argument("--audit", default="",
                    help="augment_similarity_audit.tsv from the screening step")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dataset", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.config_json, encoding="utf-8"))
    sim_cfg = cfg.get("similarity", {})
    min_cont = float(sim_cfg.get("min_containment", 0.5))
    kmer_size = int(sim_cfg.get("kmer_size", 21))
    pool = refdb_io.read_table(args.pooled)
    reps = refdb_io.read_table(args.reps)
    seqs = read_fasta(args.fasta)
    os.makedirs(args.outdir, exist_ok=True)

    def out(suffix):
        return os.path.join(args.outdir, args.dataset + suffix)

    reps["seq_key"] = [v or a for v, a in zip(reps.version, reps.accession)]
    missing = [k for k in reps.seq_key if k not in seqs]
    if missing:
        sys.exit("no sequence retrieved for: " + ", ".join(missing[:5]))
    reps["retrieved_length"] = reps.seq_key.map(lambda k: len(seqs[k]))
    bad = reps[reps.retrieved_length != reps.length]
    if len(bad):
        sys.exit("length mismatch for %d records, e.g. %s" % (len(bad), bad.accession.iloc[0]))

    def write_fasta(frame, path):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            for r in frame.itertuples():
                header = " | ".join([
                    r.seq_key, r.organism, "segment=" + (r.segment or "n/a"),
                    "family=" + r.virus_family, r.genome_status,
                    "host=" + (r.host_scientific_name or r.host or "n/a"),
                    "country=" + (r.country or "n/a"), str(r.length) + "bp"])
                fh.write(">" + header + "\n" + "\n".join(textwrap.wrap(seqs[r.seq_key], 70)) + "\n")

    write_fasta(reps[reps.tier == TIER1], out("_refdb_tier1.fasta"))
    write_fasta(reps[reps.tier == TIER2], out("_refdb_tier2.fasta"))
    write_fasta(reps, out("_refdb_all.fasta"))

    ref_tbl = tidy(reps)
    ref_tbl.to_csv(out("_reference_genomes.tsv"), sep="\t", index=False)
    curated = pool[pool.included & pool.from_host_taxon]
    tidy(curated).to_csv(out("_all_records.tsv"), sep="\t", index=False)

    excluded = pool[~pool.included]
    exc_tbl = tidy(excluded)
    exc_tbl["exclusion_reason"] = excluded.exclusion_reason.values
    exc_tbl.to_csv(out("_excluded_records.tsv"), sep="\t", index=False)

    species = (curated.groupby(["organism", "taxid"], as_index=False)
               .agg(n_records=("accession", "size"), max_length_bp=("length", "max"),
                    virus_family=("virus_family", "first"), tier=("tier", "first"),
                    probable_host_class=("probable_host_class", "first"),
                    hosts=("host_scientific_name", joined),
                    countries=("country", lambda s: joined(x.split(":")[0] for x in s)),
                    years=("collection_date", lambda s: joined(
                        re.search(r"\d{4}", x).group() for x in s if re.search(r"\d{4}", x))),
                    pubmed_ids=("pubmed", lambda s: joined(p for x in s for p in str(x).split(";"))))
               .rename(columns={"organism": "virus_name"}))
    species = (species.merge(ref_tbl[["virus_name", "accession", "genome_status", "length_bp"]]
                            .rename(columns={"accession": "representative_accession",
                                             "genome_status": "representative_genome_status",
                                             "length_bp": "representative_length_bp"}),
                            on="virus_name", how="left")
                      .drop_duplicates(["virus_name", "taxid"])
                      .sort_values(["tier", "virus_family", "virus_name"]))
    species.to_csv(out("_species_summary.tsv"), sep="\t", index=False)

    t1 = reps[reps.tier == TIER1]
    usable = ["complete genome", "coding-complete", "near-complete"]
    host_counts = (curated.host_scientific_name.replace("", "not resolved to species")
                   .value_counts().to_frame("records"))
    status_counts = reps.genome_status.value_counts().to_frame("representatives")
    exc_counts = excluded.exclusion_reason.str.slice(0, 70).value_counts().to_frame("records")
    meta = json.load(open(args.meta, encoding="utf-8"))

    n_audit, screen_txt = 0, "Augmentation was disabled, so no species-level references were considered."
    if args.audit and os.path.exists(args.audit):
        aud = pd.read_csv(args.audit, sep="\t")
        n_audit = len(aud)
        if n_audit:
            npass = int((aud.gate == "pass").sum())
            cont = (pd.to_numeric(reps.loc[~reps.from_host_taxon,
                                            "kmer_containment_to_host_record"], errors="coerce")
                    .dropna() if "kmer_containment_to_host_record" in reps.columns
                    else pd.Series(dtype=float))
            screen_txt = (
                "%d candidate reference genomes were screened for similarity to the host-derived "
                "records of their own taxid: %d cleared the %.2f containment threshold and %d were "
                "rejected, of which %d shared no %d-mer with any host-derived record and %d carried a "
                "contradictory /segment label. Where a candidate was rejected the host-derived record "
                "was kept, so those taxa are represented by what was actually sequenced from the host. "
                "The imported references that were delivered span %.2f-%.2f containment (median %.2f)."
                % (n_audit, npass, min_cont, n_audit - npass,
                   int((aud.kmer_containment == 0).sum()), kmer_size,
                   int((aud.segment_check == "conflict").sum()),
                   cont.min() if len(cont) else 0.0, cont.max() if len(cont) else 0.0,
                   cont.median() if len(cont) else 0.0))
        else:
            screen_txt = "No candidate reference genomes were available to screen."

    report = Template(textwrap.dedent(REPORT)).safe_substitute(
        screen_txt=screen_txt, n_audit=n_audit,
        dataset=args.dataset, date=datetime.date.today().isoformat(),
        query=open(args.query, encoding="utf-8").read().strip(),
        hits=meta["hit_count"], parsed=len(pool),
        curated=len(curated), excluded=len(excluded), reps=len(reps),
        n_t1=len(t1), n_t2=int((reps.tier == TIER2).sum()),
        n_nontarget=int((~reps.from_host_taxon).sum()),
        n_nontarget_taxa=int(reps.loc[~reps.from_host_taxon, "taxid"].nunique()),
        n_target=int(reps.from_host_taxon.sum()),
        t1_usable=int(t1.genome_status.isin(usable).sum()), n_species=len(species),
        terms_file=args.terms, accept_file=args.host_accept, tier2_file=args.tier2,
        taxroot=cfg["entrez"]["virus_taxid_root"],
        near_min=cfg["curation"]["completeness"]["near_complete_min_pct"],
        frag_min=cfg["curation"]["completeness"]["substantial_fragment_min_pct"],
        host_table=host_counts.to_markdown(), status_table=status_counts.to_markdown(),
        exc_table=exc_counts.to_markdown() if len(exc_counts) else "_no records rejected_")
    with open(os.path.join(args.outdir, "README_" + args.dataset + "_refdb.md"),
              "w", encoding="utf-8", newline="\n") as fh:
        fh.write(report.lstrip("\n"))

    print("representatives: %d (%d tier1 / %d tier2) | species rows: %d"
          % (len(reps), len(t1), int((reps.tier == TIER2).sum()), len(species)))


if __name__ == "__main__":
    main()
