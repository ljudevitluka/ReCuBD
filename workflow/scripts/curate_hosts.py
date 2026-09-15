"""Resolve every /host string against NCBI Taxonomy and keep host-taxon records.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import re

import refdb_entrez
import refdb_gb
import refdb_io


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--host-accept", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-host-tax", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.config_json, encoding="utf-8"))
    curation = cfg["curation"]
    accept = refdb_io.load_host_accept(args.host_accept)
    terms = refdb_io.load_search_terms(args.terms)
    name_re = refdb_io.term_regex(terms["name"])
    vern_re = refdb_io.term_regex(accept["vernacular"])
    clade_re = (re.compile(r"\b(" + "|".join(re.escape(c) for c in accept["clades"]) + r")\b")
                if accept["clades"] else re.compile(r"(?!x)x"))

    df = refdb_gb.add_viral_taxonomy(refdb_io.read_table(args.table))
    client = refdb_entrez.client_from_params(cfg["entrez"])

    def resolve(name):
        q = re.sub(r"\s*\(.*?\)", " ", name)
        q = re.sub(r"\s+(cv\.|var\.|subsp\.|strain)\s+.*$", "", q, flags=re.I).strip()
        if not q:
            return None
        hits = client.xml("esearch", db="taxonomy", term='"' + q + '"[All Names]', retmax=1).findall(".//Id")
        if not hits:
            hits = client.xml("esearch", db="taxonomy", term=q, retmax=1).findall(".//Id")
        if not hits:
            return None
        tax = client.xml("efetch", db="taxonomy", id=hits[0].text)
        return {"query": q, "taxid": hits[0].text,
                "sci_name": tax.findtext(".//Taxon/ScientificName") or "",
                "rank": tax.findtext(".//Taxon/Rank") or "",
                "lineage": tax.findtext(".//Taxon/Lineage") or ""}

    host_tax = {}
    for h in sorted({x.strip() for x in df.host if x.strip()}):
        host_tax[h] = resolve(h)
        print("host: %s -> %s" % (h, host_tax[h]["sci_name"] if host_tax[h] else "unresolved"), flush=True)

    def in_taxon(raw, info):
        if info and clade_re.search(info["lineage"] + "; " + info["sci_name"]):
            return True
        return bool(vern_re.search(raw))

    def classify(row):
        raw = row.host.strip()
        info = host_tax.get(raw)
        sci = info["sci_name"] if info else ""
        taxid = info["taxid"] if info else ""
        if raw:
            if in_taxon(raw, info):
                return "host qualifier", sci, taxid, ""
            return "none", sci, taxid, "non-target host: " + (sci or raw)
        if curation.get("accept_virus_name_match", True) and name_re.search(row.organism):
            return "virus name", "", "", ""
        if curation.get("accept_isolation_source", True) and vern_re.search(row.isolation_source):
            if curation.get("reject_symbiont_sources", True) and \
                    re.search(r"symbiont|epibiont|parasite of", row.isolation_source, re.I):
                return ("isolation source", "", "",
                        "host is a symbiont of the target taxon, not the taxon itself: " + row.isolation_source)
            return "isolation source", "", "", ""
        return "none", "", "", "no host-taxon evidence (isolation_source: " + (row.isolation_source or "n/a") + ")"

    res = df.apply(classify, axis=1, result_type="expand")
    res.columns = ["host_evidence", "host_scientific_name", "host_taxid", "host_reject_reason"]
    df = df.join(res)
    df["from_host_taxon"] = df.host_evidence != "none"
    df["refseq"] = df.accession.str.match(r"^(NC_|AC_|NG_)") | df.db_source.str.upper().str.contains("REFSEQ")
    df["placeholder_taxon"] = df.organism.str.contains(r"\bsp\.$|^unclassified", case=False, regex=True)

    refdb_io.write_table(df, args.out)
    with open(args.out_host_tax, "w", encoding="utf-8") as fh:
        json.dump(host_tax, fh, indent=2, sort_keys=True)
    print("host-taxon records: %d of %d" % (int(df.from_host_taxon.sum()), len(df)))


if __name__ == "__main__":
    main()
