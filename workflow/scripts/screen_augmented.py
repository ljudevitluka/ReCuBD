"""Screen candidate reference genomes for similarity to the host-derived records.

A candidate pulled in by augment_references shares only an NCBI taxid with the
record it would stand in for. Taxids are assigned by submitters, so for novel
virome taxa one taxid can bin sequences that are 75-90% identical or less, and a
candidate can be a different virus or a different genome segment. This step
requires a candidate to resemble the host-derived material before it is allowed
to compete for the representative slot; candidates that fail are dropped and the
host-derived record is kept.

Similarity is canonical k-mer containment: the fraction of the host record's
k-mers present in the candidate. Estimated ANI is the Mash containment estimate
1 + ln(containment)/k. At k=21, containment 0.9 ~ 99.5% ANI, 0.5 ~ 96.7%,
0.2 ~ 92.3%, 0.05 ~ 85.8%.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import refdb_entrez
import refdb_io

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical_kmers(seq, k):
    out = set()
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" in kmer:
            continue
        out.add(min(kmer, kmer.translate(COMPLEMENT)[::-1]))
    return out


def normalise_segment(label):
    """'S4', 'segment 4', 'RNA4' -> '4'; '' -> '' (unlabelled)."""
    up = re.sub(r"[^A-Z0-9]", "", str(label).upper())
    up = re.sub(r"^(SEGMENT|SEG|RNA|DNA)", "", up)
    up = re.sub(r"^S(?=\d)", "", up)
    return up


def est_ani_pct(containment, k):
    if containment <= 0:
        return float("nan")
    return round(100 * (1 + math.log(containment) / k), 1)


def definition_rank(definition):
    d = str(definition).lower()
    if re.search(r"complete genome|complete sequence|genome assembly", d):
        return 0
    if re.search(r"complete cds|polyprotein gene, complete", d):
        return 1
    return 2


def read_fasta(path):
    seqs, key = {}, None
    if not os.path.exists(path):
        return seqs
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            key = line[1:].split()[0].split("|")[0]
            seqs[key] = []
        elif key and line.strip():
            seqs[key].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def fetch_sequences(client, ids, cache_path, batch=50):
    """Fetch any ids missing from the on-disk cache, then return the cache."""
    cached = read_fasta(cache_path)
    missing = [i for i in ids if i not in cached and i.split(".")[0] not in
               {c.split(".")[0] for c in cached}]
    if missing:
        os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
        with open(cache_path, "ab") as fh:
            for start in range(0, len(missing), batch):
                fh.write(client.efetch("nuccore", missing[start:start + batch],
                                       rettype="fasta", retmode="text"))
                print("  fetched %d/%d sequences" % (min(start + batch, len(missing)), len(missing)),
                      flush=True)
        cached = read_fasta(cache_path)
    return cached


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--annotated", required=True, help="host-derived records (annotate step)")
    ap.add_argument("--augment", required=True, help="candidate genomes (augment_references step)")
    ap.add_argument("--cache", required=True, help="FASTA cache for screened sequences")
    ap.add_argument("--out-passed", required=True, help="candidates cleared to compete")
    ap.add_argument("--out-audit", required=True, help="every screened candidate, with scores")
    args = ap.parse_args()

    cfg = json.load(open(args.config_json, encoding="utf-8"))
    sim_cfg = cfg.get("similarity", {})
    k = int(sim_cfg.get("kmer_size", 21))
    min_cont = float(sim_cfg.get("min_containment", 0.5))
    max_cand = int(sim_cfg.get("max_candidates_per_taxon", 3))
    max_host = int(sim_cfg.get("max_host_records_per_taxon", 5))
    enforce_segment = bool(sim_cfg.get("enforce_segment_match", True))

    host = refdb_io.read_table(args.annotated)
    cand = refdb_io.read_table(args.augment)
    audit_cols = ["accession", "organism", "taxid", "length", "segment", "definition",
                  "best_host_accession", "best_host_length", "kmer_containment",
                  "est_ani_pct", "segment_check", "gate", "gate_reason"]

    if not len(cand):
        refdb_io.write_table(cand, args.out_passed)
        refdb_io.write_table(pd.DataFrame(columns=audit_cols), args.out_audit)
        print("no candidates to screen")
        return

    host = host[host.included & host.from_host_taxon]
    cand = cand[cand.taxid.isin(set(host.taxid))].copy()

    # Bound the work: only screen candidates that could plausibly win a slot. The
    # ranking must match the preference order used by select_representatives --
    # completeness, then RefSeq, then length -- or a RefSeq genome can be excluded
    # from screening by a longer non-RefSeq isolate and lose the slot for that reason
    # alone.
    cand["def_rank"] = cand.definition.apply(definition_rank)
    cand["segment_key"] = cand.segment.apply(normalise_segment)
    cand["is_refseq"] = (cand.accession.str.match(r"^(NC_|AC_|NG_)", na=False)
                         | cand.db_source.str.upper().str.contains("REFSEQ", na=False))
    cand = (cand.sort_values(["def_rank", "is_refseq", "length"], ascending=[True, False, False])
                .groupby(["taxid", "segment_key"], group_keys=False).head(max_cand))

    # only taxa that actually have candidates need their sequences fetched
    host_pick = (host[host.taxid.isin(set(cand.taxid))]
                 .sort_values("length", ascending=False)
                 .groupby("taxid", group_keys=False).head(max_host).copy())
    host_pick["segment_key"] = host_pick.segment.apply(normalise_segment)

    client = refdb_entrez.client_from_params(cfg["entrez"])
    ids = [(v or a) for v, a in zip(cand.version, cand.accession)] + \
          [(v or a) for v, a in zip(host_pick.version, host_pick.accession)]
    print("screening %d candidates against %d host-derived records (%d taxa)"
          % (len(cand), len(host_pick), cand.taxid.nunique()), flush=True)
    seqs = fetch_sequences(client, ids, args.cache)

    def lookup(row):
        key = row.version or row.accession
        return seqs.get(key) or seqs.get(key.split(".")[0]) or next(
            (s for kk, s in seqs.items() if kk.split(".")[0] == str(row.accession)), None)

    host_kmers = {}
    for r in host_pick.itertuples():
        s = lookup(r)
        if s and len(s) >= k:
            host_kmers[r.accession] = (canonical_kmers(s, k), r.taxid, r.segment_key, len(s))

    rows = []
    for c in cand.itertuples():
        cseq = lookup(c)
        if not cseq or len(cseq) < k:
            rows.append(dict(accession=c.accession, organism=c.organism, taxid=c.taxid,
                             length=c.length, segment=c.segment, definition=c.definition,
                             best_host_accession="", best_host_length=0, kmer_containment=0.0,
                             est_ani_pct=float("nan"), segment_check="n/a", gate="fail",
                             gate_reason="no sequence retrieved for candidate"))
            continue
        ck = canonical_kmers(cseq, k)
        comparable, conflicts = [], 0
        for acc, (hk, taxid, hseg, hlen) in host_kmers.items():
            if taxid != c.taxid:
                continue
            if enforce_segment and hseg and c.segment_key and hseg != c.segment_key:
                conflicts += 1
                continue
            comparable.append((len(hk & ck) / len(hk) if hk else 0.0, acc, hlen, hseg))
        if not comparable:
            reason = ("segment label conflicts with every host-derived record"
                      if conflicts else "no comparable host-derived sequence")
            rows.append(dict(accession=c.accession, organism=c.organism, taxid=c.taxid,
                             length=c.length, segment=c.segment, definition=c.definition,
                             best_host_accession="", best_host_length=0, kmer_containment=0.0,
                             est_ani_pct=float("nan"),
                             segment_check="conflict" if conflicts else "n/a",
                             gate="fail", gate_reason=reason))
            continue
        cont, acc, hlen, hseg = max(comparable)
        passed = cont >= min_cont
        rows.append(dict(
            accession=c.accession, organism=c.organism, taxid=c.taxid, length=c.length,
            segment=c.segment, definition=c.definition, best_host_accession=acc,
            best_host_length=hlen, kmer_containment=round(cont, 4), est_ani_pct=est_ani_pct(cont, k),
            segment_check=("match" if (hseg and c.segment_key) else "unlabelled"),
            gate="pass" if passed else "fail",
            gate_reason="" if passed else (
                "no %d-mer shared with any host-derived record" % k if cont == 0 else
                "%d-mer containment %.3f below threshold %.2f (est. ANI %.1f%%)"
                % (k, cont, min_cont, est_ani_pct(cont, k)))))

    audit = pd.DataFrame(rows, columns=audit_cols)
    refdb_io.write_table(audit, args.out_audit)

    passed_acc = set(audit.loc[audit.gate == "pass", "accession"])
    out = refdb_io.read_table(args.augment)
    out = out[out.accession.isin(passed_acc)]
    keep = audit[audit.gate == "pass"][["accession", "kmer_containment", "est_ani_pct",
                                        "best_host_accession"]]
    out = out.merge(keep, on="accession", how="left").rename(columns={
        "kmer_containment": "kmer_containment_to_host_record",
        "est_ani_pct": "est_ani_to_host_record_pct",
        "best_host_accession": "compared_with_host_accession"})
    refdb_io.write_table(out, args.out_passed)

    print("passed: %d of %d candidates (%d taxa retain an imported reference)"
          % (len(out), len(audit), out.taxid.nunique()))
    if (audit.gate == "fail").any():
        worst = audit[audit.gate == "fail"].nlargest(5, "kmer_containment")
        print("rejected %d; highest-scoring rejections:" % int((audit.gate == "fail").sum()))
        for r in worst.itertuples():
            print("  %-42s %s  containment %.3f" % (r.organism[:42], r.accession, r.kmer_containment))


if __name__ == "__main__":
    main()
