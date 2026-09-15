# host-virus-refdb

A Snakemake workflow that builds a curated, reproducible **virus reference database for a host
taxon** from NCBI GenBank/RefSeq. Nothing about the target taxon is hard-coded: three tab-separated
input files define the scope, so the same workflow builds a crayfish database, a shrimp database, or
a database for any other host group.

It was written for crayfish virome work (the shipped `config/` reproduces that build), but the only
taxon-specific knowledge lives in the input files.

## What it produces

| Output | Description |
|---|---|
| `<dataset>_refdb_tier1.fasta` | Representative genomes plausibly infecting the host taxon or its relatives |
| `<dataset>_refdb_tier2.fasta` | Viruses of plants, fungi, algae and gut bacteria recovered from the same samples — a decoy set for read screening |
| `<dataset>_refdb_all.fasta` | Both tiers |
| `<dataset>_reference_genomes.tsv` | Metadata for every sequence in the FASTA files (30 columns) |
| `<dataset>_all_records.tsv` | Every curated host-taxon record, before representative selection |
| `<dataset>_species_summary.tsv` | One row per virus taxon: hosts, countries, years, PubMed IDs, representative |
| `<dataset>_excluded_records.tsv` | Retrieved but rejected, each with its reason — the audit trail |
| `<dataset>_refdb_overview.png` | Composition figure |
| `README_<dataset>_refdb.md` | Build report: exact query, counts, method, caveats |

## Quick start

```bash
git clone <your-fork-url> && cd host-virus-refdb
conda env create -f environment.yaml && conda activate host-virus-refdb

export NCBI_API_KEY=...          # optional but recommended: 10 requests/s instead of 3
snakemake --cores 1              # add --use-conda to let Snakemake provision each rule's env
```

Sanity-check the install before the full build — this takes about 20 seconds and touches
NCBI for roughly five records:

```bash
snakemake --cores 1 --configfile config/config.test.yaml
```

Outputs land in `results/`. Rerunning is cheap: downloaded GenBank XML is cached, so changing a
curation rule re-runs only the curation and build steps.

## Retargeting it to another host taxon

Edit three files (all comment-documented, all tab-separated):

**`config/search_terms.tsv`** — the host taxon's genera, families and vernacular names.
The `use` column controls how each term is applied:

| `use` | Effect |
|---|---|
| `query` | included in the Entrez query only |
| `name` | used only to recognise host-named viruses, so they stay in tier 1 — e.g. `shrimp`, `decapod` |
| `both` | both of the above |

**`config/host_accept.tsv`** — which hosts count. Every distinct GenBank `/host` string is resolved
against NCBI Taxonomy and kept only if its lineage contains one of your `clade` rows, or the raw
string matches a `vernacular` pattern. This is what stops keyword false positives: in the crayfish
build, "marron" (a *Cherax* vernacular) also matches *Castanea sativa* cv. "Marron" chestnut
cultivars, and the taxonomy check discards those records.

**`config/tier2_groups.tsv`** — viral groups that are *not* candidate pathogens of your host.
`kind=environmental` sends a group to tier 2; `kind=phage` rejects it outright. Any taxon name that
appears in a GenBank lineage works (family, order, class).

Then set `outputs.dataset_name` in `config/config.yaml` and run. A worked second example is in
`config/examples/penaeid_shrimp/`.

## Curation logic

1. **Retrieval** — one Entrez query, full GenBank XML, so host/country/date come from source
   qualifiers rather than definition lines.
2. **Host curation** — `/host` strings resolved against NCBI Taxonomy against your accepted clades.
   Host-less records are kept if the virus name carries a host term, or the isolation source names
   the host taxon (symbiont sources are rejected — the host is the symbiont, not your taxon).
3. **Tiering** — from `tier2_groups.tsv`. A virus whose *name* carries a host term stays in tier 1
   even when its family is listed as environmental. Tier is a triage annotation inferred from
   taxonomy, not an experimental host determination.
4. **Completeness** — GenBank definition first (`complete genome`, `coding-complete`), otherwise a
   percentage of the longest definition-complete member of the same family. Percent-of-longest
   *within a taxon* is not usable as the primary criterion, because novel virome taxa usually have
   one record and would score 100% by construction.
5. **Representatives** — one per taxon per segment: complete > coding-complete > near-complete >
   fragment, then RefSeq, then host-derived, then longest. Placeholder taxa (`Picornavirales sp.`
   and friends) are kept per accession, never collapsed — one such taxid is shared by unrelated
   viruses.
6. **Species-level references** — for a virus reported from your taxon but represented only by
   fragments, the complete genome deposited from another host is added and flagged
   `from_host_taxon = False`, so you can filter it out.

## Workflow structure

```
config/          config.yaml + the three input files (+ examples/)
workflow/
  Snakefile      targets
  rules/         retrieve.smk, curate.smk, build.smk
  scripts/       one script per rule + refdb_entrez/refdb_io/refdb_gb libraries
  envs/          conda environment
results/         outputs (git-ignored)
```

Useful invocations:

```bash
snakemake --cores 1 -n                                  # dry run
snakemake --cores 1 --configfile config/config.test.yaml  # small smoke test (one genus)
snakemake --cores 1 results/tables/records_curated.tsv.gz  # stop after curation
snakemake --report report.html                          # provenance report
```

## Requirements

Python 3.9+, pandas, matplotlib, tabulate, snakemake >= 7.32. Network access to
`eutils.ncbi.nlm.nih.gov`. No NCBI account is required; an API key only raises the rate limit
from 3 to 10 requests/s.

## Running on Linux or an HPC server

The workflow is plain Python and Snakemake with no compiled or platform-specific dependencies,
and it is exercised on `ubuntu-latest` by the CI workflow in `.github/workflows/ci.yml`, which
runs both a dry run and a real end-to-end build on every push.

- **Interpreter.** Each step runs as `{python_exe} workflow/scripts/<step>.py`. Inside a conda
  env, the default `python` is correct. On a bare server where only `python3` exists, set
  `runtime.python_exe: python3` in `config/config.yaml` — no other change is needed.
- **Headless.** The figure step forces matplotlib's Agg backend, so no display or X11
  forwarding is required; `matplotlib-base` is sufficient.
- **Cores.** Retrieval is serialised by NCBI's rate limit, not by CPU, so `--cores 1` is the
  right setting; the build steps are single-threaded pandas. Memory stays well under 1 GB for a
  few thousand records.
- **Cluster submission.** Nothing in the workflow assumes local execution — add a Snakemake
  executor plugin if you want the steps submitted as jobs (`snakemake --executor slurm
  --default-resources ...`). For a database of this size the whole build takes a few minutes on
  one core, so a login-node run inside `tmux` is usually simpler than a job submission.
- **Outbound network.** Only `https://eutils.ncbi.nlm.nih.gov` is contacted. If the server
  reaches the internet through a proxy, the standard `https_proxy` / `HTTPS_PROXY` environment
  variables are honoured, since the client uses `urllib` from the standard library.
- **Resumability.** Downloaded GenBank XML is cached under `results/raw/`, so an interrupted or
  rate-limited run resumes where it stopped rather than re-fetching.

## Citation

If a database built with this workflow supports a publication, cite NCBI GenBank and the primary
studies behind the sequences you use — `<dataset>_species_summary.tsv` carries the PubMed IDs.
