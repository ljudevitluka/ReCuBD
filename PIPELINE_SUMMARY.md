# ViroCray Host Virus Reference Database Pipeline

## Overview

This Snakemake workflow builds a curated virus reference database for a chosen host taxon. The current configuration targets crayfish, but the workflow can be retargeted by editing the input files in `config/`.

The workflow is:

```text
configuration
    -> Entrez search
    -> GenBank XML retrieval
    -> record parsing
    -> host curation
    -> virus classification
    -> external-reference augmentation
    -> similarity screening
    -> representative selection
    -> sequence download
    -> final database files
```

## 1. Configuration

The main workflow is `workflow/Snakefile`. It loads `config/config.yaml`, which defines:

- input files
- NCBI Entrez settings
- accepted host logic
- completeness thresholds
- similarity-screening thresholds
- augmentation behavior
- output names and location

The workflow first writes a frozen configuration to:

```text
results/effective_config.json
```

This makes each run reproducible and allows Snakemake to detect configuration changes.

## 2. Build the NCBI query

The `build_query` rule runs `workflow/scripts/build_query.py`.

It:

1. Reads `config/search_terms.tsv`.
2. Separates search terms from host-name recognition terms.
3. Combines the search terms into one NCBI Entrez query.
4. Restricts the query to the virus taxonomy root, taxid `10239`.
5. Sends the query to NCBI.
6. Saves the matching NCBI UIDs and search metadata.

Outputs:

```text
results/query/entrez_query.txt
results/query/uids.txt
results/query/search_meta.json
```

## 3. Download GenBank records

The `fetch_genbank` rule runs `workflow/scripts/fetch_genbank.py`.

It downloads full GenBank XML records from NCBI in batches. Full records are needed to inspect qualifiers such as `/host`, `/isolation_source`, country, collection date, and segment information.

Records are stored in:

```text
results/raw/genbank/
```

The download is resumable. Existing XML files larger than 500 bytes are skipped.

## 4. Parse GenBank XML

The `parse_genbank` rule runs `workflow/scripts/parse_genbank.py`.

It:

1. Reads all downloaded XML files.
2. Extracts accession, taxonomy, host, isolation source, collection information, sequence length, and related metadata.
3. Converts the records into a pandas table.
4. Removes duplicate accessions.

Output:

```text
results/tables/records_raw.tsv.gz
```

## 5. Curate host associations

The `curate_hosts` rule uses:

- `config/host_accept.tsv`
- `config/search_terms.tsv`
- `records_raw.tsv.gz`

`workflow/scripts/curate_hosts.py` resolves host names against NCBI Taxonomy and determines whether each record is associated with the target host taxon.

A record may be accepted through:

- a valid `/host` qualifier
- a virus name containing an accepted host term
- an `/isolation_source` containing an accepted host term

Symbiont and parasite sources can be rejected so that a virus of a crayfish-associated organism is not incorrectly treated as a crayfish virus.

Outputs:

```text
results/tables/records_curated.tsv.gz
results/tables/host_taxonomy.json
```

## 6. Classify records into tiers

The `annotate` rule runs `workflow/scripts/annotate.py` and uses `config/tier2_groups.tsv`.

It determines:

- whether a record is a phage
- whether it belongs to a tier-2 environmental group
- whether the virus name is associated with the target host
- the probable host class
- whether the record is included or excluded

Typical tiers are:

```text
tier1_host_and_related
tier2_environmental_diet_microbiome
```

Phages are rejected when configured because their actual hosts are bacteria rather than the target animal.

Output:

```text
results/tables/records_annotated.tsv.gz
```

## 7. Retrieve additional references

The `augment_references` rule searches for complete genomes deposited from hosts outside the target taxon.

This helps when the target host has only a fragment of a virus but a complete genome exists from another host.

The rule:

1. Selects included tier-1 host-derived records.
2. Extracts their viral taxonomic IDs.
3. Searches NCBI for complete genomes and RefSeq records.
4. Processes taxonomic IDs in groups.
5. Downloads candidate XML records.
6. Removes accessions already present in the host-derived dataset.

Configuration:

```yaml
augmentation:
  enabled: true
  title_terms:
    - complete genome
    - complete sequence
  include_refseq: true
  max_per_chunk: 400
  taxids_per_query: 60
```

Outputs:

```text
results/raw/augment/
results/tables/augment_records.tsv.gz
```

These records are candidates only. They are not automatically accepted into the final database.

## 8. Screen augmented references by similarity

The updated workflow adds the `screen_augmented` rule, which runs `workflow/scripts/screen_augmented.py`.

An NCBI taxonomic ID alone does not guarantee that two records represent the same biological entity. A single taxonomic ID may contain divergent sequences or different genome segments. Therefore, each augmented candidate is screened before it can compete with a host-derived record.

The screening process:

1. Keeps candidates whose tax IDs also occur in host-derived records.
2. Ranks candidates by completeness, RefSeq status, and sequence length.
3. Limits the number screened per taxon and segment.
4. Selects host-derived records for comparison.
5. Downloads missing sequences into a local FASTA cache.
6. Calculates canonical k-mer sets.
7. Compares each candidate with host-derived records having the same taxon.
8. Checks segment-label compatibility.
9. Keeps only candidates that pass the similarity threshold.

Default settings:

```yaml
similarity:
  kmer_size: 21
  min_containment: 0.5
  enforce_segment_match: true
  max_candidates_per_taxon: 3
  max_host_records_per_taxon: 5
```

### K-mer containment

Containment is calculated as:

$$
C = \frac{\text{host-record k-mers also found in candidate}}{\text{host-record k-mers}}
$$

The script estimates ANI using:

$$
\widehat{\mathrm{ANI}} = 1 + \frac{\ln(C)}{k}
$$

At `k=21`, the default containment threshold of `0.5` corresponds approximately to 96.7% estimated ANI. This estimate is mainly a screening gate and should not be treated as a replacement for a full ANI analysis.

A candidate fails when:

- its sequence cannot be retrieved
- no comparable host-derived sequence exists
- its segment label conflicts with all host-derived records
- its containment is below `min_containment`

Outputs:

```text
results/tables/augment_screened.tsv.gz
results/tables/augment_similarity_audit.tsv
results/seq/screened_candidates.fasta
```

The audit table records every candidate, including its best host accession, containment, estimated ANI, segment check, pass/fail status, and rejection reason.

## 9. Select representative records

The `select_representatives` rule runs `workflow/scripts/select_representatives.py`.

It combines:

- annotated host-derived records
- similarity-screened augmented records

It assigns genome completeness categories:

1. complete genome
2. coding-complete sequence
3. near-complete sequence
4. substantial fragment
5. gene fragment
6. unscored record

For incomplete sequences, length is compared with a family-calibrated reference length. The current thresholds are:

```yaml
near_complete_min_pct: 85
substantial_fragment_min_pct: 50
```

Representative selection favors better genome status, RefSeq records, host-derived records, and longer sequences. It selects one representative per taxonomic group and segment.

Outputs:

```text
results/tables/records_pooled.tsv.gz
results/tables/representatives.tsv.gz
```

## 10. Download representative sequences

The `fetch_sequences` rule runs `workflow/scripts/fetch_sequences.py`.

It:

1. Reads the selected representatives.
2. Uses accession versions where available.
3. Downloads sequences from NCBI in batches.
4. Writes a FASTA file.
5. Verifies that the number of retrieved sequences matches the number expected.

Output:

```text
results/seq/representatives_raw.fasta
```

## 11. Write final database files

The `write_outputs` rule runs `workflow/scripts/write_outputs.py`.

With the current dataset name, the final outputs include:

```text
results/crayfish_virus_refdb_tier1.fasta
results/crayfish_virus_refdb_tier2.fasta
results/crayfish_virus_refdb_all.fasta
results/crayfish_virus_reference_genomes.tsv
results/crayfish_virus_species_summary.tsv
results/crayfish_virus_all_records.tsv
results/crayfish_virus_excluded_records.tsv
results/README_crayfish_virus_refdb.md
results/tables/augment_similarity_audit.tsv
```

When plotting is enabled, the workflow also creates:

```text
results/crayfish_virus_refdb_overview.png
```

## What changed in the updated version

Previously, external complete references were added directly to representative selection:

```text
host-derived records
    + external complete references
    -> representative selection
```

The updated workflow screens those references first:

```text
host-derived records
    + external complete references
    -> k-mer containment and segment screening
    -> representative selection
```

This prevents an external reference from replacing a host-derived sequence solely because both records share an NCBI taxonomic ID. The candidate must resemble the host-derived sequence and pass the segment-consistency check.

The process is auditable because every augmented candidate is listed in:

```text
results/tables/augment_similarity_audit.tsv
```

## Running the workflow

From the repository root:

```bash
conda activate host-virus-refdb
snakemake --cores 1
```

For a dry run:

```bash
snakemake --cores 1 -n
```

For a small test configuration:

```bash
snakemake --cores 1 --configfile config/config.test.yaml
```

The workflow requires Snakemake `>=7.32`, Python `>=3.9`, pandas, matplotlib, tabulate, and network access to NCBI Entrez. Set `NCBI_API_KEY` to increase the NCBI request rate limit.

The workflow is resumable. Existing downloaded XML and cached sequence files are reused after an interruption or failure.
