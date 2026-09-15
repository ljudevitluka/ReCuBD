# Host curation and annotation -- the taxon-specific logic, driven entirely by
# the three input files.


rule curate_hosts:
    """Resolve every /host string against NCBI Taxonomy; keep host-taxon records."""
    input:
        cfg=CONFIG_JSON,
        table=rules.parse_genbank.output.table,
        host_accept=IN["host_accept"],
        terms=IN["search_terms"],
    output:
        table=f"{RES}/tables/records_curated.tsv.gz",
        host_tax=f"{RES}/tables/host_taxonomy.json",
    log:
        f"{RES}/logs/curate_hosts.log",
    params:
        py=PY,
        script=SCRIPTS / "curate_hosts.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --table {input.table} '
        "--host-accept {input.host_accept} --terms {input.terms} "
        "--out {output.table} --out-host-tax {output.host_tax} > {log} 2>&1"


rule annotate:
    """Assign tier, probable host class and phage rejection from tier2_groups.tsv."""
    input:
        cfg=CONFIG_JSON,
        table=rules.curate_hosts.output.table,
        tier2=IN["tier2_groups"],
        terms=IN["search_terms"],
    output:
        table=f"{RES}/tables/records_annotated.tsv.gz",
    log:
        f"{RES}/logs/annotate.log",
    params:
        py=PY,
        script=SCRIPTS / "annotate.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --table {input.table} '
        "--tier2 {input.tier2} --terms {input.terms} --out {output.table} > {log} 2>&1"


rule screen_augmented:
    """Require a candidate reference genome to resemble the host-derived records of
    its own taxid before it may compete for the representative slot."""
    input:
        cfg=CONFIG_JSON,
        annotated=rules.annotate.output.table,
        augment=rules.augment_references.output.table,
    output:
        passed=f"{RES}/tables/augment_screened.tsv.gz",
        audit=f"{RES}/tables/augment_similarity_audit.tsv",
    log:
        f"{RES}/logs/screen_augmented.log",
    params:
        py=PY,
        script=SCRIPTS / "screen_augmented.py",
        cache=f"{RES}/seq/screened_candidates.fasta",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --annotated {input.annotated} '
        "--augment {input.augment} --cache {params.cache} "
        "--out-passed {output.passed} --out-audit {output.audit} > {log} 2>&1"
