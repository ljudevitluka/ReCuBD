# Retrieval from NCBI Entrez. Separate rules so that a curation change does not
# re-download anything.


rule build_query:
    """OR-join the host search terms into one Entrez query and resolve it to UIDs."""
    input:
        cfg=CONFIG_JSON,
        terms=IN["search_terms"],
    output:
        query=f"{RES}/query/entrez_query.txt",
        uids=f"{RES}/query/uids.txt",
        meta=f"{RES}/query/search_meta.json",
    log:
        f"{RES}/logs/build_query.log",
    params:
        py=PY,
        script=SCRIPTS / "build_query.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --terms {input.terms} '
        "--out-query {output.query} --out-uids {output.uids} --out-meta {output.meta} > {log} 2>&1"


rule fetch_genbank:
    """Download full GenBank records as XML, so source qualifiers are available."""
    input:
        cfg=CONFIG_JSON,
        uids=rules.build_query.output.uids,
    output:
        xmldir=directory(f"{RES}/raw/genbank"),
    log:
        f"{RES}/logs/fetch_genbank.log",
    params:
        py=PY,
        script=SCRIPTS / "fetch_genbank.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --uids {input.uids} '
        "--outdir {output.xmldir} > {log} 2>&1"


rule parse_genbank:
    input:
        xmldir=rules.fetch_genbank.output.xmldir,
    output:
        table=f"{RES}/tables/records_raw.tsv.gz",
    log:
        f"{RES}/logs/parse_genbank.log",
    params:
        py=PY,
        script=SCRIPTS / "parse_genbank.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --xmldir {input.xmldir} --out {output.table} > {log} 2>&1'


rule augment_references:
    """Optional: complete genomes of tier-1 taxa deposited from other hosts."""
    input:
        cfg=CONFIG_JSON,
        annotated=f"{RES}/tables/records_annotated.tsv.gz",
    output:
        xmldir=directory(f"{RES}/raw/augment"),
        table=f"{RES}/tables/augment_records.tsv.gz",
    log:
        f"{RES}/logs/augment_references.log",
    params:
        py=PY,
        script=SCRIPTS / "augment_references.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --annotated {input.annotated} '
        "--outdir {output.xmldir} --out {output.table} > {log} 2>&1"


rule fetch_sequences:
    input:
        cfg=CONFIG_JSON,
        reps=f"{RES}/tables/representatives.tsv.gz",
    output:
        fasta=f"{RES}/seq/representatives_raw.fasta",
    log:
        f"{RES}/logs/fetch_sequences.log",
    params:
        py=PY,
        script=SCRIPTS / "fetch_sequences.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --reps {input.reps} '
        "--out {output.fasta} > {log} 2>&1"
