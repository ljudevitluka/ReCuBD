# Representative selection and deliverables.


rule select_representatives:
    """Score completeness against a family-calibrated reference length, then keep
    one sequence per taxon per segment."""
    input:
        cfg=CONFIG_JSON,
        annotated=rules.annotate.output.table,
        augment=rules.screen_augmented.output.passed,
        tier2=IN["tier2_groups"],
        terms=IN["search_terms"],
    output:
        pooled=f"{RES}/tables/records_pooled.tsv.gz",
        reps=f"{RES}/tables/representatives.tsv.gz",
    log:
        f"{RES}/logs/select_representatives.log",
    params:
        py=PY,
        script=SCRIPTS / "select_representatives.py",
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --annotated {input.annotated} '
        "--augment {input.augment} --tier2 {input.tier2} --terms {input.terms} "
        "--out-pooled {output.pooled} --out-reps {output.reps} > {log} 2>&1"


rule write_outputs:
    input:
        cfg=CONFIG_JSON,
        pooled=rules.select_representatives.output.pooled,
        reps=rules.select_representatives.output.reps,
        fasta=rules.fetch_sequences.output.fasta,
        query=rules.build_query.output.query,
        meta=rules.build_query.output.meta,
        terms=IN["search_terms"],
        host_accept=IN["host_accept"],
        tier2=IN["tier2_groups"],
        audit=rules.screen_augmented.output.audit,
    output:
        tier1=f"{RES}/{DATASET}_refdb_tier1.fasta",
        tier2fa=f"{RES}/{DATASET}_refdb_tier2.fasta",
        allfa=f"{RES}/{DATASET}_refdb_all.fasta",
        refs=f"{RES}/{DATASET}_reference_genomes.tsv",
        species=f"{RES}/{DATASET}_species_summary.tsv",
        records=f"{RES}/{DATASET}_all_records.tsv",
        excluded=f"{RES}/{DATASET}_excluded_records.tsv",
        readme=f"{RES}/README_{DATASET}_refdb.md",
    log:
        f"{RES}/logs/write_outputs.log",
    params:
        py=PY,
        script=SCRIPTS / "write_outputs.py",
        outdir=RES,
        dataset=DATASET,
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --config-json {input.cfg} --pooled {input.pooled} '
        "--reps {input.reps} --fasta {input.fasta} --query {input.query} --meta {input.meta} "
        "--terms {input.terms} --host-accept {input.host_accept} --tier2 {input.tier2} "
        "--audit {input.audit} "
        "--outdir {params.outdir} --dataset {params.dataset} > {log} 2>&1"


rule plot_overview:
    input:
        reps=rules.write_outputs.output.refs,
        records=rules.write_outputs.output.records,
    output:
        png=f"{RES}/{DATASET}_refdb_overview.png",
    log:
        f"{RES}/logs/plot_overview.log",
    params:
        py=PY,
        script=SCRIPTS / "plot_overview.py",
        dataset=DATASET,
    conda:
        ENV
    shell:
        '"{params.py}" "{params.script}" --reps {input.reps} --records {input.records} '
        "--out {output.png} --dataset {params.dataset} > {log} 2>&1"
