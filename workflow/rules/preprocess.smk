# Stage A, part 1: one object per sample, QC, integration, clustering.
# Ported from: Jaimar_BCG_Mice_CITEseq_Project/.../hpg_code/{Seurat_preprocess.R,
# combineall_h5_harmony_20240627.R} and Jaimar_CS_mice_Project_20241202/Rcode/
# 1-QC_integration_preprocessing.R (see docs/design.md).

rule import_sample:
    input:
        lambda wc: samples.loc[wc.sample, "path"],
    output:
        f"{OUT}/01_import/{{sample}}.rds",
    params:
        mode=config["input"]["mode"],
    log:
        f"{OUT}/logs/01_import_{{sample}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/01_import_sample.R"

rule qc:
    input:
        f"{OUT}/01_import/{{sample}}.rds",
    output:
        obj=f"{OUT}/02_qc/{{sample}}.rds",
        metrics=f"{OUT}/02_qc/{{sample}}_qc_metrics.tsv",
    params:
        qc=config["qc"],
    log:
        f"{OUT}/logs/02_qc_{{sample}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/02_qc.R"

rule integrate:
    input:
        expand(f"{OUT}/02_qc/{{sample}}.rds", sample=SAMPLES),
    output:
        f"{OUT}/03_integrate/merged_harmony.rds",
    params:
        integration=config["integration"],
        samples=config["samples"],
    threads: 8
    resources:
        mem_mb=config["resources"]["integrate_mem_mb"],
        runtime=config["resources"]["integrate_minutes"],
    log:
        f"{OUT}/logs/03_integrate.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/03_integrate.R"

rule cluster:
    input:
        f"{OUT}/03_integrate/merged_harmony.rds",
    output:
        obj=f"{OUT}/04_cluster/clustered.rds",
        markers=f"{OUT}/04_cluster/markers.tsv",
        report=f"{OUT}/04_cluster/cluster_report.html",
    params:
        clustering=config["clustering"],
    threads: 8
    resources:
        mem_mb=config["resources"]["integrate_mem_mb"],
    log:
        f"{OUT}/logs/04_cluster.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/04_cluster.R"
