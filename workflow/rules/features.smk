# Stage B, part 1: per-sample feature matrices (samplewise long tables).
# Ported from: myeloid-aging-sepsis-mice-202609/02-cell-proportion (composition),
# Dayuan_lymphocyte_paper/05-pseudobulk_DE + neonatal-bcg-cs-mice-202608/04-deg
# (pseudobulk), neonatal-bcg-cs-mice-202608/02-cellchat (CellChat per sample and
# its long tables).

rule composition:
    input:
        f"{OUT}/05_annotation/cell_metadata.tsv.gz",
    output:
        f"{OUT}/06_features/composition.tsv",
    params:
        composition=config["modules"]["composition"],
    log:
        f"{OUT}/logs/10_composition.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/10_composition.R"

rule pseudobulk:
    input:
        f"{OUT}/05_annotation/annotated.rds",
    output:
        counts=f"{OUT}/06_features/pseudobulk_counts.rds",
        coldata=f"{OUT}/06_features/pseudobulk_coldata.tsv",
    params:
        expression=config["modules"]["expression"],
    resources:
        mem_mb=config["resources"]["integrate_mem_mb"],
    log:
        f"{OUT}/logs/11_pseudobulk.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/11_pseudobulk.R"

rule cellchat_per_sample:
    input:
        f"{OUT}/05_annotation/annotated.rds",
    output:
        f"{OUT}/06_features/cellchat/{{sample}}.rds",
    params:
        communication=config["modules"]["communication"],
        species=config["species"],
    resources:
        mem_mb=config["resources"]["cellchat_mem_mb"],
        runtime=config["resources"]["cellchat_minutes"],
    log:
        f"{OUT}/logs/12_cellchat_{{sample}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/12_cellchat_per_sample.R"

rule cellchat_long_tables:
    input:
        expand(f"{OUT}/06_features/cellchat/{{sample}}.rds", sample=SAMPLES),
    output:
        expand(f"{OUT}/06_features/communication_{{level}}.tsv",
               level=config["modules"]["communication"]["levels"]),
    params:
        communication=config["modules"]["communication"],
    log:
        f"{OUT}/logs/13_cellchat_long_tables.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/13_cellchat_long_tables.R"
