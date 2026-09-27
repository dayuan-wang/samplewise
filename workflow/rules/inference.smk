# Stage B, part 2: one design for every feature, then one report.
# Ported from: neonatal-bcg-cs-mice-202608/03-inference (derived features),
# Dayuan_lymphocyte_paper/05-pseudobulk_DE + 06-GSEA_pathway (expression),
# myeloid-aging-sepsis-mice-202609/05-tables (workbook layout).

rule fit_derived:
    input:
        composition=f"{OUT}/06_features/composition.tsv",
        communication=expand(f"{OUT}/06_features/communication_{{level}}.tsv",
                             level=config["modules"]["communication"]["levels"]),
    output:
        f"{OUT}/07_inference/derived_results.tsv",
    params:
        design=config["design"],
        inference=config["inference"],
        samples=config["samples"],
    log:
        f"{OUT}/logs/20_fit_derived.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/20_fit_derived.R"

rule fit_expression:
    input:
        counts=f"{OUT}/06_features/pseudobulk_counts.rds",
        coldata=f"{OUT}/06_features/pseudobulk_coldata.tsv",
    output:
        f"{OUT}/07_inference/expression_results.tsv",
    params:
        design=config["design"],
        expression=config["modules"]["expression"],
        samples=config["samples"],
    resources:
        mem_mb=config["resources"]["integrate_mem_mb"],
    log:
        f"{OUT}/logs/21_fit_expression.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/21_fit_expression.R"

rule gsea:
    input:
        f"{OUT}/07_inference/expression_results.tsv",
    output:
        f"{OUT}/07_inference/gsea_results.tsv",
    params:
        gsea=config["modules"]["gsea"],
        species=config["species"],
    log:
        f"{OUT}/logs/22_gsea.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/22_gsea.R"

rule report:
    input:
        derived=f"{OUT}/07_inference/derived_results.tsv",
        expression=f"{OUT}/07_inference/expression_results.tsv",
        gsea=f"{OUT}/07_inference/gsea_results.tsv",
    output:
        xlsx=f"{OUT}/report/results.xlsx",
        html=f"{OUT}/report/index.html",
    params:
        design=config["design"],
    log:
        f"{OUT}/logs/30_report.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/30_report.Rmd"
