# workflow/scripts/

Not written yet. Each script is ported from a script that already ran on the
sepsis data; the source is named in the rule file that calls it and in
`docs/design.md`. A script is added here only together with a run on the
16-mouse fixture.

| script | rule | ported from |
|---|---|---|
| `01_import_sample.R` | import_sample | Seurat `Read10X_h5` / `Read10X` per sample; merged-object and annotated-object entry modes |
| `02_qc.R` | qc | `Seurat_preprocess.R` (2024): mito %, gene-count bounds, optional doublet call |
| `03_integrate.R` | integrate | `combineall_h5_harmony_20240627.R`: normalise, 2,000 HVG, PCA, Harmony on the batch variable |
| `04_cluster.R` | cluster | same: SNN (k = 30) on Harmony dims, Louvain, UMAP, `FindAllMarkers` |
| `05_draft_annotation.R` | draft_annotation | marker scoring per cluster (`5-Refine_celltype_annot.R`), optional reference-based labels |
| `06_apply_annotation.R` | apply_annotation | joins the locked table to the object, writes cell metadata |
| `10_composition.R` | composition | `samplewise::build_composition()` per denominator |
| `11_pseudobulk.R` | pseudobulk | `AggregateExpression(group.by = c(celltype, sample))` as in `04-deg_pseudobulk_edger.Rmd` |
| `12_cellchat_per_sample.R` | cellchat_per_sample | `01-build_cellchat_per_mouse.R` (neonatal rebuild, 2026-09) |
| `13_cellchat_long_tables.R` | cellchat_long_tables | `02-extract_long_tables.Rmd` -> `samplewise::build_cellchat()` |
| `20_fit_derived.R` | fit_derived | `03-inference.Rmd` -> `samplewise::fit_features()` |
| `21_fit_expression.R` | fit_expression | `05.01-DEseq2_reduced_model_20250916.R` -> `samplewise::fit_expression()` |
| `22_gsea.R` | gsea | `06-GSEA_20260129.R`: clusterProfiler on signed statistics per term |
| `30_report.Rmd` | report | `10_build_tables.Rmd` + `04_report.Rmd` (myeloid task) -> `write_results()`, `plot_forest()` |
