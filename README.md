# samplewise

Sample-level statistical inference for derived single-cell features. One design
formula across cell-type composition, pseudobulk differential expression,
pathway enrichment and cell-cell communication, with the biological sample
(mouse, donor, patient) as the unit of inference, and one report at the end.

**Status: pre-alpha, 2026-09-27.** The statistical core works and is tested
(`sw_long_table()`, `build_composition()`, `fit_features()`, `write_results()`,
`plot_forest()`); the workflow is a wired rule graph whose scripts are being
ported; `fit_expression()` and `build_cellchat()` are documented interfaces
without implementations. Nothing here has run on a real dataset yet. Read
[`docs/design.md`](docs/design.md) before relying on anything.

## Where it comes from

The replicate-aware framework of Wang (2026), *Copy number variant detection
with deep learning and single-cell analysis of immune dysregulation in sepsis*,
Chapter 3, and the manuscript *Age but not sex modifies lymphoid immune
responses in murine sepsis* (57 mice, 828,744 splenocytes, GEO GSE311155).
That study tested one question, does the effect of sepsis differ by age and by
sex, across five feature types with one model. This package makes the design,
the species, the cell types and the eligibility rules configuration instead of
code.

## The idea in one call

```r
library(samplewise)

# cells: a data.frame (or Seurat object) with one row per cell
comp <- build_composition(cells, sample_col = "sample_id", celltype_col = "cell_type")

# samples: one row per sample, with the design columns
res <- fit_features(comp, samples, ~ age * sepsis + sex,
                    min_per_group = 3, min_detected = 5)

write_results(res, "composition.xlsx")
plot_forest(res, term = "ageold:sepsissepsis")
```

Every feature type passes through the same table shape (`?sw_long_table`), is
fitted by the same function with the same eligibility rules (`?fit_features`),
and lands in the same result table. Series that could not be tested keep a row
that says why.

## Layout

```
R/            the package: contract, builders, fitting, report
tests/        testthat
workflow/     Snakemake rules, environments, HiPerGator profile, scripts (being ported)
config/       config.yaml, samples.tsv, contrasts.tsv, marker list
docs/         design.md
inst/extdata/ example sample sheet (16-mouse fixture) and contrasts
```

## Installation

```r
remotes::install_github("dayuan-wang/samplewise")
```

## Citation and license

See `CITATION.cff`. GPL-3.
