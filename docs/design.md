# samplewise — design

Written 2026-09-27 by Dayuan Wang with Claude. This is the document the code is
built against. When code and this document disagree, fix one of them the same
day.

## 1. Purpose

`samplewise` is the general form of the statistical pipeline used in Wang
(2026), *Copy number variant detection with deep learning and single-cell
analysis of immune dysregulation in sepsis*, Chapter 3: a single-cell RNA-seq
analysis in which every downstream feature, not only gene expression, is tested
with one design formula and with the biological sample as the unit of
inference. The sepsis study asked one question of five feature types (cell-type
composition, pseudobulk expression, pathway enrichment, CellChat signalling,
RNA velocity, Drug2Cell scores): does the effect of sepsis differ by age and by
sex? The answer had to be computed the same way each time, and it was, but by
three generations of study-specific scripts. This repository turns that into
software that runs on other designs, other tissues and the other species.

It is two things in one repository:

- an **R package** (the repository root) that owns the statistics: the
  per-sample table every feature is reduced to, the builders that produce it
  from standard tools, one fitting function with explicit eligibility rules,
  and one report;
- a **Snakemake workflow** (`workflow/`, `config/`) that owns the plumbing:
  from Cell Ranger output to annotated object, then every module and the
  report, on a laptop or on SLURM.

What it is not: a new preprocessing method, a new communication or velocity
method, or an autonomous agent. It orchestrates validated tools and makes the
inference on top of them correct by default. The 2026-09-24 feasibility review
of the wider "sepsis platform" idea (HiPerGator,
`sepsis-framework-project-202609/`) reached the same conclusion: the defensible
contribution is the covariate-aware statistical layer, not a model.

## 2. What generalises, and what the sepsis code had hard-coded

| hard-coded in the study | in samplewise |
|---|---|
| `~ age + sex + sepsis + age:sepsis + sex:sepsis`, reference young male naive | `design.formula`, `design.reference_levels`, `design.contrasts` in `config.yaml`; every module receives the same three |
| mouse id, age, sex, group scattered across scripts | one `samples.tsv`: `sample_id`, `path`, then any design columns |
| mm10, `CellChatDB.mouse`, MSigDB mouse, ChEMBL human-to-mouse homologs | `species: mouse | human`, switched in one place |
| six lymphocyte types, APC-to-T pathways | annotation column, compartments (denominators), sender and receiver sets are configuration |
| DESeq2 for expression; lm/glm for proportions and drug scores; limma on the sqrt scale for CellChat; LMM for velocity | one `fit_features()` for every per-sample feature with a per-family transform; `fit_expression()` with a switchable engine; velocity (cell-level, mixed model) deferred to v2 |
| eligibility rules stated in three different places | `inference.min_per_group`, `inference.min_detected`, per-module `min_cells`, `min_denominator` |

Three things do not generalise by configuration and are handled by design:

1. **Annotation is a human step.** The sepsis annotation was manual, marker
   based, and re-done for the neonatal cohort. The workflow therefore has two
   stages with a checkpoint between them: stage A ends with a *draft*
   annotation (marker scores per cluster, optionally reference-based labels);
   a person writes the locked annotation file; stage B refuses to start
   without it. `accept_draft: true` exists for fixtures whose annotation is
   known, and for nothing else.
2. **The inference layer had three implementations** (lymphocyte 2025,
   myeloid 2026-08, neonatal 2026-09). They agree on the principle and differ
   in transform, engine and thresholds. Section 5 states the unified rules;
   where the study left a choice undocumented (the link for proportions and
   drug scores, defense change log E4) the default is what the published paper
   did, and the simulation study planned for v2 is what will settle it.
3. **Two of the five modules are notebooks** (scVelo, Drug2Cell with
   collaborator code). They are v2.

## 3. Architecture

```
samples.tsv + config.yaml
        │
        ▼
 stage A  import ─ qc ─ integrate(Harmony) ─ cluster ─ draft annotation
        │                                              (report + tsv)
        ║  ─── a person locks the annotation ───
        ▼
 stage B  apply annotation
        ├─ composition  ──┐
        ├─ pseudobulk ────┼──▶ per-sample tables ──▶ fit_features() / fit_expression()
        ├─ CellChat/sample┘        (the contract)              │
        │                                                      ▼
        └──────────────────────────────────────────────▶ report (xlsx + html)
```

Entry points (`input.mode`): per-sample Cell Ranger `.h5`, per-sample 10x MTX
directories, one merged Seurat object with a sample column, or an already
annotated object (stage A skipped). The last two matter for adoption: most
groups arrive with an object, not with FASTQ.

The package is usable without the workflow. `build_*()` functions accept the
objects people have (a Seurat object, a data.frame of cell metadata, a list of
per-sample CellChat objects), and `fit_features()` accepts any table in the
contract, including ones users assemble themselves.

## 4. Contracts

**Per-sample long table** (`sw_long_table()`, `sw_validate()`): one row per
sample × feature. Columns `sample_id`, `feature_type`, `unit`, `feature`,
`value`, `n_cells`, `estimable`. A value that is not estimable in a sample is
`NA` and is excluded from the fit; a zero is a measured zero. Examples:

| feature_type | unit | feature | value |
|---|---|---|---|
| composition | `B` | `proportion_of_all` | 0.41 |
| composition | `CD4 T` | `proportion_of_lymphocytes` | 0.22 |
| communication | `B>CD4 T` | `MHC-II` | 0.031 (pathway strength, sqrt-transformed at fit time) |
| drug | `CD8 T` | `bexarotene` | mean Drug2Cell score |

**Pseudobulk contract** (`fit_expression()`): a genes × pseudobulk count
matrix plus a column table with `sample_id`, `unit` (cell type) and `n_cells`.
Expression does not go through the long table (genes × cell types × samples is
too wide) but its result table has the same columns as `fit_features()`, with
`feature` = gene.

**Sample sheet**: `sample_id`, `path`, then design columns. Character columns
are factors whose reference level is the first entry of
`design.reference_levels`, or the first level alphabetically when unspecified.

**Contrasts**: `name <TAB> coefficients`, where `coefficients` is a
comma-separated list of `coefficient=weight` over the names
`sw_coef_names(samples, formula)` returns. Empty means every non-intercept
coefficient.

**Results**: `feature_type, unit, feature, term, estimate, se, statistic, df,
p, fdr, n_fit, n_detected, status`. `status` is `ok`, `ok_with_warning`, or
the name of the eligibility rule that stopped the fit. Nothing is silently
dropped: a series that was not fitted keeps a row.

## 5. Inference rules (unified)

Rules carried over from the three implementations; the source of each is
named so it can be checked.

- **Unit of inference is the sample.** Every per-sample feature is one value
  per sample; the model is a GLM on those values (`stats::glm`). A random
  intercept per sample is not identifiable for one value per sample and is
  not used; it belongs to cell-level features (velocity, v2) and to designs
  with several samples per subject (v2, `random =`).
- **Estimable vs detected** (neonatal CellChat rebuild, 2026-09-12): a
  CellChat value is estimable in a sample when sender and receiver both have
  more than `min_cells` cells there; not estimable is `NA`, never 0; detected
  means value > 0 among estimable samples. Composition: estimable when the
  denominator has at least `min_denominator` cells.
- **Eligibility to test**: at least `min_per_group = 3` estimable samples in
  every level of every factor in the formula (levels defined by the sample
  sheet, so a missing level fails the rule rather than being dropped), and
  detected in at least `min_detected = 5` samples (neonatal rule; without it,
  circuits non-zero in 2–3 samples of one group reach FDR < 0.05 under
  moderation). Pseudobulks from fewer than `min_cells = 10` cells are dropped
  (lymphocyte and neonatal DE).
- **Transforms by family**: composition identity on the proportion (as
  published; `propeller` on the same design is the sensitivity check),
  communication `sqrt` (neonatal rebuild), drug scores identity. All
  overridable per family in `inference.transform`.
- **Expression engine**: DESeq2 by default (as published), edgeR
  quasi-likelihood (neonatal step 4) and dreamlet (precision-weighted mixed
  model, for repeated samples) as options behind one interface.
- **Multiple testing**: Benjamini-Hochberg within one family, default
  `feature_type × term` (all cell types of a denominator; all pathway × edge
  combinations of a contrast); expression within cell type × term. Raw p is
  always reported.
- **Sensitivity checks** (on by default, reported beside the primary result,
  never replacing it): `propeller` for composition; log cells per sample as a
  covariate for communication (neonatal rebuild).
- **Reporting**: effect size, standard error, CI, raw p, FDR on one scale for
  every family; xlsx with a README sheet first and one sheet per term, rows
  sorted by raw p (the format collaborators in this project receive).

Requirements the sepsis design did not exercise but the lab's next data will
(feasibility addendum, 2026-09-24): a blocking factor above the sample
(litter, cage, donor) and batch as a covariate. Batch is a formula term today;
random effects above the sample are a v2 interface (`fit_features(random =
~ 1 | litter)` via lme4 or glmmTMB) and a reason dreamlet is one of the
expression engines.

## 6. Modules

| module | v1 | ported from (HiPerGator) | package function |
|---|---|---|---|
| import, QC, Harmony, clustering | yes | `Jaimar_BCG_Mice_CITEseq_Project/dropbox-upload-20260819-NS3659_data_process/All_Adult_data_analysis/hpg_code/{Seurat_preprocess.R, combineall_h5_harmony_20240627.R}`; `Jaimar_CS_mice_Project_20241202/Rcode/1-QC_integration_preprocessing.R` | workflow only |
| draft annotation | yes | `.../local R code/5-Refine_celltype_annot.R`; neonatal `01-annotation/code/annotation_update_202609.Rmd` | workflow only |
| composition | yes | `myeloid-aging-sepsis-mice-202609/02-cell-proportion/rcode/{01_extract_counts.R, 02_fit_models.R}` | `build_composition()` (done), `fit_features()` (done) |
| pseudobulk expression | yes | `Dayuan_lymphocyte_paper/05-pseudobulk_DE/code/05.01-DEseq2_reduced_model_20250916.R`; `neonatal-bcg-cs-mice-202608/04-deg/code/04-deg_pseudobulk_edger.Rmd` | `fit_expression()` (interface) |
| GSEA | yes | `Dayuan_lymphocyte_paper/06-GSEA_pathway/code/06-GSEA_20260129.R` | workflow (clusterProfiler) |
| CellChat per sample + long tables | yes | `neonatal-bcg-cs-mice-202608/02-cellchat/code/{01-build_cellchat_per_mouse.R, 02-extract_long_tables.Rmd}`; `03-inference/code/03-inference.Rmd` | `build_cellchat()` (interface), `fit_features()` |
| report | yes | `myeloid-aging-sepsis-mice-202609/05-tables/rcode/10_build_tables.Rmd`, `04-delivery/rcode/04_report.Rmd` | `write_results()`, `plot_forest()` (done) |
| RNA velocity | v2 | `Dayuan_lymphocyte_paper/04-scvelo/` (3 notebooks + LMM scripts) | cell-level mixed model |
| Drug2Cell | v2 | `Dayuan_lymphocyte_paper/08-drug2cell/` (Kenn's notebooks + R) | `build_drug()` |
| sensitivity grid across method choices | v2 | new | |

## 7. Existing tools, and why this exists anyway

Checked 2026-09-27.

- **dreamlet** (Hoffman et al., *Nature Communications* 17:9597, 2026;
  Bioconductor) and **crumblr** (Hoffman & Roussos, bioRxiv 2025.01.29.635498):
  pseudobulk expression with precision-weighted linear mixed models, and
  compositional analysis with the same machinery. Expression and composition
  only. dreamlet is one of samplewise's expression engines.
- **propeller / speckle** and **scCODA**: composition only. propeller is the
  built-in sensitivity check.
- **MultiNicheNet** (Browaeys et al., bioRxiv 2023.06.13.544751;
  multinichenetr 2.0, 2024): multi-sample, multi-condition cell-cell
  communication with sample-level modelling. The closest tool to one module;
  it does not cover the other families and is a NicheNet-based method rather
  than a wrapper around the user's chosen tool.
- **pertpy** (Heumos et al., *Nature Methods*, 2025): Python; compositional
  and pseudobulk modules among many perturbation tools.
- **nf-core/scrnaseq**: FASTQ to count matrices. Referenced, not reimplemented;
  samplewise starts at the matrices.
- **muscat**: pseudobulk DE and simulation of multi-sample designs; its
  simulator is a candidate for the v2 simulation study.

None of them puts composition, expression, pathways, communication (and, in
v2, velocity and drug scores) under one declared design with one eligibility
policy and one report. That is the claim, and it is a modest one: the value is
in doing the ordinary thing correctly every time.

## 8. Validation

1. **Unit tests** (`tests/testthat`): contract, composition arithmetic,
   recovery of a known effect, every eligibility rule, contrasts by name and
   by vector, report writers. Run by GitHub Actions on every push.
2. **Fixture** (`00-fixture/` on HiPerGator; distributed as a release asset):
   16 of the 57 sepsis mice, 2 per cell of the age × sex × sepsis design,
   1,000 cells each drawn from the cells the paper kept, written as 10x MTX
   directories with the paper's annotation for those cells. Exercises every
   rule end to end in minutes. Source: the per-sample `filtered_feature_bc_matrix.h5`
   under `Valerie_Adult_mice_all_batch_output/ALL_4_Batches_h5_files/` and
   `Dayuan_lymphocyte_paper/R_objects/seurat_object_exclude_Platelet_RBC_20250903.rds`.
3. **Full-data regression** (HiPerGator only): the 57-mouse object through
   stage B with `annotated_object` input; the composition, expression and
   CellChat results must match the paper's Figures 2, 3 and 5 tables to
   numerical tolerance. This is the test that the port changed nothing.
4. **Generalisation**: Kang et al. 2018 (GEO GSE96583; `muscData::Kang18_8vs8`),
   PBMC from 8 lupus patients, control vs 6 h IFN-β, human. Cells are
   demultiplexed to patients, so the merged-object entry mode applies, with
   `sample_id = patient × stimulation` (16 samples) and `~ stim + patient`
   (patient as a blocking term). Different species, different design, and the
   dataset muscat and dreamlet used, so published pseudobulk results exist to
   compare against. Exit criterion: config changes only, no code changes.

## 9. Roadmap

| stage | content | exit criterion |
|---|---|---|
| 0 (done 2026-09-27) | this document, package skeleton with the contract, `build_composition()`, `fit_features()`, report; workflow rule graph; config; fixture job | a reader of the README knows the inputs and outputs; tests pass |
| 1 | port the v1 modules and the workflow scripts; fixture run end to end; full-data regression | paper numbers reproduced |
| 2 | Kang18 end to end | config-only run on human data |
| 3 | velocity and Drug2Cell modules; random effects above the sample; sensitivity grid; simulation study (type I error, power, imbalance, batch, misspecified link) | the two items Chapter 4 lists as future work |

## 10. Decisions (2026-09-27)

- Name `samplewise`: free on CRAN, Bioconductor and GitHub; a PyPI project of
  that name exists (unrelated); a Python companion would be `samplewise-py`.
- License GPL-3 (CellChat, an intended dependency, is GPL-3). Can be revisited
  before the first public release.
- Home `github.com/dayuan-wang/samplewise`, private until stage 1 runs; local
  clone inside the scirc-work task folder, HiPerGator clone at
  `/orange/paefron/dayuan.wang/scrnaseq-inference-pipeline-202610/samplewise`.
- v1 excludes velocity and Drug2Cell. Expression engine default DESeq2.
- Second dataset: Kang18. Fixture: 16 mice × 1,000 cells.
- Workflow engine Snakemake (mixed R/Python, SLURM, conda per rule);
  FASTQ-to-matrix left to nf-core/scrnaseq, with an optional Cell Ranger rule
  later for groups that hold FASTQ.
