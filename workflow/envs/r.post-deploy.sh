#!/usr/bin/env bash
# Runs once after the conda env is created (Snakemake convention).
set -euo pipefail
Rscript -e 'remotes::install_github("jinworks/CellChat", upgrade = "never")'
Rscript -e 'remotes::install_github("dayuan-wang/samplewise", upgrade = "never")'
