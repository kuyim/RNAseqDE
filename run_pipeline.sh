#!/bin/bash
set -euo pipefail
export NXF_VER=25.04.6
CONFIG=conf/run.config

# Download data
nextflow run nf-core/fetchngs -r 1.12.0 \
  -c "$CONFIG" -profile docker,download -resume

# Run RNAseq
nextflow run nf-core/rnaseq -r 3.26.0 \
  -c "$CONFIG" -profile docker,quant -resume