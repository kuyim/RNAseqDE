# =============================================================================
# Dockerfile — reproducible environment for the RNA-seq DE + enrichment analysis
#
# Start from the official Bioconductor image (ships R + BiocManager +
# the system libraries the Bioc stack needs), then let renv install the EXACT
# package versions recorded in renv.lock. This reproduces the analysis
# environment on any machine.
#
# =============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_23

WORKDIR /project

# --- Restore the exact package set from the lockfile --------------------------
# Copy ONLY renv.lock first, so this expensive layer is cached and only rebuilds
# when your dependencies actually change (not when you edit a script).
COPY renv.lock renv.lock
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')" \
 && R -e "renv::restore()"

# --- Copy the analysis code ---------------------------------------------------
# The R scripts live at the project root. Data and results are mounted at
# runtime (see the docker run command), so they are NOT baked into the image —
# keeps it small and the data re-fetchable.
COPY diff_exp.R enrichment.R ./

# --- Default action: run the full analysis end to end -------------------------
# Reads the Salmon count matrix from the mounted results/ dir, writes DE +
# enrichment outputs back into it.
CMD ["sh", "-c", "Rscript diff_exp.R && Rscript enrichment.R"]
