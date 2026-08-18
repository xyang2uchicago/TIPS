# examples/config.R
#
# Resolves the TIPS repo root once, so analysis scripts under examples/
# never need to hardcode a machine-specific absolute path (e.g.
# "/Users/felixyu/Documents/GitHub/TIPS/...") in their USER INPUT block.
#
# How it works: `here::here()` walks up from the current working
# directory looking for a repo marker (.git/, an .Rproj file, etc.) and
# returns the first match. That means this works unchanged on any computer,
# as long as R is started from somewhere inside the repo checkout (the normal
# case for an RStudio project or `Rscript` run from the repo).
#
# If a script ever runs from outside the repo tree (e.g. a batch job
# whose working directory is elsewhere), set TIPS_ROOT before sourcing
# this file instead of relying on auto-detection:
#   Sys.setenv(TIPS_ROOT = "/project/xyang2/TIPS")
#
# Usage in an analysis script, replacing the old hardcoded wd:
#   source(here::here("examples", "config.R"))
#   wd          <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
#   shared_path <- paste0(shared_data_path, "/")

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required (install.packages(\"here\")) to resolve repo-relative paths.")
}

.tips_root <- Sys.getenv("TIPS_ROOT", unset = NA)
if (is.na(.tips_root)) {
  .tips_root <- here::here()
}

# Build an absolute path from repo-root-relative components, e.g.
#   tips_path("examples", "Shared_Data", "coding_genes.rds")
tips_path <- function(...) file.path(.tips_root, ...)

shared_data_path <- tips_path("examples", "Shared_Data")
