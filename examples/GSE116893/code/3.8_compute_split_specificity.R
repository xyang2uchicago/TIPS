########## BEGINNING OF USER INPUT ##########
wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
########## END OF USER INPUT ##########

# Simplifies HiG_9_hi / HiG_9_lo from the notsimplified graph list, computes
# celltype-specificity weights using sce_C9split + leiden_C9split, then patches
# the two graphs into simplified_combinedweighted.rds.
# Must be run AFTER 3.7_compute_split_DEG.R.

library(SingleCellExperiment)
library(igraph)
library(data.table)

celltype_specific_weight_version <- "10"
BioTIP_version                   <- "06232025"

source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

data_dir <- paste0(wd, "data/")
res_path <- paste0(wd, "results/")
ppi_path <- paste0(wd, "results/PPI_weight/")

# ---- load notsimplified, extract split graphs ----
ns_path <- paste0(res_path, "GSE116893_STRING_graph_perState_notsimplified.rds")
g_ns    <- readRDS(ns_path)
stopifnot("HiG_9_hi" %in% names(g_ns), "HiG_9_lo" %in% names(g_ns))

cat("HiG_9_hi (notsimplified): vcount =", vcount(g_ns[["HiG_9_hi"]]),
    "ecount =", ecount(g_ns[["HiG_9_hi"]]), "\n")
cat("HiG_9_lo (notsimplified): vcount =", vcount(g_ns[["HiG_9_lo"]]),
    "ecount =", ecount(g_ns[["HiG_9_lo"]]), "\n")

# ---- step 1: simplify (mirrors 11.2.0 logic) ----
sub_list <- list(
    HiG_9_hi = simplify(g_ns[["HiG_9_hi"]], edge.attr.comb = "max"),
    HiG_9_lo = simplify(g_ns[["HiG_9_lo"]], edge.attr.comb = "max")
)

for (nm in names(sub_list)) {
    g <- sub_list[[nm]]
    cat(nm, "after simplify: ecount =", ecount(g),
        "| multi-edges remaining:", any(which_multiple(g)),
        "| duplicates:", any(duplicated(V(g)$name)), "\n")
}

# ---- step 2: compute specificity using leiden_C9split ----
sce <- readRDS(paste0(data_dir, "sce_C9split.rds"))
if (!"logcounts" %in% assayNames(sce)) {
    library(scater)
    sce <- scater::logNormCounts(sce)
}
colData(sce)$cluster <- colData(sce)[["leiden_C9split"]]

cat("\nComputing network specificity for HiG_9_hi and HiG_9_lo...\n")
network_specificity_split <- calculate_network_specificity(
    sce,
    sub_list,
    assayName    = "logcounts",
    celltype_col = "cluster",
    method       = "pearson",
    cores        = 1,
    shrink       = TRUE,
    verbose      = TRUE
)

spec_path <- paste0(ppi_path, "network_specificity_9split.rds")
saveRDS(network_specificity_split, spec_path)
cat("Saved split specificity to:", spec_path, "\n")

# ---- step 3: update edge weights ----
weighted_sub <- update_network_weights(
    sub_list,
    network_specificity_split,
    specificity_method = "combined",
    verbose            = FALSE,
    cores              = 1
)

for (nm in names(weighted_sub)) {
    g <- weighted_sub[[nm]]
    cat(nm, "after weighting: edge attrs =", paste(edge_attr_names(g), collapse=", "),
        "| weight range", paste(round(range(E(g)$weight), 4), collapse=" - "), "\n")
}

# ---- step 4: patch into simplified_combinedweighted ----
grn_path <- paste0(ppi_path, "GSE116893_STRING_graph_perState_simplified_combinedweighted.rds")
g_sw <- readRDS(grn_path)

g_sw[["HiG_9_hi"]] <- weighted_sub[["HiG_9_hi"]]
g_sw[["HiG_9_lo"]] <- weighted_sub[["HiG_9_lo"]]

saveRDS(g_sw, file = grn_path)
cat("\nSaved updated graph_list to:", grn_path, "\n")
cat("Entries:", paste(names(g_sw), collapse=", "), "\n")
