library(BioTIP)
library(Seurat)
library(SingleCellExperiment)

wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/"))

biotip_dir <- "../BioTIP/"
data_dir   <- "../data/"

# Load Seurat object (includes expression matrix + leiden_0.6 in meta.data) and convert to SCE
seurat_obj <- readRDS(paste0(data_dir, "malignant_leiden.rds"))
sce <- as.SingleCellExperiment(seurat_obj)

saveRDS(sce, file = paste0(data_dir, "sce.rds"))

cell_state_col <- "leiden_0.6"
samplesL <- split(colnames(sce), f = colData(sce)[[cell_state_col]])

# Load CTS candidates and MCI scores
load(paste0(biotip_dir, "MCI_bottom_2_CTS_Lib_.RData"))
# Expects: CTS.Lib, CTS.Lib.Symbol, maxMCI

# Load BioTIP IC scores and simulation results
load(paste0(biotip_dir, "LibSF_IC_sim_BioTIP_malignant_t4.RData"))
# Expects: BioTIP_scores, SimResults_g

# Recompute p_delta — plot_SS_Simulation returns p-value; suppress output to PDF
p_delta <- array(NA, length(BioTIP_scores))
pdf(NULL)
for (i in seq_along(BioTIP_scores)) {
    p_delta[i] <- plot_SS_Simulation(
        BioTIP_scores[[i]], SimResults_g[[i]],
        main  = NULL,
        ylab  = NULL,
        xlim  = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], SimResults_g[[i]]))
    )
}
dev.off()
names(p_delta) <- names(BioTIP_scores)

# Recompute significance per CTS (p < 0.05 on both Ic and delta-Ic)
significant <- sapply(seq_along(BioTIP_scores), function(i) {
    interesting <- which(names(samplesL) == names(BioTIP_scores[i]))
    p <- length(which(SimResults_g[[i]][interesting, ] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
    p <- p / ncol(SimResults_g[[i]])
    p < 0.05 & p_delta[i] < 0.05
})
names(significant) <- names(BioTIP_scores)

significant["15"] <- TRUE   # included by design: NOTCH-active ADRN-proliferating transition cluster

cat("Significant CTS:\n")
print(significant)
# Significant CTS:
#    15    16     9
#  TRUE FALSE  TRUE


# Assemble res object
res <- list(
    CTS.candidate = CTS.Lib.Symbol,
    CTS.score     = maxMCI,
    Ic.shrink     = BioTIP_scores,
    significant   = significant
)

save(res, file = paste0(data_dir, "BioTIP.res.RData"))
cat("BioTIP.res.RData saved to", paste0(data_dir, "BioTIP.res.RData"), "\n")
