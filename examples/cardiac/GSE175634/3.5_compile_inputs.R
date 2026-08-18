library(BioTIP)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
# This code recreates the BioTIP.res object used in downstream analysis. Method is based off of 11.1_STRING

source(here::here("examples", "config.R"))
wd         <- tips_path("examples", "cardiac", "GSE175634/")
data_dir   <- paste0(wd, "data/")
output_dir <- paste0(wd, "data/")
########## END OF USER INPUT ##########

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sce <- readRDS(paste0(data_dir, "sce.rds"))

# Copy sce.rds to output_dir so that 24.0 can find it via data_dir
if (!file.exists(paste0(output_dir, "sce.rds"))) {
    file.copy(paste0(data_dir, "sce.rds"), paste0(output_dir, "sce.rds"))
    cat("Copied sce.rds to", output_dir, "\n")
}

# Load CTS candidates and MCI scores
load(paste0(data_dir, "CTS_Lib_Scaledata.RData"))
# Loads: CTS.Lib, CTS.Lib.Symbol, maxMCI
# names(CTS.Lib.Symbol): muscle, endoderm, CP/CF, CM, CF, CP, mix pro, endoderm.1, CP.1

# Load BioTIP IC scores (BioTIP_scores) — used only for res$Ic.shrink; no simulation rerun
load(paste0(data_dir, "LibSF_IC_sim_BioTIP_Scaledata.RData"))
# Loads: BioTIP_scores, SimResults_g

# Significance determined from prior BioTIP analysis recorded in 11.1_STRINGweighted_CTS_cardiac_network.R
# Order matches names(CTS.Lib.Symbol): muscle, endoderm, CP/CF, CM, CF, CP, mix pro, endoderm.1, CP.1
sig.DNB     <- c(T,  T,  T,    T,  T,  T, F,        F,          T)
sig.IC      <- c(T,  T,  F,    F,  T,  T, NA,       NA,         T)
sig.deltaIC <- c(T,  T,  NA,   NA, F,  T, NA,       NA,         T)

significant <- sig.DNB == TRUE & !is.na(sig.IC) & sig.IC == TRUE & !is.na(sig.deltaIC) & sig.deltaIC == TRUE
names(significant) <- names(CTS.Lib.Symbol)

cat("Significant CTS:\n")
print(significant)
# Significant CTS:
#     muscle   endoderm      CP/CF         CM         CF         CP    mix pro 
#       TRUE       TRUE      FALSE      FALSE      FALSE       TRUE      FALSE 
# endoderm.1       CP.1 
#      FALSE       TRUE 

# Assemble res object
res <- list(
    CTS.candidate = CTS.Lib.Symbol,
    CTS.score     = maxMCI,
    Ic.shrink     = BioTIP_scores,
    significant   = significant
)

save(res, file = paste0(output_dir, "BioTIP.res.RData"))
cat("BioTIP.res.RData saved to", paste0(output_dir, "BioTIP.res.RData"), "\n")
