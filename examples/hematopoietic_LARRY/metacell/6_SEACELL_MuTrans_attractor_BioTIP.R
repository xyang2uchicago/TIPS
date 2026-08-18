
## 6_SEACELL_MuTrans_attractor_BioTIP.R
## Larry metacell MuTrans h5ad -> BioTIP (attractor states)
## https://htmlpreview.github.io/?https://github.com/xyang2uchicago/BioTIP/blob/master/tutorial/Gastrulation.html
##
## July 25 results folder (all steps write here):
##   /project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/larry_BioTIP/BioTIP_attractor/
##
## Interactive use: run the CONFIG block below first, then step flags + rest of script.
## Companion: 6.0_data_MuTrans_annadata_tomatrix.R (same directory as this script on Midway code/)

##############################
## BioTIP on Midway3 — Larry SEACell / MuTrans / attractor
##############################
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\6_SEACELL_MuTrans_attractor_BioTIP.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/.
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\6.0_data_MuTrans_annadata_tomatrix.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/.
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\code\README_larry_BioTIP.md   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/.

# ssh xyang2@midway3.rcc.uchicago.edu
# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy/

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l
# squeue -p caslake --state=PD | wc -l

# sinteractive -p bigmem  --mem=500G  --account=pi-xyang2 --time=1:00:00  -c 1    # requires n_cores  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 -c 3 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 3  
# cd /project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/larry_BioTIP/
# R
## To quit your interactive job:
# exit or Ctrl-D

step0 = F  # h5ad -> Seurat QC + save seu_attractor_MuTrans_HVG.rds
step1 = F   # -> optimized_test_sd_selectionxxx.RData
step2 = F   # -> igraphL_xxx.RData
step3 = F   # -> membersL_xxx.RData, MCI bar plot
step4 = F  # -> xxxSimuMCI_xxx.pdf   done by s1  bigmem c3  500G
step5 = F   # -> SimResults_g_xx.rds  done by s2 amd c1 180G
step6 = F  # -> SimResults_s_xx.rds  bigmem or amd-hm

step7 = TRUE   # summary and plot -> LibSF_IC_sim_BioTIP_",input,"_g.pdf / _s.pdf / _2.pdf (gene + cell perm)
step8 = TRUE   # 1x4 UMAP: Cell.type.clean, attractor, entropy, land

n_cores <- 1
n_cores_permu <- 1
doParallel_permu <- FALSE

## --- CONFIG: run this block first when executing line-by-line interactively ---
LARRY_BIOTIP_ROOT <- Sys.getenv(
	"LARRY_BIOTIP_ROOT",
	"/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking"
)
results_dir <- file.path(LARRY_BIOTIP_ROOT, "results", "larry_BioTIP")
setwd(results_dir)
cat("getwd() =", getwd(), "\n")

LARRY_BIOTIP_CODE_DIR <- Sys.getenv(
	"LARRY_BIOTIP_CODE_DIR",
	file.path(LARRY_BIOTIP_ROOT, "code")
)

#setwd('F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/')
library('SingleCellExperiment')
library(Seurat)
library(dplyr) 
# devtools::install_github("xyang2uchicago/BioTIP")
library(BioTIP)
packageVersion('BioTIP')  # 1.16.0
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

library(foreach)
library(doParallel)

source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/BioTIP_update_06232025.R')  # enable n_cores  = 3
#   source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/BioTIP_update_12092024.R')
#   source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/optimize.sd_selection3.R')
#   source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/simulation_Ic2.R')

# source('E:/Git_Holly/BioTIP/R2/BioTIP_update_12092024.R')

h5ad_helper_candidates <- c(
	file.path(LARRY_BIOTIP_CODE_DIR, "6.0_data_MuTrans_annadata_tomatrix.R"),
	"F:/projects/TIPS/source/GSE140802_lineage_tracking/6.0_data_MuTrans_annadata_tomatrix.R"
)
h5ad_helper <- h5ad_helper_candidates[file.exists(h5ad_helper_candidates)][1]
if (is.na(h5ad_helper) || !nzchar(h5ad_helper)) {
	stop("Cannot find 6.0_data_MuTrans_annadata_tomatrix.R; set LARRY_BIOTIP_CODE_DIR")
}
source(h5ad_helper)
cat("sourced h5ad helpers from:", h5ad_helper, "\n")

## h5ad -> Seurat (BioTIP expects Seurat RNA layers counts + data)
## MuTrans_HVG_input = TRUE  -> 3000 HVG (larry_seacells_mutrans.h5ad)
## MuTrans_HVG_input = FALSE -> all genes; run scripts/larry_xy.py on Midway first
MuTrans_HVG_input <- TRUE

h5ad_mutrans  <- '/project/xyang2/mohsen/experiments/data/processed/larry/mutrans/larry_seacells_mutrans.h5ad'
h5ad_allgenes <- '/project/xyang2/mohsen/experiments/data/processed/larry/mutrans/larry_seacells_log1p_allgenes.h5ad'

input <- 'data'   # log1p layer "data" (not scaled; matches MuTrans preprocess)
r <- 'attractor'  # h5ad obs column for BioTIP state partition

if (MuTrans_HVG_input) {
	h5ad_obj <- h5ad_mutrans
	export_subdir <- "BioTIP_attractor"
	cat("BioTIP input: MuTrans h5ad (3000 HVG)\n")
} else {
	h5ad_obj <- h5ad_allgenes
	export_subdir <- "BioTIP_attractor_allgenes"
	cat("BioTIP input: all genes h5ad (log1p, no HVG subset)\n")
}
export_m <- paste0(
	normalizePath(file.path(results_dir, export_subdir), winslash = "/", mustWork = FALSE),
	"/"
)
dir.create(export_m, recursive = TRUE, showWarnings = FALSE)
cat("All BioTIP outputs ->", export_m, "\n")

cut.preselect = 0.1   # used in optimize.sd_selection()
getNetwork.cut.fdr=0.2 # to construct RW network and extract co-expressed gene modules
getTopMCI.gene.minsize=30 # minimum number of genes in an identified CTS. We found a number between 20 and 60 works rubustly for high-throughput scRNA-seq datasets.
# getTopMCI.n.states=5 # a number setting the number of clusters to check for the significance of DNB score. This parameter can be ignored or enlarged when inputting with more clusters. We recommend to set this number as the expected number of transitions states +1/ 
MCIbottom = ifelse( input=='Scaledata', 3, 2)  # Larry log1p data uses 2 (not 0.5)

n_iterations_mci = 500
n_iterations = 100
mci_p_cutoff = 0.05   # skip gene perm when step-4 MCI local p exceeds this in step 5
gene_p_cutoff = 0.05   # skip cell perm when step-5 local gene-perm p exceeds this in step 6

if (step0) {
	seu <- load_h5ad_as_seurat(h5ad_obj)
	print_h5ad_layer_summary(seu)
	colnames(seu@meta.data)
	rownames(seu)[1:4]
	if (!r %in% colnames(seu@meta.data)) {
		stop("Cluster column '", r, "' not in h5ad meta.data. Available: ",
		     paste(colnames(seu@meta.data), collapse = ', '))
	}
	if (input == 'Scaledata' && !"scale.data" %in% SeuratObject::Layers(seu[["RNA"]])) {
		stop("input = 'Scaledata' but h5ad has no scale.data layer; use input = 'data'")
	}
	seu_rds <- larry_seu_rds_path(export_m, MuTrans_HVG_input)
	saveRDS(seu, file = seu_rds)
	cat('step0 done: saved', seu_rds, '\n')
	if ("umap" %in% Reductions(seu)) {
		umap_file <- paste0(export_m, "umap_embeddings.rds")
		saveRDS(Embeddings(seu, "umap"), file = umap_file)
		cat('step0: saved UMAP ->', umap_file, '\n')
	}
}

if (step1 || step2 || step3 || step4 || step5 || step6 || step8) {
	if (!exists('seu')) {
		seu_rds <- larry_seu_rds_path(export_m, MuTrans_HVG_input)
		if (file.exists(seu_rds)) {
			cat('Loading seu from', seu_rds, '\n')
			seu <- readRDS(seu_rds)
		} else {
			seu <- load_h5ad_as_seurat(h5ad_obj)
		}
	}
}

#################################################################################
## step1-step3:  DO not repeat;
## SCT-normalized counts and SCT-normalized residules were calcualted
## we decide to follow up with the resutls of residules
#################################################################################


# first, we will set key parameters as followings:

if(step1) {
  #for( input in c('Scaledata','data')) {
   
	## 1) transcription preselection ##########	
	
	samplesL <- split(rownames(seu@meta.data), f = seu@meta.data[[r]])
	lengths(samplesL)
	saveRDS(samplesL, file= paste0(export_m ,'samplesL_',input,'.rds'))  # save to allow reploting in local machine

	if (input == 'data') logmat <- as.matrix( LayerData(seu, assay = "RNA", layer = "data")  )
	if (input == 'Scaledata') logmat <- as.matrix( LayerData(seu, assay = "RNA", layer = "scale.data")  )
	# Warning message:
	# In asMethod(object) :
	  # sparse->dense coercion: allocating vector of size 5.2 GiB
	dim(logmat)  # 5633 218854  
	  

	set.seed(2020)
	# pb <- txtProgressBar(min = 0, max = 100, style = 3)

	# ## this function runs longer and needs larger memory !!!!!!
	testres = optimize.sd_selection(logmat[,unlist(samplesL)], 
									samplesL, B=100,   # B=100, only  29% done within 24hours !
									cutoff=cut.preselect, 
									times=.75, 
									percent=0.8,
									n_cores = n_cores, doParallel=TRUE)  #!!!!!!!!!!!!!!!!!!
	# testres = optimize.sd_selection(logmat[,unlist(samplesL)], 
									# samplesL, B=50,   # B=100, only  29% done within 24hours !
									# cutoff=cut.preselect, 
									# times=.75, 
									# percent=0.8)
								
	save(testres, file=paste0(export_m, "optimized_test_sd_selection_",input,".RData"), compress=TRUE)
	#save(testres, file=paste0(export_m, "tested_optimized_test_sd_selection_",input,".RData"), compress=TRUE)
	class(testres)  #[1] "list"
	dim(testres[[1]])  # [1] 437 selected genes x 614 sampels of the 1st elements of samplesL
    cat('step1 done','\n')
	#} # for input
}

## 2) Network Partition

if(step2){
  #for( input in c('Scaledata','data')) {

	load(file=paste0(export_m,"optimized_test_sd_selection_",input,".RData"))

	igraphL = getNetwork(testres, fdr = getNetwork.cut.fdr)
	save(igraphL, file=paste0(export_m , 'igraphL_',input,'.RData'))
	
	cluster = getCluster_methods(igraphL)

## 3) Putative Critical Transition Signals (CTSs) Identification by MCI score

	membersL = getMCI(cluster, testres, adjust.size = FALSE, fun='BioTIP')
	save(membersL, file=paste0(export_m ,'membersL_',input,'.RData'))
 #} # for input
    
	# see the real module-size distribution per cluster
	sapply(membersL$members, function(m) sort(table(m), decreasing = TRUE))

	# specifically, the single largest module in clusters 1, 3, 13
	sapply(membersL$members[c('1','3','13')], function(m) max(table(m)))

 cat('step2 done','\n')


}


# scp -p -r   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/BioTIP_leiden_r0_6/* F:\projects\TIPS\results\GSE140802_lineage_tracking\inVitro_NMtrajectory\BioTIP_leiden_r0_6\.

#################################################################
## Holly modified the online tutorial code by
## using maxMCI rather than topMCI to allowing the 2nd highest score per cluster to be considered
#################################################################

if(step3){	
	#rowData = readRDS( 'heart_rowData.rds')
#The following plot shows MCI scores of the classified cluster samples at 19 different stages.
  #for( input in c('Scaledata','data')) {

	load(paste0(export_m ,"optimized_test_sd_selection_",input,".RData"))
	load(paste0(export_m ,"membersL_",input,".RData"))

	par(mar=c(1,1,1,1))
	ymax = ifelse(input=='data',10,8)
	## setting up higher width if you have more clusters to plot
	pdf(file=paste0(export_m ,"bar_memberL_",input,".pdf"), height=5,width=35)
	plotBar_MCI(membersL, ylim=c(0,ymax), minsize = getTopMCI.gene.minsize)
    dev.off()
	
	topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=getTopMCI.gene.minsize, n=length(membersL$members)) #### updateing by xyang2
	if(any(is.na(topMCI))) topMCI = topMCI[-which(is.na(topMCI))]			#### updateing by xyang2
	#topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=getTopMCI.gene.minsize, n=length(membersL))  ## this will only reporting top 5 !!!!!
	maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], min=getTopMCI.gene.minsize, n=2)
	maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])

	CTS.Lib = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])

	tmp <- unlist(lapply(maxMCIms[['idx']][names(topMCI)], length))
	(whoistop2nd <- names(tmp[tmp==2]))
	if(length(whoistop2nd)>0) CTS.Lib = append(CTS.Lib, maxMCIms[["2topest.members"]][whoistop2nd])

	maxMCI <- maxMCI[names(CTS.Lib)[1:length(maxMCI)]]
	maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop2nd))

	#MCIbottom = ifelse( input=='Scaledata', 3, 0.5)
	print(MCIbottom)  # 3
	#CTS.Lib <- CTS.Lib[1:length(which(topMCI >= MCIbottom))]
	#CTS.Lib.Symbol <- CTS.Lib.Symbol[1:length(which(topMCI >= MCIbottom))]
	#maxMCI <- maxMCI[1:length(which(topMCI >= MCIbottom))]
	CTS.Lib <- CTS.Lib[which(maxMCI >= MCIbottom)]  ## updated by xyang2	 
	maxMCI <- maxMCI[which(maxMCI >= MCIbottom)]  ## updated by xyang2

	names(CTS.Lib) <- make.unique(names(CTS.Lib))
	names(maxMCI) <- make.unique(names(maxMCI))

	CTS.Lib.Symbol <- CTS.Lib
	## extract the ensembel gene id #!!!
	# features = seu@assays[['RNA']]@meta.features
	# if(ncol(features)>0 & 'ENSEMBL' %in% colnames(features)){	
	# 	for (i in 1:length(CTS.Lib)) {
	# 	  CTS.Lib[[i]] <- features[CTS.Lib.Symbol[[i]],]$ENSEMBL
	# 	}
	# }
	
	save(CTS.Lib, CTS.Lib.Symbol, maxMCI, file=paste0(export_m ,'CTS_Lib_',input,'.RData'))
 # }  # for input
    cat('step3 done','\n')

## printing the toopest module size per states

}

###########################################################
# The following output CTS.Lib.Symbol has clusters ordered in deheartnding MCI score.
# From now on, using scale.data / data seperately #!!!!!!!!!
###########################################################

if(step4){	
	addingMore = FALSE
	print(input)  # 'Scaledata'  #'data'
	load(paste0(export_m, 'CTS_Lib_', input, '.RData'))  # -> CTS.Lib, CTS.Lib.Symbol, maxMCI
	names(CTS.Lib.Symbol)

	if (input == 'data') logmat <- as.matrix( LayerData(seu, assay = "RNA", layer = "data")  )  # as.matrix( seu[["RNA"]]$data )
	if (input == 'Scaledata') logmat <- as.matrix( LayerData(seu, assay = "RNA", layer = "scale.data")  )  # as.matrix( seu[["RNA"]]$scale.data )
	dim(logmat)

	samplesL <- readRDS(file = paste0(export_m, 'samplesL_',input,'.rds'))  # reuse the partition CTS.Lib was built from

	FILENAME_M <- paste0(export_m, "CTS_ShrinkM_", input, ".RData")
	if (file.exists(FILENAME_M)) {
		load(FILENAME_M)
	} else {
		M = cor.shrink(logmat[,unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)
		save(M, file = FILENAME_M, compress = TRUE)
	}

	# n_iterations_mci = 500
	FILENAME_SIM <- paste0(export_m, "SimuMCI_", n_iterations_mci, "_", input, "_", cut.preselect,
	                        "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, ".RData")

	if (addingMore && file.exists(FILENAME_SIM)) {
		load(FILENAME_SIM)          # partially-filled simuMCI from a prior (killed) run
		y <- length(simuMCI) + 1    # auto-resume, no manual number to edit
	} else {
		simuMCI <- list()
		y <- 1
	}

	set.seed(2020)
	for (i in y:length(CTS.Lib)) {
		n <- length(CTS.Lib[[i]])
		simuMCI[[i]] <- simulationMCI(n, samplesL, logmat, B = n_iterations, fun = "BioTIP", M = M)
		names(simuMCI)[i] <- names(CTS.Lib)[i]
		save(simuMCI, file = FILENAME_SIM, compress = TRUE)  # checkpoint every cluster, not just at the end
	}
	cat('MCI simulation done', '\n')

	if (any(duplicated(names(maxMCI)))) {
		names(maxMCI) <- make.unique(names(maxMCI))
	} else {
		maxMCI_ClusterIDs <- lapply(names(maxMCI), function(x) unlist(strsplit(x, split='.', fixed=T))[1]) %>% unlist()
	}

	SimMCI.figWidth = 2*length(CTS.Lib)
	pdf(file = paste0(export_m, "SimuMCI_", n_iterations, "_", input, "_", cut.preselect,
	                   "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, ".pdf"),
	    width = SimMCI.figWidth, height = 5)
	par(mfrow = c(1, length(CTS.Lib)))
	for (i in 1:length(CTS.Lib)) {
		if (grepl('.', names(maxMCI)[i], fixed = TRUE)) {
			originalClusterName = unlist(strsplit(names(maxMCI)[i], split='.', fixed=T))[1]
		} else originalClusterName = names(maxMCI)[i]

		score = maxMCI[i]; names(score) = maxMCI_ClusterIDs[i]   # restored — needed for correct highlighting
		plot_MCI_Simulation(score, simuMCI[[i]], las = 2,
		                     main = paste0("C ", originalClusterName, "; ", length(CTS.Lib[[i]]), " genes",
		                                   "\n", "vs. ", n_iterations, " times of gene-permutation"),
		                     which2point = originalClusterName)
	}
	dev.off()

	pdf(file = paste0(export_m, "SimuMCI_", n_iterations, "_", input, "_", cut.preselect,
	                   "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, "_simple.pdf"),
	    width = SimMCI.figWidth, height = 5)
	par(mfrow = c(1, length(CTS.Lib)))
	for (i in 1:length(CTS.Lib)) {
		if (grepl('.', names(maxMCI)[i], fixed = TRUE)) {
			originalClusterName = unlist(strsplit(names(maxMCI)[i], split='.', fixed=T))[1]
		} else originalClusterName = names(maxMCI)[i]

		score = maxMCI[i]; names(score) = maxMCI_ClusterIDs[i]
		tmp = matrix(simuMCI[[i]][originalClusterName,], nrow = 1)
		rownames(tmp) = originalClusterName
		plot_MCI_Simulation(score, tmp, las = 2,
		                     main = paste0("C ", originalClusterName, "; ", length(CTS.Lib[[i]]), " genes",
		                                   "\n", "vs. ", n_iterations, " times of gene-permutation"),
		                     which2point = originalClusterName)
	}
	dev.off()

	cat('step4 done', '\n')
}

###########################################################
## Tipping-Points Identification and CTSs Evaluation
# In this section, we use our Ic.shrink score to evaluate our identified CTSs. 
# In the folloinwg plots, Ic* refers to the Ic.shrink score.
# This step took too long for huge dataset (e.g., 230,786 cells)
# thus has to be run the simulation for each cell cluster seperately and then manually combine them to be continue
###########################################################


tmpFolder_name <- paste0(export_m, input,"tmp")

# Check if the folder exists, and create it if not
if (!dir.exists(tmpFolder_name)) {
  dir.create(tmpFolder_name)
  cat("Folder created: ", tmpFolder_name, "\n")
} else {
  cat("Folder already exists: ", tmpFolder_name, "\n")
}


## helpers for step 5–7 significance filters (same p logic as BioTIP examples)
get_tp_state_from_cts <- function(cts_module_name) {
	if (grepl('.', cts_module_name, fixed = TRUE)) {
		unlist(strsplit(cts_module_name, split = '.', fixed = TRUE))[1]
	} else cts_module_name
}

mci_perm_p_local <- function(maxMCI_i, simuMCI_i, cts_module_name) {
	tp_state <- get_tp_state_from_cts(cts_module_name)
	obs <- as.numeric(maxMCI_i)
	if (is.na(obs) || is.null(simuMCI_i)) return(NA_real_)
	if (!tp_state %in% rownames(simuMCI_i)) return(NA_real_)
	sim_at_tp <- as.numeric(simuMCI_i[tp_state, , drop = TRUE])
	length(which(obs <= sim_at_tp)) / length(sim_at_tp)
}

gene_perm_p_local <- function(bioTIP_scores, sim_results_g, samplesL, cts_module_name) {
	originalClusterName <- get_tp_state_from_cts(cts_module_name)
	interesting <- which(names(samplesL) == originalClusterName)
	if (length(interesting) == 0) interesting <- which(rownames(sim_results_g) == originalClusterName)
	if (length(interesting) == 0) return(NA_real_)
	obs <- bioTIP_scores[originalClusterName]
	if (length(obs) == 0 || is.na(obs)) obs <- bioTIP_scores[interesting][1]
	length(which(sim_results_g[interesting, , drop = FALSE] >= obs)) / ncol(sim_results_g)
}

gene_perm_p_global <- function(bioTIP_scores, sim_results_g, samplesL, cts_module_name) {
	originalClusterName <- get_tp_state_from_cts(cts_module_name)
	interesting <- which(names(samplesL) == originalClusterName)
	if (length(interesting) == 0) interesting <- which(rownames(sim_results_g) == originalClusterName)
	if (length(interesting) == 0) return(NA_real_)
	obs <- bioTIP_scores[originalClusterName]
	if (length(obs) == 0 || is.na(obs)) obs <- bioTIP_scores[interesting][1]
	p2 <- length(which(sim_results_g >= obs)) / ncol(sim_results_g)
	p2 / nrow(sim_results_g)
}

if(step5){

	print(input) # 'Scaledata'  #'data'  #
	load(paste0(export_m ,'CTS_Lib_',input,'.RData'))
	names(CTS.Lib.Symbol)

	samplesL <- readRDS(file= paste0(export_m ,'samplesL_',input,'.rds'))

	# n_iterations = 100
	# n_iterations_mci = 500
	# mci_p_cutoff = 0.05   # skip gene perm when step-4 MCI local p exceeds this

	FILENAME_SIM <- paste0(export_m, "SimuMCI_", n_iterations_mci, "_", input, "_", cut.preselect,
	                        "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, ".RData")
	if (!file.exists(FILENAME_SIM)) stop('Run step 4 first: ', FILENAME_SIM)
	load(FILENAME_SIM)  # simuMCI

	step5_log <- data.frame(
		cts = character(), p.MCI.local = numeric(), run_gene_perm = logical(),
		stringsAsFactors = FALSE
	)
	n_skip <- 0
	n_run <- 0
	n_already <- 0
	logmat_loaded_step5 <- FALSE

	for(i in 1:length(CTS.Lib.Symbol)){
		tmpFileName <- gsub('[ /]', '_', names(CTS.Lib.Symbol)[i])
		outFile_scores <- paste0(tmpFolder_name, "/BioTIP_scores_", tmpFileName, "_", input, ".rds")
		outFile_sim    <- paste0(tmpFolder_name, "/SimResults_g_",  tmpFileName, "_", input, ".rds")
		cts_name <- names(CTS.Lib)[i]

		if (file.exists(outFile_scores) && file.exists(outFile_sim)) {
			n_already <- n_already + 1
			next
		}

		if (length(simuMCI) < i || is.null(simuMCI[[i]])) {
			stop('Step 4 MCI simulation incomplete for module ', i, ' / ', cts_name,
			     '; finish step 4 before step 5.')
		}

		p_mci <- mci_perm_p_local(maxMCI[i], simuMCI[[i]], cts_name)
		if (!is.na(p_mci) && p_mci > mci_p_cutoff) {
			cat(paste0(cts_name, ' skipped step5 (MCI p = ', round(p_mci, 4), ' > ', mci_p_cutoff, ')\n'))
			step5_log <- rbind(step5_log, data.frame(cts = cts_name, p.MCI.local = p_mci, run_gene_perm = FALSE))
			n_skip <- n_skip + 1
			next
		}

		if (!logmat_loaded_step5) {
			if (input == 'data') logmat <- as.matrix(LayerData(seu, assay = "RNA", layer = "data"))
			if (input == 'Scaledata') logmat <- as.matrix(LayerData(seu, assay = "RNA", layer = "scale.data"))
			logmat_loaded_step5 <- TRUE
		}

		CTS <- CTS.Lib.Symbol[[i]]
		n <- length(CTS)
		BioTIP_scores <- getIc(logmat[,unlist(samplesL)], samplesL, CTS, fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none' )
		SimResults_g <- simulation_Ic(n, samplesL, logmat, B=n_iterations, fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none',
				n_cores=n_cores_permu, 
				doParallel= doParallel_permu)
		saveRDS(SimResults_g, file=outFile_sim, compress=TRUE)
		saveRDS(BioTIP_scores, file=outFile_scores, compress=TRUE)
		cat(paste0(cts_name, ' gene-perm done (MCI p = ', round(p_mci, 4), ') '))
		step5_log <- rbind(step5_log, data.frame(cts = cts_name, p.MCI.local = p_mci, run_gene_perm = TRUE))
		n_run <- n_run + 1
	}

	write.csv(step5_log,
	          file = paste0(export_m, 'CTS_step5_mciP_filter_', input, '.csv'),
	          row.names = FALSE)
	cat(sprintf('step5 done: %d gene-perm run, %d already done, %d skipped (MCI p > %s)\n',
	            n_run, n_already, n_skip, mci_p_cutoff))
}

## Cell permutation (shuffle cells within trajectory states) — matches GitHub BioTIP_IC_GSE52583 demo.
## Skips modules when step-5 local gene-perm p exceeds gene_p_cutoff (default 0.05).
if(step6){

	print(input) # 'Scaledata'  #'data'  #
	load(paste0(export_m ,'CTS_Lib_',input,'.RData'))

	samplesL <- readRDS(file= paste0(export_m ,'samplesL_',input,'.rds'))

#	n_iterations = 100
#	gene_p_cutoff = 0.05   # skip cell perm when step-5 local gene-perm p exceeds this

	step6_log <- data.frame(
		cts = character(), p.IC.local = numeric(), run_sample_perm = logical(),
		stringsAsFactors = FALSE
	)
	n_skip <- 0
	n_run <- 0
	n_already <- 0
	logmat_loaded_step6 <- FALSE

	for(i in 1:length(CTS.Lib.Symbol)){
		tmpFileName <- gsub('[ /]', '_', names(CTS.Lib.Symbol)[i])
		outFile_scores <- paste0(tmpFolder_name, "/BioTIP_scores_", tmpFileName, "_", input, ".rds")
		outFile_sim_g  <- paste0(tmpFolder_name, "/SimResults_g_",  tmpFileName, "_", input, ".rds")
		outFile_sim_s  <- paste0(tmpFolder_name, "/SimResults_s_",  tmpFileName, "_", input, ".rds")
		cts_name <- names(CTS.Lib)[i]

		if (file.exists(outFile_sim_s)) {
			n_already <- n_already + 1
			next
		}
		if (!file.exists(outFile_scores) || !file.exists(outFile_sim_g)) {
			cat(paste0(cts_name, " skipped step6 (no step-5 output: failed MCI or gene-perm filter)\n"))
			n_skip <- n_skip + 1
			next
		}

		BioTIP_scores <- readRDS(outFile_scores)
		SimResults_g  <- readRDS(outFile_sim_g)
		p_local <- gene_perm_p_local(BioTIP_scores, SimResults_g, samplesL, cts_name)

		if (!is.na(p_local) && p_local > gene_p_cutoff) {
			cat(paste0(cts_name, " skipped step6 (gene-perm p = ", round(p_local, 4), " > ", gene_p_cutoff, ")\n"))
			step6_log <- rbind(step6_log, data.frame(cts = cts_name, p.IC.local = p_local, run_sample_perm = FALSE))
			n_skip <- n_skip + 1
			next
		}

		if (!logmat_loaded_step6) {
			if (input == 'data') logmat <- as.matrix(LayerData(seu, assay = "RNA", layer = "data"))
			if (input == 'Scaledata') logmat <- as.matrix(LayerData(seu, assay = "RNA", layer = "scale.data"))
			logmat_loaded_step6 <- TRUE
		}

		CTS <- CTS.Lib.Symbol[[i]]

		SimResults_s <- matrix(nrow = length(samplesL), ncol = n_iterations)
		rownames(SimResults_s) <- names(samplesL)
# to do: enable simulation_Ic_sample for 	doParallel	
# 	Error in simulation_Ic_sample(logmat, sampleNo = length(samplesL[[j]]),  :
#   unused arguments (n_cores = n_cores_permu, doParallel = doParallel_permu)
#   sparse->dense coercion: allocating vector of size 17.4 GiB
		for (j in seq_along(samplesL)) {
			SimResults_s[j, ] <- simulation_Ic_sample(
				logmat,
				sampleNo = length(samplesL[[j]]),
				genes = CTS,
				B = n_iterations,
				fun = "BioTIP",
				shrink = TRUE,
				PCC_sample.target = 'none'
				#,n_cores=n_cores_permu 
				#, doParallel= doParallel_permu
			)
		}
		saveRDS(SimResults_s, file = outFile_sim_s, compress = TRUE)
		cat(paste0(cts_name, " cell-perm done (gene-perm p = ", round(p_local, 4), ")\n"))
		step6_log <- rbind(step6_log, data.frame(cts = cts_name, p.IC.local = p_local, run_sample_perm = TRUE))
		n_run <- n_run + 1
	}

	write.csv(step6_log,
	          file = paste0(export_m, 'CTS_step6_geneP_filter_', input, '.csv'),
	          row.names = FALSE)
	cat(sprintf("step6 done: %d cell-perm runs, %d already done, %d skipped (gene-perm p > %s)\n",
	            n_run, n_already, n_skip, gene_p_cutoff))
}
 

## Highlight the predicted transition-state point on plot_Ic_Simulation (red; others black).
highlight_tp_point <- function(ic_vec, tp_idx) {
	tp_idx <- tp_idx[1]
	pts_y <- as.numeric(ic_vec)
	pts_x <- seq_along(pts_y)
	other <- setdiff(pts_x, tp_idx)
	if (length(other) > 0) points(pts_x[other], pts_y[other], pch = 19, col = "black", cex = 0.85)
	points(tp_idx, pts_y[tp_idx], pch = 19, col = "red", cex = 1.4)
}

## Step 7: assemble IC simulation objects and plot gene + cell permutation results.
## step7: dropping insignificant orphan datatmp modules (not in CTS.Lib):  
if(step7){
	print(input) #  'Scaledata' #'data' #
	samplesL <- readRDS(file= paste0(export_m ,'samplesL_',input,'.rds'))

	load(paste0(export_m ,'CTS_Lib_',input,'.RData'))

	BioTIP_scores <- SimResults_g <- SimResults_s <- list()
	files = list.files(paste0(export_m, input,'tmp/'), pattern ='BioTIP_scores')
	for(f in files){
		tmp = readRDS(paste0(tmpFolder_name, '/', f))
		id = unlist(strsplit(f, split=paste0('_',input,'.rds')))[1]
		id = unlist(strsplit(id, split='BioTIP_scores_'))[2]
		BioTIP_scores[[id]] = tmp
	}
	files = list.files(tmpFolder_name, pattern ='SimResults_g')
	for(f in files){
		tmp = readRDS(paste0(tmpFolder_name	,'/',f))
		id = unlist(strsplit(f, split= paste0('_',input,'.rds')))[1]
		id = unlist(strsplit(id, split='SimResults_g_'))[2]
		SimResults_g[[id]] = tmp
	}
	files = list.files(tmpFolder_name, pattern ='SimResults_s')
	for(f in files){
		tmp = readRDS(paste0(tmpFolder_name	,'/',f))
		id = unlist(strsplit(f, split= paste0('_',input,'.rds')))[1]
		id = unlist(strsplit(id, split='SimResults_s_'))[2]
		SimResults_s[[id]] = tmp
	}

	names(BioTIP_scores) <- gsub('/', '_', names(BioTIP_scores))
	names(BioTIP_scores) <- gsub(' ', '_', names(BioTIP_scores))
	names(SimResults_g) <- gsub('/', '_', names(BioTIP_scores))
	names(SimResults_g) <- gsub(' ', '_', names(BioTIP_scores))
	names(SimResults_s) <- gsub('/', '_', names(BioTIP_scores))
	names(SimResults_s) <- gsub(' ', '_', names(BioTIP_scores))
	names(CTS.Lib) <- gsub('/', '_', names(CTS.Lib))
	names(CTS.Lib) <- gsub(' ', '_', names(CTS.Lib))

	## Summary table for all CTS candidates (before plotting)
	# n_iterations_mci <- 500
	FILENAME_SIM <- paste0(export_m, "SimuMCI_", n_iterations_mci, "_", input, "_", cut.preselect,
	                        "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, ".RData")
	if (file.exists(FILENAME_SIM)) load(FILENAME_SIM)  # simuMCI

	cts_summary <- data.frame(
		CTS_ID = names(CTS.Lib),
		MCI_P = NA_real_,
		IC_g_p = NA_real_,
		IC_s_p = NA_real_,
		localmax = NA_character_,
		stringsAsFactors = FALSE
	)
	for (i in seq_along(CTS.Lib)) {
		cts_id <- names(CTS.Lib)[i]
		cts_key <- cts_id

		if (!is.null(simuMCI) && length(simuMCI) >= i && !is.null(simuMCI[[i]])) {
			cts_summary$MCI_P[i] <- mci_perm_p_local(maxMCI[i], simuMCI[[i]], cts_id)
		}

		bs <- BioTIP_scores[[cts_key]]
		if (!is.null(bs)) {
			sg <- SimResults_g[[cts_key]]
			cts_summary$IC_g_p[i] <- gene_perm_p_local(bs, sg, samplesL, cts_id)

			ss <- SimResults_s[[cts_key]]
			if (!is.null(ss)) {
				tp_state <- get_tp_state_from_cts(cts_id)
				interesting <- which(names(samplesL) == tp_state)
				if (length(interesting) == 0) interesting <- which(rownames(ss) == tp_state)
				obs <- bs[tp_state]
				if (length(obs) == 0 || is.na(obs)) obs <- bs[interesting[1]]
				cts_summary$IC_s_p[i] <- length(which(ss[interesting, , drop = FALSE] >= obs)) / ncol(ss)
			}

			tp_state <- get_tp_state_from_cts(cts_id)
			ic <- as.numeric(bs)
			tp_idx <- match(tp_state, names(bs))
			if (is.na(tp_idx)) tp_idx <- which(names(samplesL) == tp_state)[1]
			if (!is.na(tp_idx) && length(tp_idx) > 0) {
				cts_summary$localmax[i] <- ifelse(ic[tp_idx] >= max(ic, na.rm = TRUE), "yes", "no")
			}
		}
	}
	write.csv(cts_summary, file = paste0(export_m, 'CTS_summary_', input, '.csv'), row.names = FALSE)
	cat("CTS summary written to ", paste0(export_m, "CTS_summary_", input, ".csv"), "\n", sep = "")
	print(cts_summary)

	all(names(BioTIP_scores) == names(SimResults_g)) #TRUE

	orphan_ids <- setdiff(names(BioTIP_scores), names(CTS.Lib))
	if (length(orphan_ids) > 0) {
		cat("step7: dropping insignificant orphan datatmp modules (not in CTS.Lib): ",
		    paste(orphan_ids, collapse = ", "), "\n")
		keep_ids <- intersect(names(BioTIP_scores), names(CTS.Lib))
		BioTIP_scores <- BioTIP_scores[keep_ids]
		SimResults_g <- SimResults_g[keep_ids]
		if (length(SimResults_s) > 0) {
			SimResults_s <- SimResults_s[keep_ids]
		}
	}
	missing_ids <- setdiff(names(CTS.Lib), names(BioTIP_scores))
	if (length(missing_ids) > 0) {
		warning("step7: CTS.Lib modules without datatmp files (plots skip these): ",
		        paste(missing_ids, collapse = ", "))
	}

	BioTIP_scores <- BioTIP_scores[names(CTS.Lib)]
	SimResults_g <- SimResults_g[names(CTS.Lib)]
	if (length(SimResults_s) > 0) {
		SimResults_s <- SimResults_s[names(CTS.Lib)]
	}
	save(SimResults_g, SimResults_s, BioTIP_scores,
	     file=paste0(export_m ,"LibSF_IC_sim_BioTIP_",input,".RData"), compress=TRUE)
	cat("IC simulation objects assembled (", sum(!sapply(BioTIP_scores, is.null)),
	    " / ", length(CTS.Lib), " CTS modules with gene-perm results)\n", sep = "")

	n_with_sample_perm <- sum(!sapply(SimResults_s, is.null))
	if (n_with_sample_perm == 0) {
		warning("No SimResults_s files found; only gene-perm plots will be produced (cell perm skipped).")
	} else {
		cat(sprintf("step7: %d / %d CTS modules have cell-perm results\n",
		            n_with_sample_perm, length(CTS.Lib)))
	}

	load(paste0(export_m ,"LibSF_IC_sim_BioTIP_",input,".RData"))

	idx_g <- as.integer(which(!sapply(BioTIP_scores, is.null) & !sapply(SimResults_g, is.null)))
	if (length(idx_g) == 0) stop("No step-5 gene-perm results found for step 7 plotting.")

	n_iterations <- ncol(SimResults_g[[idx_g[1]]])
	n_iterations_s <- if (n_with_sample_perm > 0) {
		ncol(SimResults_s[[which(!sapply(SimResults_s, is.null))[1]]])
	} else NA

	p.IC.local = p.IC = p.IC.local.s = p.IC.s = array()
	for(i in idx_g){
		n = length(CTS.Lib[[i]])
	    if(grepl('.', names(BioTIP_scores[i]), fix=TRUE))  {
			originalClusterName = unlist(strsplit(names(BioTIP_scores[i]), split='.', fixed=T))[1]
		} else originalClusterName = names(BioTIP_scores[i])

		interesting = which(names(samplesL) == originalClusterName)

		p = length(which(SimResults_g[[i]][interesting,] >= BioTIP_scores[[i]][originalClusterName]))
		p.IC.local[i] = p/ncol(SimResults_g[[i]])
		p2 = length(which(SimResults_g[[i]] >= BioTIP_scores[[i]][originalClusterName]))
		p2 = p2/ncol(SimResults_g[[i]])
		p.IC[i] = p2/nrow(SimResults_g[[i]])

		if (!is.null(SimResults_s[[i]])) {
			ps = length(which(SimResults_s[[i]][interesting,] >= BioTIP_scores[[i]][originalClusterName]))
			p.IC.local.s[i] = ps/ncol(SimResults_s[[i]])
			p2s = length(which(SimResults_s[[i]] >= BioTIP_scores[[i]][originalClusterName]))
			p2s = p2s/ncol(SimResults_s[[i]])
			p.IC.s[i] = p2s/nrow(SimResults_s[[i]])
		} else {
			p.IC.local.s[i] = p.IC.s[i] = NA
		}
	}
	names(p.IC.local) = names(p.IC) = names(p.IC.local.s) = names(p.IC.s) = names(BioTIP_scores)

   if(any(grepl(' ', names(BioTIP_scores[[1]])))) {
	for(i in 1:length(BioTIP_scores)) names(BioTIP_scores[[i]]) <- gsub(' ', '_', names(BioTIP_scores[[i]]))
	names(SimResults_g) =  gsub(' ', '_', names(SimResults_g))
	names(SimResults_s) =  gsub(' ', '_', names(SimResults_s))
	for(i in 1:length(SimResults_g)) rownames(SimResults_g[[i]]) <- gsub(' ', '_', rownames(SimResults_g[[i]]) )
	for(i in 1:length(SimResults_s)) rownames(SimResults_s[[i]]) <- gsub(' ', '_', rownames(SimResults_s[[i]]) )
	}
   if(any(grepl('/', names(BioTIP_scores[[1]])))) {
	for(i in 1:length(BioTIP_scores)) names(BioTIP_scores[[i]]) <- gsub('/', '_', names(BioTIP_scores[[i]]))
	names(SimResults_g) =  gsub('/', '_', names(SimResults_g))
	names(SimResults_s) =  gsub('/', '_', names(SimResults_s))
	for(i in 1:length(SimResults_g)) rownames(SimResults_g[[i]]) <- gsub('/', '_', rownames(SimResults_g[[i]]) )
	for(i in 1:length(SimResults_s)) rownames(SimResults_s[[i]]) <- gsub('/', '_', rownames(SimResults_s[[i]]) )
	}

	## plotIcSignificance (BioTIP_update_06232025.R) uses which2point = names(BioTIP_scores)[i],
	## which must be a tipping-point attractor in rownames(SimResults), not the CTS module ID
	## (e.g. module "11.1" -> attractor "11"). Same logic as _2.pdf and Gastrulation tutorial:
	## https://github.com/xyang2uchicago/BioTIP/blob/master/tutorial/Gastrulation.html
	prep_plot_ic_significance_lists <- function(idx, bs, cts, sim) {
		idx <- as.integer(idx)
		mod_ids <- names(bs)[idx]
		tp_names <- vapply(mod_ids, get_tp_state_from_cts, character(1))
		out_bs <- lapply(idx, function(i) bs[[i]])
		out_cts <- lapply(idx, function(i) cts[[i]])
		out_sim <- lapply(idx, function(i) sim[[i]])
		names(out_bs) <- names(out_cts) <- names(out_sim) <- tp_names
		list(bs = out_bs, cts = out_cts, sim = out_sim)
	}
	plot_g <- prep_plot_ic_significance_lists(idx_g, BioTIP_scores, CTS.Lib, SimResults_g)

	plotIcSignificance(filename=paste0(export_m ,"LibSF_IC_sim_BioTIP_",input,"_g.pdf"),
				BioTIP_scores = plot_g$bs,
				CTS.candidate = plot_g$cts,
				SimResults_g = plot_g$sim,
                width=21, height=5, nc=length(idx_g))

	idx_s <- as.integer(intersect(idx_g, which(!sapply(SimResults_s, is.null))))
	if (length(idx_s) > 0) {
		plot_s <- prep_plot_ic_significance_lists(idx_s, BioTIP_scores, CTS.Lib, SimResults_s)

		plotIcSignificance(filename=paste0(export_m ,"LibSF_IC_sim_BioTIP_",input,"_s.pdf"),
					BioTIP_scores = plot_s$bs,
					CTS.candidate = plot_s$cts,
					SimResults_g = plot_s$sim,
	                width=21, height=5, nc=length(idx_s))
	}

	ylim_g = 1
	pdf(file=paste0(export_m ,"LibSF_IC_sim_BioTIP_",input,"_2.pdf"), width=10, height=8)
	  for(i in idx_g){
		n <- length(CTS.Lib[[i]])
	    if(grepl('.', names(BioTIP_scores[i]), fix=TRUE))  {
			originalClusterName = unlist(strsplit(names(BioTIP_scores[i]), split='.', fixed=T))[1]
		} else originalClusterName = names(BioTIP_scores[i])

		interesting = which(names(samplesL) == originalClusterName)
		if (length(interesting) == 0) interesting = which(rownames(SimResults_g[[i]]) == originalClusterName)

		if (!is.null(SimResults_s[[i]])) {
			ylim_s = max(c(BioTIP_scores[[i]], SimResults_s[[i]]), na.rm=TRUE) * 1.05
			par(mfrow=c(2,2))
			plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab="Ic*", ylim=c(0, ylim_g),
							main=paste0("Cluster ", names(CTS.Lib)[i], " (", n, " genes)\n",
									   "vs. ", n_iterations, " gene permutations"),
							fun="matplot", which2point= interesting)
			highlight_tp_point(BioTIP_scores[[i]], interesting)
			plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
							   main = paste0("Delta Ic* (gene perm, ", n, " genes)"), ylab=NULL,
							   xlim=range(c(BioTIP_scores[[i]][originalClusterName], SimResults_g[[i]])))
			plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_s[[i]], las = 2, ylab="Ic*", ylim=c(0, ylim_s),
							main=paste0("Cluster ", names(CTS.Lib)[i], " (", n, " genes)\n",
									   "vs. ", n_iterations_s, " cell permutations"),
							fun="matplot", which2point= interesting)
			highlight_tp_point(BioTIP_scores[[i]], interesting)
			plot_SS_Simulation(BioTIP_scores[[i]], SimResults_s[[i]],
							   main = paste0("Delta Ic* (cell perm, ", n, " genes)"), ylab=NULL,
							   xlim=range(c(BioTIP_scores[[i]][originalClusterName], SimResults_s[[i]])))
		} else {
			par(mfrow=c(1,2))
			plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], las = 2, ylab="Ic*", ylim=c(0, ylim_g),
							main=paste0("Cluster ", names(CTS.Lib)[i], " (", n, " genes)\n",
									   "vs. ", n_iterations, " gene permutations"),
							fun="matplot", which2point= interesting)
			highlight_tp_point(BioTIP_scores[[i]], interesting)
			plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]],
							   main = paste0("Delta Ic* (gene perm; cell perm skipped, ", n, " genes)"), ylab=NULL,
							   xlim=range(c(BioTIP_scores[[i]][originalClusterName], SimResults_g[[i]])))
		}
	  }
	dev.off()

	cat("step7 done\n")

}


## Step 8: UMAP panels — Cell.type.clean, attractor, entropy, land (1 x 4)
if (step8) {
	seu_rds <- larry_seu_rds_path(export_m, MuTrans_HVG_input)
	if (exists("seu") && "umap" %in% Reductions(seu)) {
		cat("step8: using seu in memory (with umap)\n")
	} else if (file.exists(seu_rds)) {
		cat("step8: loading seu from", seu_rds, "\n")
		seu <- readRDS(seu_rds)
	} else {
		cat("step8: loading seu from h5ad\n")
		seu <- load_h5ad_as_seurat(h5ad_obj)
	}
	if (!"umap" %in% Reductions(seu)) {
		cat("step8: attaching UMAP from h5ad obsm...\n")
		seu <- attach_umap_to_seurat(seu, h5ad_obj)
	}
	if (!"umap" %in% Reductions(seu)) {
		stop("No UMAP in Seurat object. Install rhdf5 (BiocManager::install('rhdf5')), ",
		     "re-run step0, then step8. h5ad: ", h5ad_obj)
	}
	cat("step8: Reductions =", paste(Reductions(seu), collapse = ", "),
	    "; UMAP dim =", paste(dim(Embeddings(seu, "umap")), collapse = "x"), "\n")

	need_cols <- c("Cell.type.clean", "attractor", "entropy", "land")
	missing_cols <- setdiff(need_cols, colnames(seu@meta.data))
	if (length(missing_cols) > 0) {
		stop("Missing meta.data columns for step8: ", paste(missing_cols, collapse = ", "),
		     "\nAvailable: ", paste(colnames(seu@meta.data), collapse = ", "))
	}

	seu$entropy <- as.numeric(as.character(seu$entropy))
	seu$land <- as.numeric(as.character(seu$land))
	seu$attractor <- as.factor(as.character(seu$attractor))
	seu$Cell.type.clean <- as.factor(as.character(seu$Cell.type.clean))

	if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for step8")
	if (!requireNamespace("patchwork", quietly = TRUE)) {
		stop("Install patchwork for step8: install.packages('patchwork')")
	}
	if (!requireNamespace("ggrastr", quietly = TRUE)) {
		stop("Install ggrastr for step8: install.packages('ggrastr')")
	}
	library(ggrastr)

	umap_w <- 12   # inches (4 panels)
	umap_h <- 3
	umap_dpi <- 300

	theme_umap <- ggplot2::theme_classic(base_size = 10) +
		ggplot2::theme(
			plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
			axis.title = ggplot2::element_blank(),
			axis.text = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			legend.text = ggplot2::element_text(size = 6),
			legend.key.size = ggplot2::unit(0.25, "cm")
		)

	p_celltype <- DimPlot(seu, reduction = "umap", group.by = "Cell.type.clean",
	                      pt.size = 1, label = TRUE, label.size = 2.5, repel = TRUE,
	                      raster = FALSE) +
		ggplot2::ggtitle("Cell.type.clean") + theme_umap
	p_celltype <- ggrastr::rasterise(p_celltype, dpi = umap_dpi)

	p_attractor <- DimPlot(seu, reduction = "umap", group.by = "attractor",
	                       pt.size = 1, label = TRUE, label.size = 2.5, repel = TRUE,
	                       raster = FALSE) +
		ggplot2::ggtitle("attractor") + theme_umap
	p_attractor <- ggrastr::rasterise(p_attractor, dpi = umap_dpi)

	p_entropy <- FeaturePlot(seu, features = "entropy", reduction = "umap",
	                         pt.size = 1, order = TRUE, raster = FALSE) +
		ggplot2::scale_colour_gradient(low = "lightgrey", high = "darkred") +
		ggplot2::ggtitle("entropy") + theme_umap +
		ggplot2::theme(legend.position = "right")
	p_entropy <- ggrastr::rasterise(p_entropy, dpi = umap_dpi)

	p_land <- FeaturePlot(seu, features = "land", reduction = "umap",
	                      pt.size = 1, order = TRUE, raster = FALSE) +
		ggplot2::scale_colour_gradient(low = "lightgrey", high = "darkred") +
		ggplot2::ggtitle("land score") + theme_umap +
		ggplot2::theme(legend.position = "right")
	p_land <- ggrastr::rasterise(p_land, dpi = umap_dpi)

	umap_pdf <- paste0(export_m, "metacell_umap_celltype_attractor_entropy_land.pdf")
	pdf(umap_pdf, width = umap_w, height = umap_h, useDingbats = FALSE)
	print(p_celltype | p_attractor | p_entropy | p_land)
	dev.off()
	cat("step8 done ->", umap_pdf, "\n")
}

# scp -p -r   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/larry_BioTIP/* F:\projects\TIPS\results\GSE140802_lineage_tracking\inVitro_NMtrajectory\larry_BioTIP\
