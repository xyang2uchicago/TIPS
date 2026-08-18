rebuild_mat <- FALSE
source(here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_STRING", "code", "24.0_acat_load_input_clean.R"))
seed_TF # "MEF2C", "GATA4", "MSX2"
names(graph_list)
# [1] "HiG_blood"                  "HiG_cardiac.b"              "HiG_cardiac.c"
#  [4] "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"    "HiG_mixedMesoderm.a"
# [10] "HiG_mixedMesoderm.b"        "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"        "HiG_endothelial.b"
# [16] "HiG_cardiac.a"              "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"
names(DEG)
#   [1] "blood"                  "cardiac.b"              "cardiac.c"              "endothelial.a"
#  [5] "endothelial.c"          "endothelial.d"          "extraembryonicMesoderm" "mesodermProgenitors"
#  [9] "mixedMesoderm.a"        "mixedMesoderm.b"        "pharyngealMesoderm"     "presomiticMesoderm.a"
# [13] "presomiticMesoderm.b"   "somiticMesoderm"        "endothelial.b"          "cardiac.a"


celltype_col # [1] "subcelltype"
CP_cluster # 'cardiac.a'
CM_cluster # 'cardiac.c'
CF_cluster # 'extraembryonicMesoderm'
CMES_cluster # 'mixedMesoderm.a'

CTS_ID <- "cardiac.a"
CTS_name <- paste0("CTS_", CTS_ID)

lengths(CTS)
#   cardiac.a endothelial.b
#          37            33

class(sce) # [1] "SingleCellExperiment"

(updir <- getwd())
# "/Users/felixyu/Documents/GitHub/TIPS/examples/cardiac/IbarraSoria2018/results/GSE181346_heart_scATAC/ChIPseq_predicted_cardiac.a"
# create the directory if it doesn't exist
dir.create(file.path(updir, paste0("ChIPseq_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/ChIPseq_predicted_", CTS_ID))


key <- key_TFs <- "ISL1"


########################################################
##  input 6 -- data-driven --- binary glag CTS genes
### understand CTS.CP.1 ######
# for(CTS_name in c('CTS_CP.1','CTS_CP'))

# Step1 binary annotation for genes expresion, accessibility (see 24.0xxx.R
fileName <- paste0("../binary_annot_", CTS_name, "_scATAC_Maven2023_gene_ISL1_v3.tsv")
mat <- read.table(fileName, sep = "\t", header = T)
dim(mat) # [1] 37 29

print(colnames(mat))
#  [1] "CMES_hi"                          "CP_hi"                            "CM_hi"
#  [4] "CF_hi"                            "PCW6CP_access"                    "PCW8_CM_access"
#  [7] "PCW19_CM_access"                  "PCW8_CF_access"                   "PCW19_CF_access"
# [10] "PCW8_SMC_access"                  "PCW19_SMC_access"                 "PCW6_CM_access"
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"                  "iEPC_access"
# [16] "CTS_cardiac.a"                    "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"         "Maven2023_gene_ISL1_dn_T"
# [22] "Maven2023_gene_ISL1_dn_L"         "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound" "ISL1_CP_bound"                    "ISL1_CP_candidate"
# [28] "ISL1_CM_candidate"                "ISL1_CF_candidate"


####################################################
### extract subnetworks and add ISL1-bound links ###
####################################################

source(here::here("R", "celltype_specific_weight_v10.R"))

# Step 1.3) heatmap confirming key TF’s self impact by checking its targets among the CTS_CP, candidates to be the highest pagerank TFs !!!!!!!!!!!!!1
p <- heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
    key = key, coding_genes = coding_genes, TF = TF_human,
    chip_targets = TRUE, show_SMC_access = TRUE
)
#   candidate genes:  16
pdf(file = paste0("heatmap_blocked_", CTS_name, "_scATAC_chIPseq_", key, "_v3_coding_target.pdf"), height = 4)
print(p)
dev.off()





##########################################################
## -- subset of CTS[['cardiac.a']] that are
## ISL1-target
## exclusively highly expressed (HiG) in CM or CF
## further narrow the subset to be open at CP, but the bone PPIN of CTS at CP does NOT require CP accessibility, considering Isl1's pioneering role
library(BioNet)
packageVersion("BioNet") # '1.56.0'
library(igraph)
library(tibble)
library(gplots)
# for(CTS_name in c('CTS_CP.1','CTS_CP'))

{
    # Step 1.4) subPPI for targets of the key_TF
    mat <- read.table(paste0("../binary_annot_", CTS_name, "_scATAC_Maven2023_gene_ISL1_v3.tsv"), sep = "\t", header = T)
    dim(mat) # [1] 37  29

    key_target_column <- "ISL1_CP_bound"
    graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list, CTS_name, CHD,
        key = key,
        keep_selfloop = TRUE, # whether to keep the self-loop of the key
        TF_bound_column_name = key_target_column,
        TF_appendix = key,
        edge_colored_by_Maven2023_ISL1KO = TRUE,
        key_in_TFfamily = key
    )

    saveRDS(graph_TF_list, file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_v3.rds"))
    lengths(graph_TF_list)
    #      CTSHiG.CP_TF.target     CTS.CP_TF.target_HiGCM     CTS.CP_TF.target_HiGCF CTSHiG.CP_TF.target_CPopen
    #                    12                          6                          2                          3


    # Step 1.5) venn diagram comparing the keyTF targets at CP stage and two terminal lineage-specific stages
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_v3.rds"))

    plot_TF_targeted_pull_candidate(graph_TF_list, key, CTS_name, saveFigure = TRUE)
    # => PPI_graph_ISL1_GRN_prediction_CTS_cardiac.a_v3.pdf
}


##########################################################
#### step 2.1) quantify the edge weight changes and validate by Maven Isk1_KO in CP results

if (grepl(".", key, fixed = T)) warnning(paste("pls assign gene symbol to ", key))

# ## heler to paly around the layout to be beautful
make_layout <- function(g, seed) {
    set.seed(seed)
    layout_with_fr(g)
}


graph_TF_list <- readRDS(file = paste0("PPI_graph_ISL1_GRN_prediction_", CTS_name, "_v3.rds"))
names(graph_TF_list)
# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"
# [3] "CTS.CP_TF.target_HiGCF"     "CTSHiG.CP_TF.target_CPopen"
edge_attr_names(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])
# [1] "weight"         "corexp_sign"    "coexp_target"   "norm_PPI_score" "color"          "lty"


res <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CM", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CM_cluster, TF_symbol = key,
    HVG = rownames(sce)
)
cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
# ISL1 CM vcount(res[["g_CT_sub"]]): 6

if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_CM_final.rds"))
names(res)
# [1] "g_CT_sub"              "g_descendant_sub"


res <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CF", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CF_cluster, TF_symbol = key,
    HVG = rownames(sce)
)
cat(paste0(key, ' CF vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
# ISL1 CF vcount(res[["g_CT_sub"]]): 2

if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_CF_final.rds"))


(files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds"))
# [1] "PPI_graph_ISL1_GRN_prediction_CTS_cardiac.a_CF_final.rds"
# [2] "PPI_graph_ISL1_GRN_prediction_CTS_cardiac.a_CM_final.rds"


final_table <- NULL
for (pull in c("CM", "CF")) {
    res <- readRDS(file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_", pull, "_final.rds"))
    g1 <- res[["g_CT_sub"]]
    g2 <- res[["g_descendant_sub"]]

    if (vcount(g1) > 0 & vcount(g2) > 0 & vcount(g1) == vcount(g2)) {
        change_df <- edge_change_table(g1 = g1, g2 = g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)
    } else {
        cat(paste0("No edges in ", pull, "-pull subnetwork\n"))
    }


    predict <- prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5, title = paste0(pull, "_pull_subnetwork_", key))
    # => TIPS_delta_edge_reweighting_CF_pull_subnetwork.pdf

    tmp <- cbind(linkage = pull, change_df, ChIP = key)

    final_table <- rbind(final_table, tmp)
}

dim(final_table) # 11 11
####### generateing final table of the predicted subnetwork #################

write.table(final_table,
    file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
