rebuild_mat <- FALSE
source("/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_STRING/code/24.0_acat_load_input_clean.R")

seed_TF # "ISL1"
names(graph_list)
# [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"
# [7] "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"
# [13] "HiG_19"      "HiG_7"       "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"
# [19] "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"
# [25] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"
# [31] "CTS_8"
names(DEG)
#  [1] "1"  "2"  "3"  "4"  "5"  "6"  "9"  "10" "12" "14" "17" "18" "19" "7"  "11" "15" "16" "13" "8"

celltype_col # [1] "label"
CP_cluster # '8'
CM_cluster # '17'
CF_cluster # '16'
CMES_cluster # '4'


lengths(CTS)
#   7 11 15 16  8  16.1
#  32 52 67 40  54  79

class(sce) # [1] "SingleCellExperiment"

(updir <- getwd())
# "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/results/GSE181346_heart_scATAC/ChIPseq_predicted_8"
# create the directory if it doesn't exist
dir.create(file.path(updir, paste0("ChIPseq_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/ChIPseq_predicted_", CTS_ID))

seed_TF # 'ISL1'
key <- key_TFs <- "ISL1"


########################################################
##  input 6 -- data-driven --- binary glag CTS genes
### understand CTS.CP.1 ######
# for(CTS_name in c('CTS_CP.1','CTS_CP'))

# Step1 binary annotation for genes expresion, accessibility (see 24.0xxx.R
fileName <- paste0("../binary_annot_", CTS_name, "_scATAC_Maven2023_gene_ISL1_v3.tsv")
mat <- read.table(fileName, sep = "\t", header = T)
dim(mat) # [1] 54 28

colnames(mat)
# [1] "CP_hi"                            "CM_hi"                            "CF_hi"
# [4] "PCW6CP_access"                    "PCW8_CM_access"                   "PCW19_CM_access"
# [7] "PCW8_CF_access"                   "PCW19_CF_access"                  "PCW8_SMC_access"
# [10] "PCW19_SMC_access"                 "PCW6_CM_access"                   "PCW6_CF_access"
# [13] "PCW6_SMC_access"                  "iEPC_access"                      "CTS_8"
# [16] "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"         "Maven2023_gene_ISL1_up_L"
# [19] "Maven2023_gene_ISL1_dn_E"         "Maven2023_gene_ISL1_dn_T"         "Maven2023_gene_ISL1_dn_L"
# [22] "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"   "Gao2019_gene_Isl1.iCPC_CPC.bound"
# [25] "ISL1_CP_bound"                    "ISL1_CP_candidate"                "ISL1_CM_candidate"
# [28] "ISL1_CF_candidate"


####################################################
### extract subnetworks and add ISL1-bound links ###
####################################################

source("../../../../../../R/celltype_specific_weight_v10.R")

# Step 1.3) heatmap confirming key TF’s self impact by checking its targets among the CTS_CP, candidates to be the highest pagerank TFs !!!!!!!!!!!!!1
p <- heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
    key = key, coding_genes = coding_genes, TF = TF_human,
    chip_targets = TRUE, show_SMC_access = TRUE
)
#   candidate genes:  11
pdf(file = paste0("heatmap_blocked_", CTS_name, "_scATAC_chIPseq_", key, "_v3_coding_target.pdf"), height = 4)
print(p)
dev.off()


##########################################################
## -- subset of CTS[['CP.1']] that are
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
    # mat = read.table(paste0('../../binary_annot_',CTS_name,'_scATAC_Maven2023_gene_ISL1_v3.tsv'), sep='\t', header=T)
    dim(mat) # [1] 54  29

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
    #        CTSHiG.CP_TF.target     CTS.CP_TF.target_HiGCM     CTS.CP_TF.target_HiGCF CTSHiG.CP_TF.target_CPopen
    #                   10                          4                          3                          5


    # Step 1.5) venn diagram comparing the keyTF targets at CP stage and two terminal lineage-specific stages
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_v3.rds"))

    plot_TF_targeted_pull_candidate(graph_TF_list, key, CTS_name, saveFigure = TRUE)
    # => PPI_graph_ISL1_GRN_prediction_CTS_8_v3.pdf
}


##########################################################
#### step 2.1) quantify the edge weight changes and validate by Maven Isk1_KO in CP results

if (grepl(".", key, fixed = T)) warnning(paste("pls assign gene symbol to ", key))

# ## helper to play around with the layout
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
# ISL1 CM vcount(res[["g_CT_sub"]]): 5

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
# ISL1 CF vcount(res[["g_CT_sub"]]): 3

if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_CF_final.rds"))


(files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds"))
# [1] "PPI_graph_ISL1_GRN_prediction_CTS_8_CF_final.rds"
# [2] "PPI_graph_ISL1_GRN_prediction_CTS_8_CM_final.rds"


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

dim(final_table) # 8 11
####### generateing final table of the predicted subnetwork #################

write.table(final_table,
    file = paste0("PPI_graph_", key, "_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
