library(igraph)
library(dplyr)
library(data.table)

########## BEGINNING OF USER INPUT ##########

code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(ppi_path)

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
# (33 networks: HiGCTS_13 added via lfc_HiGCTS=0.5 in 11.1)

##################################################
## PageRank and EigenCentrality
##################################################
page <- lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)
ic   <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)

df_pr <- lapply(page, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
    rbindlist(., idcol = names(.))
colnames(df_pr)[1] <- "signature"

df_ic <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(df_ic)[1] <- "signature"

df <- merge(df_pr, df_ic, by = c("signature", "gene"))
df$PPI_cat <- sapply(df$signature, function(x) strsplit(x, "_")[[1]][1])
df$PPI_cat <- factor(df$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

(dim(df)[1])
# 10549

##################################################
## Strength and normalized strength
##################################################
V_strength <- lapply(graph_list, function(g) {
    str_val  <- strength(g, weights = E(g)$weight)
    norm_str <- str_val / (vcount(g) - 1)
    data.frame(strength = str_val, normalized.strength = norm_str, gene = names(str_val))
}) %>% rbindlist(., idcol = "signature")

df <- merge(df, V_strength, by = c("signature", "gene"))

##################################################
## ANND (average nearest neighbor strength; no rewiring in code_core)
##################################################
annd_obs <- list()
for (i in names(graph_list)) {
    G <- graph_list[[i]]
    V_isolated <- which(igraph::degree(G) == 0)
    G <- delete_vertices(G, V_isolated)
    if (vcount(G) > 1) {
        annd_obs[[i]] <- knn(G, weights = E(G)$weight)$knn
    } else {
        annd_obs[[i]] <- setNames(numeric(0), character(0))
    }
}

df_annd <- lapply(annd_obs, function(x) data.frame(annd = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(df_annd)[1] <- "signature"

df <- merge(df, df_annd, by = c("signature", "gene"), all.x = TRUE)

##################################################
## Rank columns  (p-values set NA — no rewiring in code_core)
##################################################
df <- df %>%
    group_by(signature) %>%
    mutate(
        rank_by_PR                  = rank(-PageRank,            na.last = "keep"),
        rank_by_strength            = rank(-strength,            na.last = "keep"),
        rank_by_normalized.strength = rank(-normalized.strength, na.last = "keep"),
        rank_by_ANND                = rank(-annd,                na.last = "keep"),
        p.PageRank                  = NA_real_,
        rank_by_p.PR                = NA_real_,
        p.annd                      = NA_real_,
        rank_by_p.ANND              = NA_real_
    ) %>%
    ungroup()

(dim(df))
# [1] 10549   16

saveRDS(df, file = "df_PAGERANK_strength_ANND.rewiring.P.rds")
write.table(df, file = "df_PAGERANK_strength_ANND.rewiring.P.tsv", sep = "\t", row.names = FALSE)
cat("Saved df_PAGERANK_strength_ANND.rewiring.P files\n")

##################################################
## Betweenness Centrality
##################################################
betweenness_list <- lapply(graph_list, function(x) betweenness(x, weights = 1 / E(x)$weight))

for (i in seq_along(betweenness_list)) {
    betweenness_list[[i]] <- data.frame(
        BetweennessCentrality = betweenness_list[[i]],
        gene                  = names(betweenness_list[[i]])
    ) %>% mutate(rank_by_BC = rank(-BetweennessCentrality, na.last = "keep"))
}

df_BC <- betweenness_list %>% rbindlist(., idcol = names(.))
colnames(df_BC)[1] <- "signature"
df_BC$PPI_cat <- sapply(df_BC$signature, function(x) strsplit(x, "_")[[1]][1])

(dim(df_BC))
# [1] 10549    5

write.table(df_BC, file = "df_betweeness.tsv", sep = "\t", row.names = FALSE)
cat("Saved df_betweeness.tsv\n")
