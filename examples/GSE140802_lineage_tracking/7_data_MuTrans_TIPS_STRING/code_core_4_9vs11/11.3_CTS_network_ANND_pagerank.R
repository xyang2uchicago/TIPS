## 11.3_CTS_network_ANND_pagerank.R — Felix TIPS_core (LARRY MuTrans)

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(data.table)
})

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(ppi_path)
s <- "combined"
########## END OF USER INPUT ##########

graph_list <- readRDS(paste0(db, "_", ppi_tag, "_graph_perState_simplified_", s, "weighted.rds"))

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

V_strength <- lapply(graph_list, function(g) {
  str_val  <- strength(g, weights = E(g)$weight)
  norm_str <- str_val / (vcount(g) - 1)
  data.frame(strength = str_val, normalized.strength = norm_str, gene = names(str_val))
}) %>% rbindlist(., idcol = "signature")

df <- merge(df, V_strength, by = c("signature", "gene"))

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

df <- df %>%
  group_by(signature) %>%
  mutate(
    rank_by_PR = rank(-PageRank, na.last = "keep"),
    rank_by_strength = rank(-strength, na.last = "keep"),
    rank_by_normalized.strength = rank(-normalized.strength, na.last = "keep"),
    rank_by_ANND = rank(-annd, na.last = "keep"),
    p.PageRank = NA_real_, rank_by_p.PR = NA_real_,
    p.annd = NA_real_, rank_by_p.ANND = NA_real_
  ) %>%
  ungroup()

saveRDS(df, file = "df_PAGERANK_strength_ANND.rewiring.P.rds")
write.table(df, file = "df_PAGERANK_strength_ANND.rewiring.P.tsv", sep = "\t", row.names = FALSE)

betweenness_list <- lapply(graph_list, function(x) betweenness(x, weights = 1 / E(x)$weight))
for (i in seq_along(betweenness_list)) {
  betweenness_list[[i]] <- data.frame(
    BetweennessCentrality = betweenness_list[[i]],
    gene = names(betweenness_list[[i]])
  ) %>% mutate(rank_by_BC = rank(-BetweennessCentrality, na.last = "keep"))
}
df_BC <- betweenness_list %>% rbindlist(., idcol = names(.))
colnames(df_BC)[1] <- "signature"
df_BC$PPI_cat <- sapply(df_BC$signature, function(x) strsplit(x, "_")[[1]][1])
write.table(df_BC, file = "df_betweeness.tsv", sep = "\t", row.names = FALSE)

message("[11.3] saved PPI_weight/df_PAGERANK_strength_ANND.rewiring.P.rds")
message("[11.3] next: source('12.0_rank_by_PageRank_BC.R')")
