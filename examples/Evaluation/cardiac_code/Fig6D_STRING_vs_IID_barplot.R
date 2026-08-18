#  source('F:/projects/scRNA/source/cardiac_CTS_GRN/evaluation_reproducibility/1_STRING_vs_IID_jaccard_heatmap.R')

setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_reproducibility/STRING_vs_IID/")

## =========================
## 1) Packages
## =========================
library("ggplot2")
library("pheatmap")
library("openxlsx")
library("dplyr")


########################################################
## bar plot compare STRING vs IID for #of nodes and # of edges
########################################################
df = data.frame(Database = c("STRING_hs","STRING_ms","IID_hs_ms"),
                Species = c("hs", "ms", 'hs_ms'),
                Nodes = c(19699, 22048, 17182 ),
                Edges = c(7533072, 6859804, 1517715))

library(ggplot2)
library(gridExtra)

# Reshape data to long format for easier plotting
df_long <- reshape2::melt(df,
                          id.vars = c("Database", "Species"),
                          measure.vars = c("Nodes", "Edges"),
                          variable.name = "Metric",
                          value.name = "Count")

# Plot for Nodes
p_nodes <- ggplot(subset(df_long, Metric == "Nodes"), aes(x = Database, y = Count, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.8)) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width=0.8)) +
  labs(x = "Database", y = "# of Nodes", title = "Number of Nodes") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "lightblue", "orange"))

# Plot for Edges
p_edges <- ggplot(subset(df_long, Metric == "Edges"), aes(x = Database, y = Count, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.8)) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width=0.8)) +
  labs(x = "Database", y = "# of Edges", title = "Number of Edges") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "lightblue", "orange"))

pdf(file = "bar_plot_STRING_vs_IID.pdf")  # !!! fig 5D left
print(grid.arrange(p_nodes, p_edges, ncol = 2))
dev.off()

