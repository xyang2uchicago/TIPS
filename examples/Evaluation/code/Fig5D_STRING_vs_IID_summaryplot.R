setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_reproducibility/STRING_vs_IID/")

## =========================
## 1) Packages
## =========================
library("ggplot2")
library("pheatmap")
library("openxlsx")
library("dplyr")



################
## summary, string vs iid per dataset 
## =========================
jac_mat = read.csv("pairwise_jaccard_lineage_lean_matrix.csv", row.names = 1)
jac_mat <- as.matrix(jac_mat)
mode(jac_mat) <- "numeric"
jac_mat_dual = read.csv("pairwise_jaccard_dual_biased_matrix.csv", row.names = 1)
jac_mat_dual <- as.matrix(jac_mat_dual)
mode(jac_mat_dual) <- "numeric"

## Helper: convert Jaccard matrix to long format
## Keep only within-dataset STRING vs IID comparisons
get_STRING_IID_pairs <- function(mat, group_name) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  nm <- rownames(mat)
  out <- list()

  for (d in c("Ibarra", "Pijuan", "Elorbany")) {

    if (group_name == "Dual-pull") {
      s <- paste0(d, "_STRING")
      i <- paste0(d, "_IID")

      if (s %in% nm && i %in% nm) {
        out[[length(out) + 1]] <- data.frame(
          jaccard = mat[s, i],
          Group = "Dual-pull",
          Dataset = d
        )
      }

    } else {
      for (lin in c("CM", "CF")) {
        s <- paste0(lin, "_", d, "_STRING")
        i <- paste0(lin, "_", d, "_IID")

        if (s %in% nm && i %in% nm) {
          out[[length(out) + 1]] <- data.frame(
            jaccard = mat[s, i],
            Group = paste0(lin, "-lean"),
            Dataset = d
          )
        }
      }
    }
  }

  do.call(rbind, out)
}


lean_df <- get_STRING_IID_pairs(jac_mat, "Lineage-lean")
dual_df <- get_STRING_IID_pairs(jac_mat_dual, "Dual-pull")

plot_df <- rbind(dual_df, lean_df)

plot_df$Group <- factor(
  plot_df$Group,
  levels = c("Dual-pull", "CM-lean", "CF-lean")
)

plot_df$Dataset <- factor(
  plot_df$Dataset,
  levels = c("Ibarra", "Pijuan", "Elorbany")
)

print(plot_df)
plot_df$Dataset <- factor(plot_df$Dataset, levels = c("Ibarra", "Pijuan", "Elorbany"))

p_STRING_IID <- ggplot(plot_df, aes(Group, jaccard)) +
  geom_point(aes(color = Dataset), size = 4,
             position = position_jitter(width = 0.12, height = 0)) +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 5, color = "black") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Robustness to PPI backbone",
    subtitle = "STRING vs IID overlap of TIPS signature gene categories",
    x = "TIPS signature gene category",
    y = "Jaccard score",
    color = "Dataset"
  ) +
  theme_classic(base_size = 14)


pdf(file = "STRING_vs_IID_dual_pull_and_lineage_lean_summary.pdf",  width = 7, height = 4)
print(p_STRING_IID)
dev.off()


#############################
## cross-dataset summary plot 
## =========================  

