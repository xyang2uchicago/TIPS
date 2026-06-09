## Questions: Do individual cells simultaneously activate both CTS-derived pulls, and how balanced are those two pulls?
## Method:
# twin-pull strength = how much both distinct modules are active in the same cell
# balance = whether the cell is CM-biased or EMT-biased


sig_map <- data.frame(
  db = c( "IbarraSoria", "Pijuan_Sala", "Elorbany"),
  signature = c("cardiac.a", "8", "CP.1"),
  cf_col = c("CF_pullcardiac.a", "CF_pull8", "CF_pullCP.1"),
  cm_col = c("CM_pullcardiac.a", "CM_pull8_wincreased_nodes", "CM_pullCP.1_wincreased_nodes"),
  source_path = c('IbarraSoria2018_E8.25_v9', 'GSE87038_Pijuan2019_v9', 'GSE175634_iPSC_CM_weighted_v9'),
  stringsAsFactors = FALSE
)

db = 'IbarraSoria'
rebuild_mat = FALSE
source('F:/projects/scRNA/source/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/24.0_acat_load_input_clean.R')

seed_TF  #"MEF2C", "GATA4", "MSX2"
names(graph_list)
#[1] "HiG_blood"                  "HiG_cardiac.b"              "HiG_cardiac.c"             
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


celltype_col  #[1] "subcelltype"
CP_cluster # 'cardiac.a'
CM_cluster # 'cardiac.c'
CF_cluster # 'extraembryonicMesoderm'
CMES_cluster # 'mixedMesoderm.a'

CTS_name # 'cardiac.a'
CTS_ID #  'CTS_cardiac.a'

lengths(CTS)
#       endothelial.b            cardiac.a  
#                  33                   37        
class(sce) # [1] "SingleCellExperiment"

dim(mat) # dim(mat)[1] 37  29

updir = "F:/projects/scRNA/results/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/GSE181346_heart_scATAC"
setwd(updir)
setwd('../validation')

library(igraph)
library(dplyr)
library(igraph)
library(openxlsx)
library(gplots)
library(ggExtra)


pull_path = 'F:/projects/TIPS/doc/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026.xlsx'), sheet=15)
dim(pull_df)  # 122  15
head(pull_df, 3)
     # Database    CTS_ID linkeage     x     y                 w_CP
# 1 IbarraSoria cardiac.a       CM GATA6  TBX5  0.14273739867444099
# 2 IbarraSoria cardiac.a       CM  CBFB GATA6  3.65114162434605E-2
# 3 IbarraSoria cardiac.a       CM  CBFB  TBX5 -1.28698510207887E-2
             # w_decendent                  delta             abs_delta
# 1  5.3873131230046403E-2 -8.8864267444394796E-2 8.8864267444394796E-2
# 2 -5.6903109139440902E-3 -4.2201727157404598E-2 4.2201727157404598E-2
# 3 -2.1015662268983001E-4    1.26596943980988E-2   1.26596943980988E-2
  # direction  status rank               TF_highConf      motif  NES
# 1  decrease changed    1                      <NA>       <NA> <NA>
# 2  decrease changed    2 CBFB (directAnnotation).  hdpi__CBFB 4.67
# 3  increase changed    3 CBFB (directAnnotation).  hdpi__CBFB 4.67
table(pull_df$Database)
   # Elorbany IbarraSoria Pijuan_Sala 
         # 79          31          12 



## for validation, we narow to validate the linege-specific (not shared) nodese !!!!!!!!!!!!!!!!!!
table(pull_df$Database, pull_df$direction)
             # decrease increase unchanged
  # Elorbany          69       10         0
  # IbarraSoria       23        8         0
  # Pijuan_Sala        8        3         1
pull_df = subset(pull_df, direction %in% c('decrease' ,'increase'))
pull_df = subset(pull_df, Database == db)
table(pull_df$direction, pull_df$linkeage)
           # CF CM
  # decrease  2  6
  # increase  0  3


cols <- setNames(
  c("orange", "blue", "lightgreen"),
  c(CP_cluster, CM_cluster, CF_cluster)
)

CTS_name # '8'
dim(sce) # [1] 10938  7240

## ----------------------------
## 1) calcualte colMean score for lineage-specfic genes
## ----------------------------

     cf_pull_name = paste0('CF_pull', CTS_name)			#!!!ADDED
	 cm_pull_name = paste0('CM_pull', CTS_name)			
	 dual_pull_name = paste0('dual_pull', CTS_name)	


   # Step 1 — Compute per-cell CF/ECM signature score 
	cf_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CF'  & CTS_ID==CTS_name)[, c('x','y')] %>%
					unlist() %>% unique()
	cm_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CM'  & CTS_ID==CTS_name)[, c('x','y')] %>%
				unlist() %>% unique()
	dual_pull_genes = intersect(cf_pull_genes_0, cm_pull_genes_0 )		
	cf_pull_genes = setdiff(cf_pull_genes_0, dual_pull_genes)  #!!!ADDED
	cm_pull_genes = setdiff(cm_pull_genes_0, dual_pull_genes)
    # toadd = case_when(db == 'Elorbany' ~ c('FGF10','LRRTM1'),  #'FGF10', 
    #             db == 'IbarraSoria' ~ NA, 
    #             db == 'Pijuan_Sala' ~ NA)  
	#  if(add_increased_nodes & !is.na(toadd))  cm_pull_genes = c(cm_pull_genes, toadd)

	colData(sce)[[cf_pull_name]] <- Matrix::colMeans(logcounts(sce)[cf_pull_genes, , drop = FALSE])
	colData(sce)[[cm_pull_name]] <- Matrix::colMeans(logcounts(sce)[cm_pull_genes, , drop = FALSE])	
	colData(sce)[[dual_pull_name]] <- Matrix::colMeans(logcounts(sce)[dual_pull_genes, , drop = FALSE])				

## for the hiPSC dataset, add the score of CM_pull that adding the increased nodes (FGF10, LRRTM1)
	if(db == 'Elorbany') {
	  cm_pull_genes = c(cm_pull_genes, 'FGF10', 'LRRTM1')
	  colData(sce)[[paste0(cm_pull_name,'_wincreased_nodes')]] <- Matrix::colMeans(logcounts(sce)[cm_pull_genes, , drop = FALSE])	
	}  
	if(db == 'Pijuan_Sala') {
	  cm_pull_genes = c(cm_pull_genes, 'HMGA2')
	  colData(sce)[[paste0(cm_pull_name,'_wincreased_nodes')]] <- Matrix::colMeans(logcounts(sce)[cm_pull_genes, , drop = FALSE])	
	}  


 
## ----------------------------
## helpers
## ----------------------------

# true 0-1 rank scaling
rank01 <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  n_ok <- sum(ok)

  if (n_ok == 0L) return(out)
  if (n_ok == 1L) {
    out[ok] <- 0.5
    return(out)
  }

  rx <- rank(x[ok], ties.method = "average")
  out[ok] <- (rx - 1) / (n_ok - 1)
  out
}

## ----------------------------
## 2) build meta and choose score columns
## ----------------------------
meta = colData(sce) %>% as.data.frame()
cf_col = sig_map$cf_col[which(sig_map$db == db)]
cm_col = sig_map$cm_col[which(sig_map$db == db)]
 

 ## ----------------------------
## 3) percentile scaling + TwinPull + Balance
## ----------------------------
meta <- meta %>%
  dplyr::mutate(
    CF   = rank01(.data[[cf_col]]),
    CM   = rank01(.data[[cm_col]]),
    Dual = rank01(.data[[dual_pull_name]]),

    # IMPORTANT: pmin high only when both are high
    TwinPull = pmin(CF, CM),

    # alternative smoother score:
    # TwinPull = sqrt(CF * CM),

    Balance = (CM - CF) / (CM + CF + 1e-6)
  )


meta$TwinIntensity <- (meta$CF + meta$CM) / 2
meta$TwinStrength   <- sqrt(meta$CF * meta$CM)

## ----------------------------
## 4) choose annotation/state column for panel K
## ----------------------------
## CHANGE THIS to your actual annotation column
celltype_col  # 'label'

stopifnot(celltype_col %in% colnames(meta))

meta_plot <- meta %>%
  dplyr::filter(!is.na(.data[[celltype_col]])) %>%
  dplyr::mutate(State = .data[[celltype_col]])

## optional: keep only selected states
## edit these labels to your real labels
# plot_states <- c(CP_cluster, CM_cluster, CF_cluster)
# meta_plot <- meta_plot %>% dplyr::filter(State %in% plot_states)
 
meta_plot <- meta_plot %>%
  dplyr::mutate(State = factor(State))


## ----------------------------
## 5) state centroids
## ----------------------------

centroids <- meta_plot %>%
  filter(!is.na(Balance), !is.na(TwinPull)) %>%
  filter(!is.na(Balance), !is.na(TwinPull)) %>%  group_by(State) %>%
  summarise(
    Balance  = median(Balance, na.rm = TRUE),
    TwinPull = median(TwinPull, na.rm = TRUE),
    TwinIntensity = median(TwinIntensity, na.rm = TRUE),
    TwinStrength  = median(TwinStrength, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

## ----------------------------
## 6) panel K: TwinPull vs Balance
## ----------------------------
y_lab = list(TwinPull = "pmin(CM, CF)", TwinIntensity = "psum(CM, CF)/2", TwinStrength = "sqrt(CF * CM)")

pK_list = list()
for(y_var in c('TwinPull', 'TwinIntensity', 'TwinStrength')) {
pK_list[[y_var]] <- ggplot(meta_plot, aes(x = Balance, y = .data[[y_var]])) +
   geom_bin2d(bins = 45) +
  scale_fill_gradient(low = "grey95", high = "grey35") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, color = "grey40") +


  geom_point(
    data = centroids,
    aes(x = Balance, y = .data[[y_var]], color = State),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x = Balance, y = .data[[y_var]], label = State, color = State),
    inherit.aes = FALSE,
    size = 3.3,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2
  ) +
 
  guides(color = guide_legend(
    override.aes = list(),
    title = "State",
    label.theme = element_text(size = 10)
  )) +
  scale_color_manual(
    values = cols,
    guide = "legend",
    na.value = "black",
    labels = function(x) {
      x <- as.character(x)
      x[x == CP_cluster] <- "CP"
      x[x == CM_cluster] <- "CM"
      x[x == CF_cluster] <- "EMT/mes-like"
      x
    }
  ) +
  annotate("text", x = -0.98, y = 0.98, label = "CF-biased", hjust = 0, size = 3.5) +
  annotate("text", x =  0.00, y = 0.98, label = "balanced",  hjust = 0.5, size = 3.5) +
  annotate("text", x =  0.98, y = 0.98, label = "CM-biased", hjust = 1, size = 3.5) +

  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
  labs(
    x = "Balance (CF/mesenchymal \u2190 0 \u2192 CM)",
    y = y_lab[[y_var]],
    title = paste0(db, " ", y_var)
  ) +
  theme_classic(base_size = 11)
}


#pK <- do.call(gridExtra::grid.arrange, pK_list)
p_final <- wrap_plots(pK_list, ncol = 1, guides = "collect") &
  theme(legend.position = "right")


## optional save
pdf(file= paste0("panelK_", db, "_TwinPull_Balance.pdf"), width = 20, height = 9)
grid.draw(p_final)
dev.off()

