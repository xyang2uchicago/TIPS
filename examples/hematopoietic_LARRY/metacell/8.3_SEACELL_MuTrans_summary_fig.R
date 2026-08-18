# Cell/metacell counts per attractor and annotated lineage
# Marker heatmap supporting all 14 attractor annotations
# jaccard score comparing STRINGto IID for 4_9vs11

## --- CONFIG: run this block first when executing line-by-line interactively ---
LARRY_ROOT <- Sys.getenv(
	"LARRY_BIOTIP_ROOT",
	"F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/data/processed/larry/mutrans"
)

LARRY_BIOTIP_CODE_DIR <- Sys.getenv(
	"LARRY_BIOTIP_CODE_DIR",
	file.path("F:/projects/TIPS/source/GSE140802_lineage_tracking/")
)

results_dir = 'F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/'
setwd(results_dir)
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


# source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/BioTIP_update_06232025.R')  # enable n_cores  = 3
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

h5ad_mutrans  <- file.path(LARRY_ROOT, 'larry_seacells_mutrans.h5ad')
h5ad_allgenes <- file.path(LARRY_ROOT, 'larry_seacells_org_mutrans.h5ad')
seu_mutrans <- file.path(results_dir, 'BioTIP_attractor/seu_attractor_MuTrans_HVG.rds')

r <- 'attractor'  # h5ad obs column for BioTIP state partition

# Cell/metacell counts per attractor and annotated lineage

# UMAP reference: step 8 in 6_SEACELL_MuTrans_attractor_BioTIP.R writes
#   BioTIP_attractor/metacell_umap_celltype_attractor_entropy_land.pdf
# Summary bar/heatmap colors below are taken from the same DimPlot() calls.

# Load Seurat object (same seu as step 8)
seu <- readRDS(seu_mutrans)

seu$attractor <- factor(as.character(seu$attractor))
seu$Cell.type.clean <- factor(as.character(seu$Cell.type.clean))

## Extract discrete colours from step-8 DimPlot (before ggrastr rasterise)
tips_umap_dimplot_cols <- function(seu, group.by) {
	levs <- levels(factor(as.character(seu@meta.data[[group.by]])))
	p <- DimPlot(
		seu, reduction = "umap", group.by = group.by,
		pt.size = 1, label = FALSE, raster = FALSE
	)
	df <- p$data
	grp_col <- if (group.by %in% names(df)) group.by else "ident"
	col_col <- intersect(c("colour", "color"), names(df))[1]
	if (!is.na(col_col) && grp_col %in% names(df)) {
		cols <- vapply(
			split(df[[col_col]], df[[grp_col]]),
			function(x) as.character(x[1]),
			character(1)
		)
		return(setNames(cols[levs], levs))
	}
	pal <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$palette
	if (is.function(pal)) {
		return(setNames(pal(length(levs)), levs))
	}
	if (requireNamespace("scales", quietly = TRUE)) {
		return(setNames(scales::hue_pal()(length(levs)), levs))
	}
	setNames(
		grDevices::hcl(
			h = seq(15, 375, length.out = length(levs) + 1L)[seq_along(levs)],
			c = 100, l = 65
		),
		levs
	)
}

lineage_cols <- tips_umap_dimplot_cols(seu, "Cell.type.clean")
attractor_cols <- tips_umap_dimplot_cols(seu, "attractor")
step8_grad <- c(low = "lightgrey", high = "darkred")   # FeaturePlot entropy/land

## Figure output: save PDFs under larry_BioTIP/ (not BioTIP_attractor/)
## Input seu still from BioTIP_attractor/seu_attractor_MuTrans_HVG.rds
out_dir <- results_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Attractor order: "numerical" (0..13) or "trajectory" (MuTrans pseudotime order)
attractor_order_mode <- "trajectory"
attractor_order_trajectory <- as.character(0:12)
attractor_highlight <- c("4", "12")   # highlight on x-axis (e.g. A4, A12)

order_attractors <- function(ids, mode = attractor_order_mode) {
	ids <- as.character(ids)
	if (mode == "numerical") {
		return(as.character(sort(as.integer(ids))))
	}
	if (mode == "trajectory") {
		ord <- intersect(attractor_order_trajectory, ids)
		return(c(ord, setdiff(ids, ord)))
	}
	stop("Unknown attractor_order_mode: ", mode)
}

attractor_label <- function(id) paste0("A", id)

# Summarize cell counts per attractor
attractor_counts <- table(seu$attractor)
print(attractor_counts)
#   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
# 111 138  88  73  93  80  49  92 121 102  70 109  16  58 

# Summarize cell counts per annotated lineage (Cell.type.clean)
lineage_counts <- table(seu$Cell.type.clean)
print(lineage_counts)
    #  Baso   Ccr7_DC       Eos Erythroid  Lymphoid      Mast       Meg 
    #    93        15        20        21         1        63       107 
    #  Mono      Neut   Unddiff       pDC 
    #   116       192       571         1 

# Summarize cell counts by attractor and annotated lineage (joint)
joint_counts <- table(seu$attractor, seu$Cell.type.clean)
print(joint_counts)

attr_order <- order_attractors(rownames(joint_counts), attractor_order_mode)
joint_counts <- joint_counts[attr_order, , drop = FALSE]

if (!requireNamespace("ggplot2", quietly = TRUE)) {
	stop("Install ggplot2 for proportion plots: install.packages('ggplot2')")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
	stop("Install patchwork for aligned panels: install.packages('patchwork')")
}
suppressPackageStartupMessages({
	library(ggplot2)
	library(patchwork)
})

joint_long <- as.data.frame(joint_counts, stringsAsFactors = FALSE)
names(joint_long) <- c("attractor", "Cell.type.clean", "n")
joint_long$attractor <- factor(joint_long$attractor, levels = attr_order)
joint_long$Cell.type.clean <- factor(joint_long$Cell.type.clean, levels = names(lineage_cols))
attractor_lab_levels <- attractor_label(attr_order)
joint_long$attractor_lab <- factor(
	attractor_label(as.character(joint_long$attractor)),
	levels = attractor_lab_levels
)
joint_long <- dplyr::group_by(joint_long, attractor, attractor_lab)
joint_long <- dplyr::mutate(
	joint_long,
	total_n = sum(n),
	prop = ifelse(total_n > 0, n / total_n, 0)
)
joint_long <- dplyr::ungroup(joint_long)

attr_totals <- dplyr::distinct(joint_long, attractor, attractor_lab, total_n)
attr_totals$attractor <- factor(attr_totals$attractor, levels = attr_order)
attr_totals$attractor_lab <- factor(as.character(attr_totals$attractor_lab), levels = attractor_lab_levels)
attr_totals$highlight <- as.character(attr_totals$attractor) %in% attractor_highlight

## axis styling: colours must follow attractor_lab_levels (not alphabetical A1,A10,...)
axis_text_colors <- setNames(
	ifelse(attr_order %in% attractor_highlight, "#C00000", "black"),
	attractor_lab_levels
)
axis_text_face <- setNames(
	ifelse(attr_order %in% attractor_highlight, "bold", "plain"),
	attractor_lab_levels
)

## Main figure: 100% stacked proportions + metacell counts (aligned)
p_prop <- ggplot(joint_long, aes(x = attractor_lab, y = prop, fill = Cell.type.clean)) +
	geom_col(position = "fill", width = 0.78, color = "white", linewidth = 0.15) +
	scale_x_discrete(limits = attractor_lab_levels) +
	scale_fill_manual(values = lineage_cols, name = "Cell.type.clean") +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
	labs(x = NULL, y = "Cell-type proportion", title = "Attractor composition (100% stacked)") +
	theme_classic(base_size = 10) +
	theme(
		axis.text.x = element_text(
			angle = 90, vjust = 0.5, hjust = 1,
			color = unname(axis_text_colors[attractor_lab_levels]),
			face = unname(axis_text_face[attractor_lab_levels])
		),
		legend.position = "right",
		legend.text = element_text(size = 7),
		legend.key.size = unit(0.3, "cm"),
		plot.title = element_text(face = "bold", hjust = 0.5)
	)

p_n <- ggplot(attr_totals, aes(x = attractor_lab, y = total_n, fill = highlight)) +
	geom_col(width = 0.78, color = "grey30", linewidth = 0.2) +
	geom_text(aes(label = total_n), vjust = -0.25, size = 2.8) +
	scale_x_discrete(limits = attractor_lab_levels) +
	scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#F4A582"), guide = "none") +
	scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
	labs(x = "Attractor", y = "Metacell count", title = "Metacells per attractor") +
	theme_classic(base_size = 10) +
	theme(
		axis.text.x = element_text(
			angle = 90, vjust = 0.5, hjust = 1,
			color = unname(axis_text_colors[attractor_lab_levels]),
			face = unname(axis_text_face[attractor_lab_levels])
		),
		plot.title = element_text(face = "bold", hjust = 0.5)
	)

p_aligned <- p_prop / p_n + plot_layout(heights = c(2.2, 1))
ggsave(
	file.path(out_dir, "metacell_attractor_lineage_prop100_and_counts.pdf"),
	p_aligned, width = 12, height = 7, bg = "white"
)
message("[8.3] wrote ", file.path(out_dir, "metacell_attractor_lineage_prop100_and_counts.pdf"))



# Marker heatmap supporting all 14 attractor annotations
## using the intersection of identified DEG per attractor and the published_genes
## for Undiff, using the top 5 Undiff markers

pub_xlsx <- Sys.getenv(
	"WEINREB_TABLE_S3",
	"F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx"
)
tips_data_dir <- file.path(
	LARRY_BIOTIP_CODE_DIR, "7_data_MuTrans_TIPS_STRING", "data"
)
deg_rdata <- file.path(
	tips_data_dir, "DEG_perState_min.prop0.25_lfc0.5_FDFR0.05.rds"
)
markers_rdata <- file.path(
	tips_data_dir, "markers.up_ttest_min.prop0.25_lfc0.5_FDFR0.05.rds"
)
heatmap_n_undiff <- 5L
heatmap_n_max <- 8L
heatmap_cluster_attractors <- TRUE
heatmap_cluster_method <- "lineage"   # "lineage" = Cell.type.clean composition; "markers" = marker-gene expression

lineage_to_weinreb <- c(
	Baso = "Basophil",
	Ccr7_DC = "Dendritic cell",
	Eos = "Eosinophil",
	Lymphoid = "Lymphoid",
	Mast = "Mast cell",
	Meg = "Megakaryocyte",
	Mono = "Monocyte",
	Neut = "Neutrophil",
	pDC = "Dendritic cell"
)

match_genes_seu <- function(genes, seu_rownames) {
	genes <- unique(as.character(genes))
	genes <- genes[nzchar(genes)]
	if (!length(genes)) return(character())
	idx <- match(toupper(genes), toupper(seu_rownames))
	hits <- genes[!is.na(idx)]
	seu_rn <- seu_rownames[idx[!is.na(idx)]]
	stats::setNames(seu_rn, toupper(hits))
}

rank_deg_genes <- function(attr, markers.up, deg_attr) {
	deg_attr <- unique(toupper(as.character(deg_attr)))
	tab <- markers.up[[as.character(attr)]]
	if (is.null(tab) || !nrow(tab)) return(deg_attr)
	tab <- tab[order(-tab$summary.logFC), , drop = FALSE]
	ranked <- toupper(rownames(tab))
	c(intersect(ranked, deg_attr), setdiff(deg_attr, ranked))
}

## Seurat v5-safe mean expression matrix (genes x attractor labels)
mean_expr_by_attractor <- function(seu, genes, group_levels) {
	mat <- tryCatch(
		Seurat::GetAssayData(seu, assay = "RNA", layer = "data"),
		error = function(e) {
			Seurat::GetAssayData(seu, assay = "RNA", slot = "data")
		}
	)
	if (inherits(mat, "dgCMatrix") || inherits(mat, "Matrix")) {
		mat <- as.matrix(mat)
	}
	genes <- intersect(as.character(genes), rownames(mat))
	if (!length(genes)) {
		return(matrix(numeric(0), nrow = 0, ncol = length(group_levels),
		              dimnames = list(NULL, group_levels)))
	}
	grp <- attractor_label(as.character(seu$attractor))
	out <- matrix(
		NA_real_,
		nrow = length(genes),
		ncol = length(group_levels),
		dimnames = list(genes, group_levels)
	)
	for (lab in group_levels) {
		cells <- which(grp == lab)
		if (length(cells)) {
			out[, lab] <- rowMeans(mat[genes, cells, drop = FALSE])
		}
	}
	out
}

expr_mat_to_long <- function(avg_mat) {
	avg_mat <- as.matrix(avg_mat)
	genes <- rownames(avg_mat)
	labs <- colnames(avg_mat)
	data.frame(
		gene_seu = rep(genes, times = ncol(avg_mat)),
		attractor_lab = rep(labs, each = nrow(avg_mat)),
		expression = as.vector(avg_mat),
		stringsAsFactors = FALSE
	)
}

## Cluster attractors (for heatmap column / facet order + dendrogram)
cluster_attractors <- function(
    joint_counts,
    attr_order,
    avg_mat = NULL,
    method = c("lineage", "markers")
) {
	method <- match.arg(method)
	attr_ids <- as.character(attr_order)
	if (method == "lineage") {
		comp <- joint_counts[attr_ids, , drop = FALSE]
		comp <- comp / pmax(rowSums(comp), 1)
		rownames(comp) <- attr_ids
		d <- dist(comp, method = "euclidean")
	} else {
		if (is.null(avg_mat) || !ncol(avg_mat)) {
			stop("avg_mat required for method = 'markers'", call. = FALSE)
		}
		mat <- t(scale(t(avg_mat)))
		mat[is.na(mat)] <- 0
		d <- dist(t(mat), method = "euclidean")
	}
	hc <- stats::hclust(d, method = "ward.D2")
	list(
		hc = hc,
		order = attr_ids[hc$order],
		labels = attractor_label(attr_ids[hc$order])
	)
}

plot_attractor_dendrogram <- function(hc, attr_order_clust, axis_text_colors, axis_text_face) {
	lab_levels <- attractor_label(attr_order_clust)
	if (requireNamespace("ggdendro", quietly = TRUE)) {
		dd <- ggdendro::dendro_data(hc)
		leaf_map <- stats::setNames(lab_levels, seq_along(lab_levels))
		lab_df <- dd$labels
		lab_df$attractor_lab <- leaf_map[as.character(lab_df$label)]
		lab_df$leaf_color <- unname(axis_text_colors[lab_df$attractor_lab])
		lab_df$leaf_face <- unname(axis_text_face[lab_df$attractor_lab])
		p <- ggplot2::ggplot() +
			ggplot2::geom_segment(
				data = dd$segments,
				ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
			) +
			ggplot2::geom_text(
				data = lab_df,
				ggplot2::aes(x = x, y = y, label = attractor_lab, colour = leaf_color),
				hjust = 1.05, size = 2.5, angle = 90, fontface = lab_df$leaf_face
			) +
			ggplot2::scale_colour_identity() +
			ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
			ggplot2::labs(title = "Attractor clustering (lineage composition)", x = NULL, y = NULL) +
			ggplot2::theme_void(base_size = 8) +
			ggplot2::theme(
				plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(b = 2)),
				plot.margin = margin(t = 2, r = 12, b = 0, l = 28)
			) +
			ggplot2::coord_cartesian(clip = "off")
		return(p)
	}
	ggplot2::ggplot() +
		ggplot2::annotate(
			"text", x = 0.5, y = 0.5,
			label = "Install ggdendro for dendrogram plot\n(cluster order still applied)",
			size = 3, hjust = 0.5
		) +
		ggplot2::theme_void()
}

if (!file.exists(deg_rdata)) {
	warning("[8.3] marker heatmap skipped — run TIPS 11.1 first: ", deg_rdata, call. = FALSE)
} else if (!file.exists(pub_xlsx)) {
	warning("[8.3] marker heatmap skipped — missing Weinreb Table S3: ", pub_xlsx, call. = FALSE)
} else {
	if (!requireNamespace("readxl", quietly = TRUE)) {
		stop("Install readxl for marker heatmap: install.packages('readxl')", call. = FALSE)
	}
	DEG <- readRDS(deg_rdata)
	DEG <- lapply(DEG, function(x) unique(toupper(as.character(x))))
	markers.up <- if (file.exists(markers_rdata)) readRDS(markers_rdata) else NULL

	pub_df <- readxl::read_excel(pub_xlsx, sheet = "DGE of progenitors in vitro")
	published_genes <- split(
		toupper(as.character(pub_df[["Gene symbol"]])),
		as.character(pub_df$Lineage)
	)
	published_genes <- lapply(published_genes, unique)

	attr_dominant_lineage <- apply(joint_counts, 1, function(x) {
		colnames(joint_counts)[which.max(x)]
	})
	attr_dominant_lineage <- attr_dominant_lineage[attr_order]

	pick_attractor_markers <- function(attr, dom_lineage) {
		deg_attr <- DEG[[as.character(attr)]]
		if (is.null(deg_attr)) return(character())
		if (identical(dom_lineage, "Undiff")) {
			ranked <- rank_deg_genes(attr, markers.up, deg_attr)
			return(head(ranked, heatmap_n_undiff))
		}
		weinreb <- unname(lineage_to_weinreb[dom_lineage])
		if (is.na(weinreb) || !weinreb %in% names(published_genes)) {
			ranked <- rank_deg_genes(attr, markers.up, deg_attr)
			return(head(ranked, heatmap_n_max))
		}
		hits <- intersect(deg_attr, published_genes[[weinreb]])
		ranked <- rank_deg_genes(attr, markers.up, hits)
		if (!length(ranked)) {
			ranked <- rank_deg_genes(attr, markers.up, deg_attr)
		}
		head(ranked, heatmap_n_max)
	}

	marker_by_attr <- lapply(attr_order, function(a) {
		pick_attractor_markers(a, attr_dominant_lineage[[a]])
	})
	names(marker_by_attr) <- attr_order

	marker_tbl <- do.call(rbind, lapply(attr_order, function(a) {
		genes <- marker_by_attr[[a]]
		if (!length(genes)) return(NULL)
		data.frame(
			attractor = a,
			attractor_lab = attractor_label(a),
			dominant_lineage = attr_dominant_lineage[[a]],
			weinreb_lineage = ifelse(
				attr_dominant_lineage[[a]] == "Undiff",
				"top Undiff DEG",
				unname(lineage_to_weinreb[attr_dominant_lineage[[a]]])
			),
			gene = genes,
			stringsAsFactors = FALSE
		)
	}))
	if (is.null(marker_tbl) || !nrow(marker_tbl)) {
		warning("[8.3] marker heatmap skipped — no marker genes selected", call. = FALSE)
	} else {
		marker_tbl$gene_seu <- vapply(
			marker_tbl$gene,
			function(g) {
				h <- match_genes_seu(g, rownames(seu))
				if (length(h)) unname(h[1]) else NA_character_
			},
			character(1)
		)
		marker_tbl <- marker_tbl[!is.na(marker_tbl$gene_seu), , drop = FALSE]

		marker_tbl$gene_block <- factor(
			marker_tbl$attractor_lab,
			levels = attractor_lab_levels
		)
		gene_order <- marker_tbl$gene_seu[
			order(marker_tbl$gene_block, marker_tbl$gene)
		]
		gene_order <- unique(gene_order)

		avg_mat <- mean_expr_by_attractor(seu, gene_order, attractor_lab_levels)
		if (!ncol(avg_mat)) {
			warning("[8.3] marker heatmap skipped — empty expression matrix", call. = FALSE)
		} else {
		attr_order_heat <- attr_order
		attractor_lab_levels_heat <- attractor_lab_levels
		hc_attr <- NULL
		if (isTRUE(heatmap_cluster_attractors)) {
			cl <- cluster_attractors(
				joint_counts, attr_order, avg_mat = avg_mat,
				method = heatmap_cluster_method
			)
			attr_order_heat <- cl$order
			attractor_lab_levels_heat <- cl$labels
			hc_attr <- cl$hc
			message(
				"[8.3] heatmap attractor order (hclust): ",
				paste(attr_order_heat, collapse = ", ")
			)
		}
		axis_text_colors_heat <- setNames(
			ifelse(attr_order_heat %in% attractor_highlight, "#C00000", "black"),
			attractor_lab_levels_heat
		)
		axis_text_face_heat <- setNames(
			ifelse(attr_order_heat %in% attractor_highlight, "bold", "plain"),
			attractor_lab_levels_heat
		)

		avg_mat <- avg_mat[, attractor_lab_levels_heat, drop = FALSE]
		expr_long <- expr_mat_to_long(avg_mat)
		expr_long$attractor_lab <- factor(expr_long$attractor_lab, levels = attractor_lab_levels_heat)
		expr_long$gene_seu <- factor(expr_long$gene_seu, levels = rev(gene_order))
		expr_long$z_score <- ave(
			expr_long$expression,
			expr_long$gene_seu,
			FUN = function(x) as.numeric(scale(x))
		)
		expr_long$z_score[is.na(expr_long$z_score)] <- 0

		facet_info <- unique(marker_tbl[, c("attractor", "attractor_lab", "gene_block", "dominant_lineage")])
		facet_info <- facet_info[order(match(facet_info$attractor, attr_order_heat)), , drop = FALSE]
		facet_info$panel_title <- paste0(
			facet_info$attractor_lab, "  \u2022  ", facet_info$dominant_lineage
		)
		facet_attr_order <- intersect(attr_order_heat, facet_info$attractor)
		marker_genes_by_attr <- split(marker_tbl$gene_seu, marker_tbl$attractor)

		make_marker_panel <- function(attr_id, expr_long, show_x = FALSE) {
			lab <- attractor_label(attr_id)
			info <- facet_info[facet_info$attractor == attr_id, , drop = FALSE]
			if (!nrow(info)) return(NULL)
			lin <- info$dominant_lineage[[1]]
			genes_this <- unique(as.character(marker_genes_by_attr[[as.character(attr_id)]]))
			genes_this <- genes_this[nzchar(genes_this)]
			if (!length(genes_this)) return(NULL)
			sub <- expr_long[as.character(expr_long$gene_seu) %in% genes_this, , drop = FALSE]
			if (!nrow(sub)) return(NULL)
			genes_rev <- rev(genes_this)
			sub$gene_seu <- factor(as.character(sub$gene_seu), levels = genes_rev)
			sub$highlight_home <- as.character(sub$attractor_lab) == lab
			lin_col <- if (lin %in% names(lineage_cols)) lineage_cols[[lin]] else "grey70"

			p <- ggplot(sub, aes(x = attractor_lab, y = gene_seu, fill = z_score)) +
				geom_tile(color = "white", linewidth = 0.15)
			home_sub <- sub[sub$highlight_home, , drop = FALSE]
			if (nrow(home_sub) > 0L) {
				p <- p + geom_tile(
					data = home_sub,
					aes(x = attractor_lab, y = gene_seu),
					inherit.aes = FALSE,
					fill = NA, color = "black", linewidth = 0.55
				)
			}
			p <- p +
				scale_x_discrete(limits = attractor_lab_levels_heat) +
				scale_fill_gradient2(
					low = step8_grad["low"], mid = "white", high = step8_grad["high"],
					midpoint = 0, name = "Z-scored\navg expr"
				) +
				labs(
					title = info$panel_title[[1]],
					x = if (show_x) "Attractor (hclust order)" else NULL,
					y = NULL
				) +
				theme_classic(base_size = 9) +
				theme(
					plot.title = element_text(
						face = "bold", size = 9, hjust = 0,
						color = lin_col,
						margin = margin(b = 1, t = 2)
					),
					plot.margin = margin(t = 2, r = 4, b = 0, l = 4),
					axis.text.x = if (show_x) {
						element_text(
							angle = 90, vjust = 0.5, hjust = 1, size = 7,
							color = unname(axis_text_colors_heat[attractor_lab_levels_heat]),
							face = unname(axis_text_face_heat[attractor_lab_levels_heat])
						)
					} else {
						element_blank()
					},
					axis.text.y = element_text(size = 7),
					axis.ticks.x = if (show_x) element_line() else element_blank(),
					panel.grid = element_blank(),
					legend.position = "none"
				)
			p
		}

		panels <- lapply(seq_along(facet_attr_order), function(i) {
			make_marker_panel(
				facet_attr_order[[i]], expr_long,
				show_x = identical(facet_attr_order[[i]], facet_attr_order[[length(facet_attr_order)]])
			)
		})
		panels <- panels[!vapply(panels, is.null, logical(1))]

		p_heat_body <- patchwork::wrap_plots(panels, ncol = 1, guides = "collect") &
			theme(legend.position = "right")

		if (!is.null(hc_attr)) {
			p_dend <- plot_attractor_dendrogram(
				hc_attr, attr_order_heat,
				axis_text_colors_heat, axis_text_face_heat
			)
			p_heat <- (p_dend / p_heat_body) +
				patchwork::plot_layout(heights = c(0.14, 1))
		} else {
			p_heat <- p_heat_body
		}

		p_heat <- p_heat +
			patchwork::plot_annotation(
				title = "Attractor marker heatmap",
				subtitle = paste0(
					"Each panel lists markers selected for that attractor ",
					"(DEG \u2229 Weinreb progenitor signature; Undiff attractors: top ",
					heatmap_n_undiff, " DEG). ",
					"Panels and columns ordered by attractor dendrogram ",
					"(", heatmap_cluster_method, " similarity). ",
					"Panel title colour = dominant Cell.type.clean; ",
					"black tile outline = home attractor."
				),
				theme = theme(
					plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
					plot.subtitle = element_text(hjust = 0.5, size = 8, lineheight = 1.1)
				)
			)

		heat_h <- max(7, 0.55 * length(gene_order) + 2)
		ggsave(
			file.path(out_dir, "metacell_attractor_marker_heatmap.pdf"),
			p_heat, width = 12, height = heat_h, bg = "white"
		)
		if (!is.null(hc_attr)) {
			saveRDS(
				list(
					hc = hc_attr,
					order = attr_order_heat,
					method = heatmap_cluster_method
				),
				file.path(out_dir, "metacell_attractor_marker_hclust.rds")
			)
		}
		write.table(
			marker_tbl[, c("attractor", "attractor_lab", "dominant_lineage", "weinreb_lineage", "gene", "gene_seu")],
			file.path(out_dir, "metacell_attractor_marker_genes.tsv"),
			sep = "\t", quote = FALSE, row.names = FALSE
		)
		message("[8.3] wrote ", file.path(out_dir, "metacell_attractor_marker_heatmap.pdf"))
		message("[8.3] wrote ", file.path(out_dir, "metacell_attractor_marker_genes.tsv"))
		}
	}
}


########################
## Jaccard score: STRING vs IID (TAG 4_9vs11 only)
## Logic mirrors cardiac 1_STRING_vs_IID + 1.2 summary plot
########################

if (!requireNamespace("ggplot2", quietly = TRUE)) {
	stop("Install ggplot2 for Jaccard summary plot", call. = FALSE)
}

tips_tag <- "4_9vs11"
tips_cp <- "4"
tips_cm <- "9"   # Baso
tips_cf <- "11"  # Meg

tips_results_dir <- function(ppi) {
	file.path(
		LARRY_BIOTIP_CODE_DIR,
		paste0("7_data_MuTrans_TIPS_", ppi),
		paste0("results_core_", tips_tag)
	)
}

jaccard_score <- function(x, y) {
	x <- unique(na.omit(as.character(x)))
	y <- unique(na.omit(as.character(y)))
	u <- union(x, y)
	if (!length(u)) return(NA_real_)
	length(intersect(x, y)) / length(u)
}

genes_from_dualpull_table <- function(path) {
	if (!file.exists(path)) {
		return(list(CM = character(), CF = character()))
	}
	df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
	if (!all(c("linkage", "from", "to") %in% colnames(df))) {
		stop("Unexpected dualpull table columns: ", path, call. = FALSE)
	}
	out <- list(CM = character(), CF = character())
	for (lk in c("CM", "CF")) {
		sub <- df[df$linkage == lk, , drop = FALSE]
		if (nrow(sub)) {
			out[[lk]] <- unique(c(sub$from, sub$to))
		}
	}
	out
}

build_string_iid_gene_sets <- function(cm_bias, cf_bias) {
	dual <- intersect(cm_bias, cf_bias)
	list(
		`Dual-pull` = dual,
		`CM-lean` = setdiff(cm_bias, dual),
		`CF-lean` = setdiff(cf_bias, dual)
	)
}

dualpull_path <- function(ppi) {
	file.path(
		tips_results_dir(ppi),
		paste0("cisTarget_predicted_", tips_cp),
		paste0("PPI_graph_GRN_prediction_CTS_", tips_cp, "_dualpull_final_table.tsv")
	)
}

string_path <- dualpull_path("STRING")
iid_path <- dualpull_path("IID")
if (!file.exists(string_path) || !file.exists(iid_path)) {
	warning(
		"[8.3] Jaccard plot skipped — run STRING and IID TIPS 24.1 for ",
		tips_tag, " first",
		call. = FALSE
	)
} else {
	string_lists <- genes_from_dualpull_table(string_path)
	iid_lists <- genes_from_dualpull_table(iid_path)

	string_sets <- build_string_iid_gene_sets(string_lists$CM, string_lists$CF)
	iid_sets <- build_string_iid_gene_sets(iid_lists$CM, iid_lists$CF)

	group_levels <- c("Dual-pull", "CM-lean", "CF-lean")
plot_df <- do.call(rbind, lapply(group_levels, function(grp) {
	data.frame(
		Group = grp,
		jaccard = jaccard_score(string_sets[[grp]], iid_sets[[grp]]),
		n_STRING = length(string_sets[[grp]]),
		n_IID = length(iid_sets[[grp]]),
		stringsAsFactors = FALSE
	)
}))
rownames(plot_df) <- NULL

plot_df$Group <- factor(plot_df$Group, levels = group_levels)
plot_df$label <- ifelse(
		is.na(plot_df$jaccard),
		paste0("n/a (STRING=", plot_df$n_STRING, ", IID=", plot_df$n_IID, ")"),
		sprintf("%.3f", plot_df$jaccard)
	)

write.table(
	plot_df,
	file.path(out_dir, paste0("jaccard_STRING_vs_IID_", tips_tag, ".tsv")),
	sep = "\t", quote = FALSE, row.names = FALSE
)

plot_sub <- plot_df[!is.na(plot_df$jaccard), , drop = FALSE]
if (!nrow(plot_sub)) {
	warning(
		"[8.3] Jaccard plot skipped — no overlapping gene categories ",
		"(CM dual-pull empty for ", tips_tag, "?)",
		call. = FALSE
	)
} else {
	p_jaccard <- ggplot(plot_sub, aes(x = Group, y = jaccard)) +
		geom_point(
			color = "#4C97C9", size = 4,
			position = position_jitter(width = 0.08, height = 0)
		) +
		stat_summary(
			fun = mean, geom = "point",
			shape = 18, size = 5, color = "black"
		) +
		scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
		labs(
			title = paste0("Robustness to PPI backbone (", tips_tag, ")"),
			subtitle = paste0(
				"STRING vs IID overlap of TIPS signature genes ",
				"(CP=A", tips_cp, ", CM=A", tips_cm, " Baso, CF=A", tips_cf, " Meg)"
			),
			x = "TIPS signature gene category",
			y = "Jaccard score STRING vs IID"
		) +
		theme_classic(base_size = 14) +
		theme(
			axis.text.x = element_text(angle = 35, hjust = 1),
			plot.title = element_text(face = "bold", hjust = 0.5),
			plot.subtitle = element_text(hjust = 0.5, size = 10)
		)

	ggsave(
		file.path(out_dir, paste0("jaccard_STRING_vs_IID_", tips_tag, ".pdf")),
		p_jaccard, width = 5.5, height = 4.5, bg = "white"
	)
    
	message("[8.3] wrote ", file.path(out_dir, paste0("jaccard_STRING_vs_IID_", tips_tag, ".pdf")))
	message(
		"[8.3] Jaccard ", tips_tag, ": ",
		paste(plot_df$Group, plot_df$label, sep = "=", collapse = "; ")
	)
}
}