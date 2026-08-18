## Larry MuTrans metacell h5ad -> Seurat (BioTIP input)
## Sourced by SEACELL_MuTrans_attractor_BioTIP.R
##
## Larry h5ad layers (supp_larry.py):
##   counts = summed linear total-count-normalized metacell expression
##   data   = log1p(counts)  -> BioTIP input = 'data' (not Scaledata)
## Do NOT re-normalize (MuTrans skipped normalize_total; only log1p was applied).
##
## Midway: BiocManager::install(c("zellkonverter", "rhdf5"))


if (!exists("LARRY_H5AD_HELPERS_LOADED")) {

ensure_cell_names <- function(sce) {
	if (!is.null(colnames(sce)) && all(nzchar(colnames(sce)))) return(sce)
	meta <- as.data.frame(SummarizedExperiment::colData(sce))
	id_cols <- c("cell", "barcode", "Cell.barcode", "cell_id", "cell_index0", "seacell")
	hit <- intersect(id_cols, colnames(meta))
	if (length(hit) > 0) {
		colnames(sce) <- as.character(meta[[hit[1]]])
	} else {
		colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
	}
	sce
}

read_umap_from_h5ad <- function(path, cell_names = NULL) {
	if (!requireNamespace("rhdf5", quietly = TRUE)) return(NULL)
	candidates <- c("obsm/X_umap", "obsm/umap", "obsm/UMAP")
	emb <- NULL
	for (p in candidates) {
		if (rhdf5::h5exists(path, p)) {
			emb <- tryCatch(as.matrix(rhdf5::h5read(path, p)), error = function(e) NULL)
			if (!is.null(emb)) break
		}
	}
	if (is.null(emb) || length(emb) == 0) return(NULL)
	if (nrow(emb) == 2L && ncol(emb) > 2L) emb <- t(emb)
	if (ncol(emb) < 2L) return(NULL)
	emb <- emb[, 1:2, drop = FALSE]
	colnames(emb) <- c("UMAP_1", "UMAP_2")
	if (is.null(cell_names)) {
		for (idx_path in c("obs/_index", "obs/__index")) {
			if (rhdf5::h5exists(path, idx_path)) {
				cell_names <- rhdf5::h5read(path, idx_path)
				break
			}
		}
	}
	if (!is.null(cell_names)) {
		cell_names <- as.character(cell_names)
		if (length(cell_names) == nrow(emb)) rownames(emb) <- cell_names
	}
	emb
}

attach_umap_to_seurat <- function(seu, path) {
	if ("umap" %in% Reductions(seu)) return(seu)
	emb <- read_umap_from_h5ad(path, colnames(seu))
	if (is.null(emb)) return(seu)
	if (is.null(rownames(emb)) || !all(colnames(seu) %in% rownames(emb))) {
		if (nrow(emb) == ncol(seu)) rownames(emb) <- colnames(seu) else {
			stop("UMAP rows (", nrow(emb), ") do not match Seurat cells (", ncol(seu), ")")
		}
	} else {
		emb <- emb[colnames(seu), , drop = FALSE]
	}
	seu[["umap"]] <- SeuratObject::CreateDimReducObject(
		embeddings = emb, key = "UMAP_", assay = DefaultAssay(seu)
	)
	seu
}

load_h5ad_as_seurat <- function(path) {
	if (!file.exists(path)) stop("h5ad not found: ", path)
	if (!requireNamespace("zellkonverter", quietly = TRUE)) {
		stop("Install zellkonverter: BiocManager::install('zellkonverter')")
	}
	sce <- zellkonverter::readH5AD(path)
	assay_names <- SummarizedExperiment::assayNames(sce)
	if (length(assay_names) == 0) stop("No assays found in h5ad")
	cat("readH5AD assays:", paste(assay_names, collapse = ", "), "\n")

	sce <- ensure_cell_names(sce)
	rd_names <- SingleCellExperiment::reducedDimNames(sce)
	if (length(rd_names) > 0) {
		cat("readH5AD reducedDims:", paste(rd_names, collapse = ", "), "\n")
	} else cat("readH5AD reducedDims: (none)\n")

	meta <- as.data.frame(SummarizedExperiment::colData(sce))
	pick_assay <- function(candidates) {
		hit <- intersect(candidates, assay_names)
		if (length(hit)) hit[1] else NA_character_
	}
	counts_name <- pick_assay(c("counts", "raw_counts", "raw"))
	data_name   <- pick_assay(c("data", "logcounts", "log1p", "normalized"))
	scale_name  <- pick_assay(c("scale.data", "scaled", "scaled.data"))
	main_name   <- assay_names[1]
	as_mat <- function(nm) as.matrix(SummarizedExperiment::assay(sce, nm))

	if (!is.na(counts_name) && !is.na(data_name)) {
		seu <- CreateSeuratObject(counts = as_mat(counts_name), meta.data = meta)
		SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = as_mat(data_name))
	} else if (!is.na(counts_name)) {
		seu <- CreateSeuratObject(counts = as_mat(counts_name), meta.data = meta)
		seu <- NormalizeData(seu, verbose = FALSE)
	} else if (!is.na(data_name)) {
		mat <- as_mat(data_name)
		seu <- CreateSeuratObject(counts = mat, meta.data = meta)
		SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = mat)
	} else {
		mat <- as_mat(main_name)
		seu <- CreateSeuratObject(counts = mat, meta.data = meta)
		if (max(mat, na.rm = TRUE) <= 50) {
			SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = mat)
		} else seu <- NormalizeData(seu, verbose = FALSE)
	}
	if (inherits(seu[["RNA"]], "Assay5")) seu <- JoinLayers(seu)
	if (!is.na(scale_name)) {
		SeuratObject::SetAssayData(seu, assay = "RNA", layer = "scale.data", new.data = as_mat(scale_name))
	}
	if (!"data" %in% SeuratObject::Layers(seu[["RNA"]])) seu <- NormalizeData(seu, verbose = FALSE)

	umap_name <- rd_names[grepl("umap", rd_names, ignore.case = TRUE)]
	if (length(umap_name) == 0) umap_name <- intersect(c("X_umap", "UMAP", "umap"), rd_names)
	if (length(umap_name) > 0) {
		emb <- as.matrix(SingleCellExperiment::reducedDim(sce, umap_name[1]))
		emb <- emb[, 1:2, drop = FALSE]
		rownames(emb) <- colnames(sce)
		colnames(emb) <- c("UMAP_1", "UMAP_2")
		seu[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(seu))
		cat("UMAP loaded from reducedDim:", umap_name[1], "\n")
	} else {
		seu <- attach_umap_to_seurat(seu, path)
		if ("umap" %in% Reductions(seu)) cat("UMAP loaded from h5ad obsm via rhdf5\n")
	}
	seu
}

larry_seu_rds_path <- function(export_m, mu_trans_hvg = TRUE) {
	paste0(export_m, if (mu_trans_hvg) "seu_attractor_MuTrans_HVG.rds" else "seu_attractor_allgenes.rds")
}

print_h5ad_layer_summary <- function(seu) {
	cat("Seurat:", ncol(seu), "cells x", nrow(seu), "genes\n")
	for (layer_name in SeuratObject::Layers(seu[["RNA"]])) {
		mat <- SeuratObject::LayerData(seu, assay = "RNA", layer = layer_name)
		cat("  ", layer_name, " dim", paste(dim(mat), collapse = "x"),
		    " range", paste(signif(range(mat, na.rm = TRUE), 4), collapse = "-"), "\n")
	}
}

LARRY_H5AD_HELPERS_LOADED <- TRUE
}
