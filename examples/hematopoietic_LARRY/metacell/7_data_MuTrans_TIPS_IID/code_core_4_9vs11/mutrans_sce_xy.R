## mutrans_sce_xy.R — Larry MuTrans SCE from Seurat (metacell helpers for TIPS 11.x / 24.x)

tips_sce_assay_name <- function() {
  get0("tips_sce_assay", ifnotfound = "logcounts")
}

DEFAULT_LARRY_SEURAT_RDS <- file.path(
  here::here("examples", "hematopoietic_LARRY", "data"),
  "seu_attractor_MuTrans_HVG.rds"
)

resolve_larry_seurat_path <- function(seurat_path = NULL) {
  candidates <- character()
  if (!is.null(seurat_path) && nzchar(seurat_path)) {
    candidates <- c(candidates, seurat_path)
  }
  env_path <- Sys.getenv("SEURAT_RDS", "")
  if (nzchar(env_path)) candidates <- c(candidates, env_path)
  if (exists("seurat_rds", envir = .GlobalEnv, inherits = FALSE)) {
    candidates <- c(candidates, get("seurat_rds", envir = .GlobalEnv))
  }
  if (exists("larry_biotip_dir", envir = .GlobalEnv, inherits = FALSE)) {
    candidates <- c(
      candidates,
      file.path(get("larry_biotip_dir", envir = .GlobalEnv), "seu_attractor_MuTrans_HVG.rds")
    )
  }
  candidates <- c(candidates, DEFAULT_LARRY_SEURAT_RDS)
  candidates <- unique(candidates[nzchar(candidates)])
  found <- candidates[file.exists(candidates)]
  if (length(found)) {
    return(normalizePath(found[1], winslash = "/", mustWork = TRUE))
  }
  list(candidates = unique(normalizePath(candidates, winslash = "/", mustWork = FALSE)))
}

normalize_cluster_labels <- function(x) {
  sub("\\.0+$", "", as.character(x))
}

tips_sce_assay_names <- function(sce) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required", call. = FALSE)
  }
  names(SummarizedExperiment::assays(sce))
}

assert_tips_sce_assay <- function(sce, assay = NULL) {
  if (is.null(assay)) assay <- tips_sce_assay_name()
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment is required", call. = FALSE)
  }
  avail <- tips_sce_assay_names(sce)
  if (!assay %in% avail) {
    stop(
      "SCE missing assay '", assay, "' (expected Seurat RNA data layer). ",
      "Found: ", paste(avail, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(assay)
}

get_seurat_rna_data_matrix <- function(.seu, assay = "RNA") {
  if (!assay %in% Seurat::Assays(.seu)) {
    stop("Seurat missing assay: ", assay, call. = FALSE)
  }
  mat <- tryCatch(
    Seurat::GetAssayData(.seu, assay = assay, layer = "data"),
    error = function(e) {
      tryCatch(
        Seurat::GetAssayData(.seu, assay = assay, slot = "data"),
        error = function(e2) {
          stop("Cannot read Seurat ", assay, " data layer: ", e2$message, call. = FALSE)
        }
      )
    }
  )
  if (length(mat) == 0L || max(dim(mat)) == 0L) {
    stop("Seurat ", assay, " data layer is empty", call. = FALSE)
  }
  mat
}

assert_sce_orientation <- function(sce) {
  nr <- nrow(sce)
  nc <- ncol(sce)
  rn_head <- paste(head(rownames(sce), 6L), collapse = ", ")
  if (nr < 200L) {
    stop(
      "SCE has too few rows (", nr, ") — rownames should be genes (~HVG 1000+), not metacells.\n",
      "  dim: ", nr, " x ", nc, "\n",
      "  rownames head: ", rn_head, "\n",
      "  Fix: rebuild from seu_attractor_MuTrans_HVG.rds via load_larry_sce().",
      call. = FALSE
    )
  }
  if (nc > nr) {
    stop(
      "SCE looks transposed: ncol(", nc, ") > nrow(", nr, "). Expected genes x metacells.\n",
      "  rownames head: ", rn_head,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

assert_sce_gene_overlap <- function(
    sce,
    reference_genes,
    min_overlap = 100L,
    label = "reference genes"
) {
  ref <- unique(as.character(reference_genes))
  ref <- ref[nzchar(ref)]
  ov <- sum(ref %in% rownames(sce))
  if (ov < min_overlap) {
    stop(
      "SCE rownames barely overlap ", label, " (", ov, " / ", length(ref), ").\n",
      "  sce dim: ", nrow(sce), " x ", ncol(sce), "\n",
      "  rownames head: ", paste(head(rownames(sce), 8L), collapse = ", "), "\n",
      "  ref head: ", paste(head(ref, 8L), collapse = ", "), "\n",
      "  Fix: sce <- load_larry_sce()",
      call. = FALSE
    )
  }
  message("[assert_sce_gene_overlap] ", label, ": ", ov, " / ", length(ref))
  invisible(ov)
}

validate_tips_sce <- function(sce, context = "SCE", deg_path = NULL) {
  assert_sce_orientation(sce)
  assert_tips_sce_assay(sce)
  if (is.null(deg_path) && exists("deg_rdata", inherits = TRUE)) {
    deg_path <- get("deg_rdata", inherits = TRUE)
  }
  if (!is.null(deg_path) && nzchar(deg_path) && file.exists(deg_path)) {
    deg_obj <- readRDS(deg_path)
    assert_sce_gene_overlap(
      sce, unlist(deg_obj, use.names = FALSE),
      label = paste0(context, " vs DEG")
    )
  }
  invisible(sce)
}

build_mutrans_sce <- function(seurat_path = NULL, save_path = NULL) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment is required", call. = FALSE)
  }

  resolved <- resolve_larry_seurat_path(seurat_path)
  if (is.list(resolved)) {
    stop(
      "Larry Seurat not found. Tried:\n  ",
      paste(resolved$candidates, collapse = "\n  "),
      "\nSet SEURAT_RDS or pass seurat_path = '.../seu_attractor_MuTrans_HVG.rds'",
      call. = FALSE
    )
  }
  seurat_path <- resolved
  message("[mutrans_sce] Seurat: ", seurat_path)

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to load Larry Seurat", call. = FALSE)
  }
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(SingleCellExperiment))

  assay_name <- tips_sce_assay_name()
  label_col <- if (exists("group_col", inherits = TRUE)) {
    get("group_col", inherits = TRUE)
  } else {
    "attractor"
  }

  .seu <- readRDS(seurat_path)
  if (inherits(.seu[["RNA"]], "Assay5")) .seu <- JoinLayers(.seu)
  if (!label_col %in% colnames(.seu@meta.data)) {
    stop(
      "Seurat missing group column '", label_col, "' in meta.data. ",
      "Wrong Seurat object? Columns: ", paste(colnames(.seu@meta.data), collapse = ", "),
      call. = FALSE
    )
  }

  data_mat <- get_seurat_rna_data_matrix(.seu, assay = "RNA")
  attractor_vals <- normalize_cluster_labels(.seu@meta.data[[label_col]])
  sce_obj <- SingleCellExperiment(
    assays = stats::setNames(list(data_mat), assay_name),
    colData = S4Vectors::DataFrame(.seu@meta.data, row.names = colnames(.seu))
  )
  colData(sce_obj)[[label_col]] <- attractor_vals
  sce_obj$label <- attractor_vals
  rm(.seu)
  gc()

  assert_tips_sce_assay(sce_obj, assay_name)
  validate_tips_sce(sce_obj, context = "Seurat build")
  message(
    "[mutrans_sce] assay=", assay_name, " (Seurat RNA data layer)",
    " | metacells: ", ncol(sce_obj),
    " | attractors: ", length(unique(sce_obj$label)),
    " | genes: ", nrow(sce_obj)
  )
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(sce_obj, save_path)
    message("[mutrans_sce] saved -> ", save_path)
  }
  invisible(sce_obj)
}

sync_sce_cluster_labels <- function(sce, group_col = NULL, label_col = "label") {
  if (is.null(group_col)) {
    group_col <- if (exists("group_col", inherits = TRUE)) {
      get("group_col", inherits = TRUE)
    } else {
      "attractor"
    }
  }
  cd <- SummarizedExperiment::colData(sce)
  if (!group_col %in% colnames(cd)) {
    stop(
      "sce missing group column '", group_col, "'. Columns: ",
      paste(colnames(cd), collapse = ", "),
      call. = FALSE
    )
  }
  vals <- normalize_cluster_labels(cd[[group_col]])
  SummarizedExperiment::colData(sce)[[label_col]] <- vals
  if (group_col != label_col) {
    SummarizedExperiment::colData(sce)[[group_col]] <- vals
  }
  sce
}

remap_graph_symbols_to_reference <- function(graph, reference_symbols) {
  if (!igraph::vcount(graph) || is.null(igraph::V(graph)$name)) {
    return(graph)
  }
  ref <- unique(as.character(reference_symbols))
  ref <- ref[nzchar(ref)]
  if (!length(ref)) return(graph)
  mapped_syms <- igraph::V(graph)$name
  if (all(mapped_syms %in% toupper(ref))) {
    igraph::V(graph)$name <- ref[match(toupper(mapped_syms), toupper(ref))]
  }
  graph
}

harmonize_graph_names_to_sce <- function(graph_list, sce, verbose = TRUE) {
  ref <- rownames(sce)
  ref_by_upper <- stats::setNames(ref, toupper(ref))
  n_remapped <- 0L
  for (nm in names(graph_list)) {
    g <- graph_list[[nm]]
    if (!igraph::vcount(g)) next
    vnames <- igraph::V(g)$name
    if (is.null(vnames)) next
    hits <- ref_by_upper[toupper(vnames)]
    ok <- !is.na(hits)
    if (any(ok)) {
      vnames[ok] <- unname(hits[ok])
      igraph::V(g)$name <- vnames
      graph_list[[nm]] <- g
      n_remapped <- n_remapped + sum(ok)
    }
  }
  if (verbose) {
    message("[harmonize_graph_names_to_sce] remapped ", n_remapped, " vertices to SCE casing")
  }
  graph_list
}

report_hig_sce_overlap <- function(graph_list, sce, min_n = 5L) {
  hig <- grep("^HiG_", names(graph_list), value = TRUE)
  hig <- hig[!grepl("^HiGCTS_", hig)]
  ov <- vapply(hig, function(nm) {
    length(intersect(rownames(sce), igraph::V(graph_list[[nm]])$name))
  }, numeric(1))
  message(
    "[report_hig_sce_overlap] HiG overlap with SCE: ",
    paste(sprintf("%s=%d", names(ov), as.integer(ov)), collapse = ", ")
  )
  bad <- names(ov)[ov < min_n]
  if (length(bad)) {
    warning(
      "Low SCE overlap for: ", paste(bad, collapse = ", "),
      " (need >= ", min_n, "). Check rownames(sce) vs graph symbols.",
      call. = FALSE
    )
  }
  invisible(ov)
}

sce_upper_rownames <- function(sce, assay = NULL) {
  if (is.null(assay)) assay <- tips_sce_assay_name()
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required", call. = FALSE)
  }
  rn <- toupper(rownames(sce))
  rownames(sce) <- rn
  mat_a <- SummarizedExperiment::assay(sce, assay)
  rownames(mat_a) <- rn
  SummarizedExperiment::assay(sce, assay, withDimnames = FALSE) <- mat_a
  sce
}

prepare_dualpull_case <- function(graph_list, sce, mat, DEG) {
  graph_list <- harmonize_graph_names_to_sce(graph_list, sce, verbose = FALSE)
  rownames(mat) <- toupper(rownames(mat))
  for (nm in names(graph_list)) {
    g <- graph_list[[nm]]
    if (igraph::vcount(g) > 0 && !is.null(igraph::V(g)$name)) {
      igraph::V(g)$name <- toupper(igraph::V(g)$name)
      graph_list[[nm]] <- g
    }
  }
  list(
    graph_list = graph_list,
    mat = mat,
    DEG = lapply(DEG, toupper),
    sce = sce_upper_rownames(sce)
  )
}

assert_dualpull_graphs <- function(graph_list, CM_cluster, CF_cluster, CTS_ID, CTS_name) {
  req <- c(
    paste0("HiG_", CM_cluster), paste0("HiG_", CF_cluster),
    CTS_name, paste0("HiGCTS_", CTS_ID)
  )
  missing <- setdiff(req, names(graph_list))
  if (length(missing)) {
    stop(
      "graph_list missing for dual-pull: ", paste(missing, collapse = ", "),
      "\n  Re-run 11.2.0 (HiG_", CM_cluster, " / HiG_", CF_cluster,
      " must be in weighted RDS).",
      call. = FALSE
    )
  }
  empty <- req[vapply(req, function(nm) igraph::vcount(graph_list[[nm]]), numeric(1)) == 0]
  if (length(empty)) {
    stop("Empty dual-pull graphs: ", paste(empty, collapse = ", "), call. = FALSE)
  }
  invisible(req)
}
