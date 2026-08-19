## 11.1.1_check_vertex_duplication.R — STRING duplicate symbol cleanup (Felix / lung pattern)
##
## Writes correct_n_edges_HiG_STRING2.14.0.rds for 11.2.0 when HiG graphs have
## multiple STRING proteins mapped to the same gene symbol.
## Also writes correct_n_edges_HiG_STRING2.14.0_auto_pick_review.txt (KEPT/DROPPED).
##
## After one successful biomaRt run, writes data_dir/unique_STRING_mapping_correction.txt.
## Re-runs read that table and skip biomaRt (local STRINGdb only).

suppressPackageStartupMessages({
  library(igraph)
  library(STRINGdb)
})

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "single_cell"), "code_core_11_10vs17"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
########## END OF USER INPUT ##########

string_mapping_correction_file <- file.path(data_dir, "unique_STRING_mapping_correction.txt")

read_mapping_correction_cache <- function(path) {
  tab <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  required <- c(
    "graph_id", "input_symbol", "STRING_id", "gene_id",
    "gene_symbol_biomart", "decision", "pick_reason"
  )
  if (!all(required %in% names(tab))) {
    stop(
      "[11.1.1] ", basename(path), " missing columns: ",
      paste(setdiff(required, names(tab)), collapse = ", "),
      call. = FALSE
    )
  }
  tab
}

apply_cached_dup_resolution <- function(
  mapped, graph_id, input_symbol, mapping_cache, dup_review = NULL
) {
  dn <- input_symbol
  dup_entries <- mapped[mapped$symbol %in% dn, , drop = FALSE]
  sub <- mapping_cache[
    mapping_cache$graph_id == graph_id &
      toupper(mapping_cache$input_symbol) == toupper(dn),
    ,
    drop = FALSE
  ]
  if (!nrow(sub)) {
    warning(
      "[11.1.1] no cached rows for ", graph_id, " / ", dn,
      " — keeping all duplicate STRING_ids"
    )
    return(list(mapped = mapped, dup_review = dup_review))
  }
  missing_ids <- setdiff(dup_entries$STRING_id, sub$STRING_id)
  if (length(missing_ids)) {
    warning(
      "[11.1.1] cache missing STRING_id(s) for ", graph_id, " / ", dn, ": ",
      paste(missing_ids, collapse = ", ")
    )
  }
  drop_ids <- sub$STRING_id[sub$decision == "DROPPED"]
  if (length(drop_ids)) {
    mapped <- mapped[!(mapped$STRING_id %in% drop_ids), , drop = FALSE]
  }
  review_block <- sub[sub$STRING_id %in% dup_entries$STRING_id, required_review_cols(sub), drop = FALSE]
  dup_review <- if (is.null(dup_review)) review_block else rbind(dup_review, review_block)
  list(mapped = mapped, dup_review = dup_review)
}

required_review_cols <- function(tab) {
  c(
    "graph_id", "input_symbol", "STRING_id", "gene_id",
    "gene_symbol_biomart", "decision", "pick_reason"
  )
}

if (!file.exists(graph_notsimp_file)) {
  stop("Run 11.1 first — missing: ", graph_notsimp_file)
}
if (!file.exists(deg_rdata)) {
  stop("Missing DEG: ", deg_rdata)
}

graph_list <- readRDS(graph_notsimp_file)
graph_simp <- lapply(graph_list, igraph::simplify, edge.attr.comb = "max")

graphs_with_duplicates <- sapply(graph_simp, function(g) {
  nm <- V(g)$name
  if (is.null(nm)) nm <- as.character(V(g))
  any(duplicated(nm))
})
dup_graphs <- names(which(graphs_with_duplicates))
dup_hig <- dup_graphs[grepl("^HiG_", dup_graphs)]

if (!length(dup_hig)) {
  message("[11.1.1] no duplicate vertex names after simplify — skipping correct_n_edges")
  if (file.exists(correct_edges_file)) {
    message("[11.1.1] existing ", basename(correct_edges_file), " left in place")
  }
  invisible(NULL)
} else {
  message("[11.1.1] duplicate vertices in: ", paste(dup_hig, collapse = ", "))

  DEG <- readRDS(deg_rdata)
  DEG <- lapply(DEG, function(x) data.frame(symbol = x, stringsAsFactors = FALSE))

  string_db <- STRINGdb$new(
    version = "12.0", species = db_species,
    score_threshold = 200, network_type = "full",
    input_directory = string_ppin_dir
  )

  use_mapping_cache <- file.exists(string_mapping_correction_file)
  mapping_cache <- NULL
  if (use_mapping_cache) {
    message("[11.1.1] using cached correction table (skip biomaRt): ", string_mapping_correction_file)
    mapping_cache <- read_mapping_correction_cache(string_mapping_correction_file)
  }

  open_mmusculus_mart <- function() {
    hosts <- c("www.ensembl.org", "useast.ensembl.org", "asia.ensembl.org")
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(300, old_timeout))
    last_err <- NULL
    for (host in hosts) {
      mart <- tryCatch({
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = "mmusculus_gene_ensembl",
          host = host
        )
      }, error = function(e) {
        last_err <<- conditionMessage(e)
        NULL
      })
      if (!is.null(mart)) {
        message("[11.1.1] biomaRt host: ", host)
        return(mart)
      }
      message("[11.1.1] biomaRt host failed (", host, "): ", last_err)
    }
    stop("[11.1.1] biomaRt failed on all Ensembl mirrors", call. = FALSE)
  }

  # Shared across metacell/ and single_cell/ arms (see hematopoietic_LARRY/data/)
  peptide_cache_file <- file.path(dirname(wd), "data", "ensembl_peptide_to_gene_mmusculus.rds")
  load_peptide_gene_map <- function(peptide_ids) {
    peptide_ids <- unique(peptide_ids[nzchar(peptide_ids)])
    cached <- if (file.exists(peptide_cache_file)) {
      readRDS(peptide_cache_file)
    } else {
      data.frame(
        ensembl_peptide_id = character(),
        ensembl_gene_id = character(),
        external_gene_name = character(),
        stringsAsFactors = FALSE
      )
    }
    missing <- setdiff(peptide_ids, cached$ensembl_peptide_id)
    if (length(missing)) {
      message("[11.1.1] biomaRt lookup for ", length(missing), " peptide id(s)")
      fetched <- tryCatch({
        mart <- open_mmusculus_mart()
        biomaRt::getBM(
          attributes = c("ensembl_peptide_id", "ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_peptide_id",
          values = missing,
          mart = mart
        )
      }, error = function(e) {
        message(
          "[11.1.1] biomaRt fetch failed — using edge-count / first-row fallback only: ",
          conditionMessage(e)
        )
        data.frame(
          ensembl_peptide_id = character(),
          ensembl_gene_id = character(),
          external_gene_name = character(),
          stringsAsFactors = FALSE
        )
      })
      if (nrow(fetched)) {
        cached <- unique(rbind(cached, fetched))
        dir.create(dirname(peptide_cache_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(cached, peptide_cache_file)
        message("[11.1.1] cached ", nrow(cached), " peptide mappings -> ", peptide_cache_file)
      }
    }
    cached
  }

  map_protein_to_gene <- function(string_id, peptide_map) {
    peptide_id <- sub("10090\\.", "", string_id)
    hit <- peptide_map[match(peptide_id, peptide_map$ensembl_peptide_id), , drop = FALSE]
    if (!nrow(hit)) {
      data.frame(STRING_id = string_id, gene_id = NA, gene_symbol = NA, stringsAsFactors = FALSE)
    } else {
      data.frame(
        STRING_id = string_id,
        gene_id = hit$ensembl_gene_id[1],
        gene_symbol = hit$external_gene_name[1],
        stringsAsFactors = FALSE
      )
    }
  }

  correct_n_edges <- NULL
  dup_review <- NULL
  hig_clusters <- sub("^HiG_", "", dup_hig)
  peptide_map <- NULL

  if (!use_mapping_cache) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop(
        "[11.1.1] biomaRt required when ", basename(string_mapping_correction_file),
        " is missing — install.packages('biomaRt')",
        call. = FALSE
      )
    }
    dup_string_ids <- character()
    for (i in hig_clusters) {
      if (!i %in% names(DEG)) next
      mapped_pre <- string_db$map(DEG[[i]], "symbol", removeUnmappedRows = TRUE)
      dup_name_pre <- unique(mapped_pre$symbol[duplicated(mapped_pre$symbol)])
      if (!length(dup_name_pre)) next
      for (dn in dup_name_pre) {
        dup_string_ids <- c(
          dup_string_ids,
          mapped_pre$STRING_id[mapped_pre$symbol %in% dn]
        )
      }
    }
    peptide_map <- load_peptide_gene_map(unique(sub("10090\\.", "", dup_string_ids)))
  }

  for (i in hig_clusters) {
    if (!i %in% names(DEG)) next
    diff_exp <- DEG[[i]]
    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)
    dup_name <- unique(mapped$symbol[duplicated(mapped$symbol)])
    if (!length(dup_name)) next

    message("[11.1.1] HiG_", i, " duplicated symbols: ", paste(dup_name, collapse = ", "))
    graph_id <- paste0("HiG_", i)

    for (dn in dup_name) {
      if (use_mapping_cache) {
        res <- apply_cached_dup_resolution(
          mapped, graph_id = graph_id, input_symbol = dn,
          mapping_cache = mapping_cache, dup_review = dup_review
        )
        mapped <- res$mapped
        dup_review <- res$dup_review
      } else {
        dup_entries <- mapped[mapped$symbol %in% dn, , drop = FALSE]
        annotation <- do.call(
          rbind,
          lapply(dup_entries$STRING_id, map_protein_to_gene, peptide_map = peptide_map)
        )
        annotated_dup <- cbind(dup_entries, annotation[, c("gene_id", "gene_symbol")])

        keep_idx <- which(toupper(annotated_dup$gene_symbol) == toupper(dn))
        pick_reason <- "biomaRt_symbol_match"
        if (!length(keep_idx)) {
          keep_idx <- 1L
          pick_reason <- "fallback_first_row"
        }
        kept_id <- annotated_dup$STRING_id[keep_idx[1]]
        drop_ids <- setdiff(annotated_dup$STRING_id, kept_id)

        review_block <- data.frame(
          graph_id = graph_id,
          input_symbol = dn,
          STRING_id = annotated_dup$STRING_id,
          gene_id = annotated_dup$gene_id,
          gene_symbol_biomart = annotated_dup$gene_symbol,
          decision = ifelse(annotated_dup$STRING_id %in% kept_id, "KEPT", "DROPPED"),
          pick_reason = ifelse(
            annotated_dup$STRING_id %in% kept_id, pick_reason, "not_selected"
          ),
          stringsAsFactors = FALSE
        )
        dup_review <- if (is.null(dup_review)) review_block else rbind(dup_review, review_block)

        if (length(drop_ids)) {
          mapped <- mapped[!(mapped$STRING_id %in% drop_ids), , drop = FALSE]
        }
      }
    }

    hits <- mapped$STRING_id
    graph <- string_db$get_subnetwork(hits)
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol

    for (dn in dup_name) {
      if (!dn %in% V(graph)$name) next
      edges <- incident(graph, dn, mode = "all")
      n <- get.edgelist(graph)[edges, , drop = FALSE]
      res <- data.frame(
        graph_id = graph_id,
        names = dn,
        n_edge = nrow(n),
        STRING_id = subset(mapped, symbol == dn)$STRING_id,
        stringsAsFactors = FALSE
      )
      correct_n_edges <- if (is.null(correct_n_edges)) res else rbind(correct_n_edges, res)
    }
  }

  if (is.null(correct_n_edges) || !nrow(correct_n_edges)) {
    stop(
      "[11.1.1] duplicates detected but could not build correct_n_edges — manual fix needed ",
      "(see lung code_core_lung2026/11.1.1_check_vertex_duplication.R)",
      call. = FALSE
    )
  }

  graph_list <- readRDS(graph_notsimp_file)
  correct_n_edges$vertex_index_to_remove <- vector("list", nrow(correct_n_edges))

  for (g_name in dup_hig) {
    g <- graph_list[[g_name]]
    vertex_names <- V(g)$name
    duplicated_names <- unique(vertex_names[duplicated(vertex_names)])
    if (!length(duplicated_names)) next

    for (dup_name in duplicated_names) {
      dup_indices <- which(vertex_names == dup_name)
      n <- vapply(dup_indices, function(idx) length(incident(g, idx, mode = "all")), integer(1))
      x <- which(
        correct_n_edges$graph_id == g_name &
          toupper(correct_n_edges$names) %in% toupper(dup_name)
      )
      if (!length(x)) {
        warning("[11.1.1] no correct_n_edges row for ", dup_name, " in ", g_name)
        next
      }
      expected <- as.numeric(correct_n_edges$n_edge[x])
      observed <- as.numeric(n)
      if (length(expected) > length(observed)) expected <- expected[seq_along(observed)]
      keep_index <- which.min(abs(observed - expected))
      correct_n_edges$vertex_index_to_remove[x] <- list(dup_indices[-keep_index])
    }
  }

  zero_indices <- vapply(
    correct_n_edges$vertex_index_to_remove,
    function(x) length(x) == 1L && x == 0L,
    logical(1)
  )
  if (any(zero_indices)) {
    correct_n_edges$vertex_index_to_remove[zero_indices] <- NA
  }

  correct_n_edges_out <- correct_n_edges
  correct_n_edges_out$vertex_index_to_remove <- vapply(
    correct_n_edges_out$vertex_index_to_remove,
    function(x) {
      if (is.null(x) || all(is.na(x))) return(NA_character_)
      paste(x, collapse = ",")
    },
    character(1)
  )

  write.table(
    correct_n_edges_out,
    file = sub("\\.rds$", ".txt", basename(correct_edges_file)),
    sep = "\t", row.names = FALSE
  )
  saveRDS(correct_n_edges_out, file = correct_edges_file)
  message("[11.1.1] wrote ", correct_edges_file, " (", nrow(correct_n_edges_out), " rows)")

  if (!is.null(dup_review) && nrow(dup_review)) {
    write.table(
      dup_review,
      file = correct_edges_review_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    message("[11.1.1] wrote auto-pick review -> ", correct_edges_review_file)

    if (!use_mapping_cache) {
      dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
      write.table(
        dup_review,
        file = string_mapping_correction_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      message("[11.1.1] wrote shared correction cache -> ", string_mapping_correction_file)
    }
  }
}
# --- manual run log (backup; CP4 arms) ---
# [11.1.1] duplicate vertices in: HiG_0, HiG_1, HiG_10, HiG_11, HiG_13, HiG_2, HiG_3, HiG_5, HiG_7, HiG_8, HiG_9, HiG_12, HiG_6, HiG_4
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_0 duplicated symbols: NOV, HIST1H2BC
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_1 duplicated symbols: NOV, HIST1H2BC
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_10 duplicated symbols: AKAP2, HIST1H2BC
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_11 duplicated symbols: HIST1H2BC
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_13 duplicated symbols: TRAC, HIST1H2BC
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_2 duplicated symbols: NOV, HIST1H2BC
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_3 duplicated symbols: NOV
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_5 duplicated symbols: NOV, HIST1H2BC
# Warning:  we couldn't map to STRING 2% of your identifiers[11.1.1] HiG_7 duplicated symbols: NOV
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_8 duplicated symbols: NOV
# Warning:  we couldn't map to STRING 1% of your identifiers[11.1.1] HiG_9 duplicated symbols: NOV, AKAP2
# Warning:  we couldn't map to STRING 1% of your identifiers[11.1.1] HiG_12 duplicated symbols: AKAP2
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_6 duplicated symbols: AKAP2
# Warning:  we couldn't map to STRING 3% of your identifiers[11.1.1] HiG_4 duplicated symbols: AKAP2, HIST1H2BC
# [11.1.1] wrote F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4//correct_n_edges_HiG_STRING2.14.0.rds (20 rows)
# [11.1.1] wrote auto-pick review -> F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4//correct_n_edges_HiG_STRING2.14.0_auto_pick_review.txt
# [11.1.1] wrote shared correction cache -> data_dir/unique_STRING_mapping_correction.txt (after first successful biomaRt run)
