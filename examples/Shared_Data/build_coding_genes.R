## Build human protein-coding gene list and save to Shared_Data/coding_genes.rds
## Output: character vector of HGNC gene symbols for protein-coding genes
## Run once; downstream scripts read coding_genes.rds directly.

shared_path <- here::here("examples", "Shared_Data/")

########################################################
## Method 1: org.Hs.eg.db (offline, primary)
########################################################
library(org.Hs.eg.db)
library(AnnotationDbi)

cols <- columns(org.Hs.eg.db)
cat("Available columns in org.Hs.eg.db:\n")
print(cols)

if ("GENETYPE" %in% cols) {
    all_entries <- select(org.Hs.eg.db,
        keys    = keys(org.Hs.eg.db, keytype = "ENTREZID"),
        columns = c("SYMBOL", "GENETYPE"),
        keytype = "ENTREZID"
    )
    coding_genes <- unique(all_entries$SYMBOL[
        !is.na(all_entries$GENETYPE) & all_entries$GENETYPE == "protein-coding"
    ])
    coding_genes <- coding_genes[!is.na(coding_genes) & coding_genes != ""]
    cat("Method 1 (org.Hs.eg.db GENETYPE): ", length(coding_genes), "protein-coding genes\n")
} else {
    message("GENETYPE column not available in this version of org.Hs.eg.db — falling back to Method 2")
    coding_genes <- NULL
}

########################################################
## Method 2: biomaRt (online fallback, or if Method 1 result seems too small)
########################################################
if (is.null(coding_genes) || length(coding_genes) < 15000) {
    message("Trying biomaRt for human protein-coding genes...")
    library(biomaRt)
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    res <- getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        mart       = mart
    )
    coding_mart <- unique(res$hgnc_symbol[res$gene_biotype == "protein_coding"])
    coding_mart <- coding_mart[!is.na(coding_mart) & coding_mart != ""]
    cat("Method 2 (biomaRt): ", length(coding_mart), "protein-coding genes\n")

    if (is.null(coding_genes)) {
        coding_genes <- coding_mart
    } else {
        ## take union of both sources
        coding_genes <- unique(c(coding_genes, coding_mart))
        cat("Union of Method 1 + Method 2: ", length(coding_genes), "genes\n")
    }
}

########################################################
## Summary and save
########################################################
cat("\nFinal human coding gene list:\n")
cat("  Total:", length(coding_genes), "\n")
cat("  First 10:", head(coding_genes, 10), "\n")

## cross-check: overlap with mouse coding_genes_mouse.rds via toupper
coding_mouse_rds <- paste0(shared_path, "coding_genes_mouse.rds")
if (file.exists(coding_mouse_rds)) {
    coding_mouse <- unique(readRDS(coding_mouse_rds))
    overlap <- length(intersect(coding_genes, toupper(coding_mouse)))
    cat("  Overlap with mouse list (via toupper):", overlap,
        "(", round(100 * overlap / length(coding_genes), 1), "% of human)\n")
}

saveRDS(coding_genes, file = paste0(shared_path, "coding_genes.rds"))
cat("\nSaved to:", paste0(shared_path, "coding_genes.rds"), "\n")
