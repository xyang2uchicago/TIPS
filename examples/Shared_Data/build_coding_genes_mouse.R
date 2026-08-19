## Build mouse protein-coding gene list and save to Shared_Data/coding_genes_mouse.rds
## Output: character vector of MGI gene symbols for protein-coding genes
## Run once; downstream scripts read coding_genes_mouse.rds directly.

shared_path <- here::here("examples", "Shared_Data/")

########################################################
## Method 1: org.Mm.eg.db (offline, primary)
########################################################
library(org.Mm.eg.db)
library(AnnotationDbi)

cols <- columns(org.Mm.eg.db)
cat("Available columns in org.Mm.eg.db:\n")
print(cols)

if ("GENETYPE" %in% cols) {
    all_entries <- select(org.Mm.eg.db,
        keys    = keys(org.Mm.eg.db, keytype = "ENTREZID"),
        columns = c("SYMBOL", "GENETYPE"),
        keytype = "ENTREZID"
    )
    coding_genes_mouse <- unique(all_entries$SYMBOL[
        !is.na(all_entries$GENETYPE) & all_entries$GENETYPE == "protein-coding"
    ])
    coding_genes_mouse <- coding_genes_mouse[!is.na(coding_genes_mouse) & coding_genes_mouse != ""]
    cat("Method 1 (org.Mm.eg.db GENETYPE): ", length(coding_genes_mouse), "protein-coding genes\n")
    # Method 1 (org.Mm.eg.db GENETYPE):  26310 protein-coding genes
} else {
    message("GENETYPE column not available in this version of org.Mm.eg.db — falling back to Method 2")
    coding_genes_mouse <- NULL
}

########################################################
## Method 2: biomaRt (online fallback, or if Method 1 result seems too small)
########################################################
if (is.null(coding_genes_mouse) || length(coding_genes_mouse) < 15000) {
    message("Trying biomaRt for mouse protein-coding genes...")
    library(biomaRt)
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    res <- getBM(
        attributes = c("mgi_symbol", "gene_biotype"),
        mart       = mart
    )
    coding_mart <- unique(res$mgi_symbol[res$gene_biotype == "protein_coding"])
    coding_mart <- coding_mart[!is.na(coding_mart) & coding_mart != ""]
    cat("Method 2 (biomaRt): ", length(coding_mart), "protein-coding genes\n")

    if (is.null(coding_genes_mouse)) {
        coding_genes_mouse <- coding_mart
    } else {
        ## take union of both sources
        coding_genes_mouse <- unique(c(coding_genes_mouse, coding_mart))
        cat("Union of Method 1 + Method 2: ", length(coding_genes_mouse), "genes\n")
    }
}

########################################################
## Summary and save
########################################################
cat("\nFinal mouse coding gene list:\n")
cat("  Total:", length(coding_genes_mouse), "\n")
# Total: 26310
cat("  First 10:", head(coding_genes_mouse, 10), "\n")
# First 10: Pzp2 Aanat Aatk Abca1 Abca4 Abca2 Abcb7 Abcg1 Abi1 Abl1

## cross-check: overlap with human coding_genes.rds via toupper
coding_human <- unique(readRDS(paste0(shared_path, "coding_genes.rds")))
overlap <- length(intersect(toupper(coding_genes_mouse), coding_human))
cat("  Overlap with human list (via toupper):", overlap,
    "(", round(100 * overlap / length(coding_genes_mouse), 1), "% of mouse)\n")
# Overlap with human list (via toupper): 16400 ( 62.3 % of mouse)

saveRDS(coding_genes_mouse, file = paste0(shared_path, "coding_genes_mouse.rds"))
cat("\nSaved to:", paste0(shared_path, "coding_genes_mouse.rds"), "\n")
