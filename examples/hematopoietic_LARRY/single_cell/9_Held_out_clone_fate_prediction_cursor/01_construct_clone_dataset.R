## 01_construct_clone_dataset.R
## Step 1: clone-level dataset + feasibility counts (no predictive performance).
##
## For every clone with C11 cells at day 2:
##   n_d2_c11, n_d6, n_meg_d6, n_baso_d6, well-specific day-6 counts.
## Fate counts use the full public metadata (130,887 rows). Expression/state
## use the working NM-trajectory object (92,527).

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_configure()

if (!file.exists(cfg$public_meta)) stop("Missing public metadata: ", cfg$public_meta)
if (!file.exists(cfg$public_clone)) stop("Missing clone matrix: ", cfg$public_clone)

message("[01] reading public metadata")
pub <- read_public_meta(cfg$public_meta)
message("[01] public n=", nrow(pub), " libraries=", length(unique(pub$Library)))

message("[01] reading clone matrix")
clone_mat <- read_clone_matrix(cfg$public_clone)
if (nrow(clone_mat) != nrow(pub)) {
  stop("Clone matrix rows (", nrow(clone_mat), ") != public metadata rows (", nrow(pub), ")")
}
message("[01] clone matrix ", nrow(clone_mat), " x ", ncol(clone_mat))

message("[01] reading working-object metadata")
work <- read_working_meta(cfg)
work$clone_id <- assign_clone_ids(clone_mat, work$cell_index0)
work$n_clones <- attr(work$clone_id, "n_clones")
work$multi_clone <- attr(work$clone_id, "multi_clone")
n_multi <- attr(work$clone_id, "n_multi_clone_cells")
message("[01] working cells with >1 clone: ", n_multi, " (excluded from clone_id, not assigned to the first clone)")

pub$clone_id <- assign_clone_ids(clone_mat, pub$cell_index0)
n_multi_pub <- attr(pub$clone_id, "n_multi_clone_cells")
message("[01] public cells with >1 clone: ", n_multi_pub, " (excluded from fate counts)")

state <- as.character(work[[cfg$state_col]])
is_c11 <- state == as.character(cfg$state_id)
is_d2 <- work$Time.point == cfg$day_early
d2_c11 <- work[is_c11 & is_d2, , drop = FALSE]

n_d2_c11 <- nrow(d2_c11)
n_d2_c11_multi <- sum(d2_c11$multi_clone, na.rm = TRUE)
n_d2_c11_cloned <- sum(!is.na(d2_c11$clone_id))
c11_clones <- sort(unique(d2_c11$clone_id[!is.na(d2_c11$clone_id)]))

message("[01] C11 day-2 cells: ", n_d2_c11,
        "; with clone assignment: ", n_d2_c11_cloned,
        "; distinct C11 clones: ", length(c11_clones))

## clone-level rows for C11 clones
agg_cells <- function(df, clone_ids) {
  df <- df[!is.na(df$clone_id) & df$clone_id %in% clone_ids, , drop = FALSE]
  if (!nrow(df)) {
    return(data.frame(clone_id = integer(), n = integer(), stringsAsFactors = FALSE))
  }
  as.data.frame(table(df$clone_id), stringsAsFactors = FALSE)
}

d2_c11_n <- as.data.frame(table(d2_c11$clone_id[!is.na(d2_c11$clone_id)]), stringsAsFactors = FALSE)
names(d2_c11_n) <- c("clone_id", "n_d2_c11")
d2_c11_n$clone_id <- as.integer(d2_c11_n$clone_id)

d2_all <- work[work$Time.point == cfg$day_early & !is.na(work$clone_id) & work$clone_id %in% c11_clones, ]
d2_all_n <- as.data.frame(table(d2_all$clone_id), stringsAsFactors = FALSE)
names(d2_all_n) <- c("clone_id", "n_d2_any")
d2_all_n$clone_id <- as.integer(d2_all_n$clone_id)

late <- pub[!is.na(pub$clone_id) & pub$clone_id %in% c11_clones & pub$Time.point == cfg$day_late, ]
late$is_meg <- as.character(late$Cell.type.annotation) == cfg$meg_label
late$is_baso <- as.character(late$Cell.type.annotation) == cfg$baso_label
late$is_undiff <- is_undiff_label(late$Cell.type.annotation)
late$is_mature <- !late$is_undiff

summarise_late <- function(df) {
  if (!nrow(df)) {
    return(data.frame(
      clone_id = integer(), n_d6 = integer(), n_meg_d6 = integer(), n_baso_d6 = integer(),
      n_undiff_d6 = integer(), n_mature_d6 = integer(),
      stringsAsFactors = FALSE
    ))
  }
  stats::aggregate(
    cbind(
      n_d6 = 1L, n_meg_d6 = as.integer(df$is_meg), n_baso_d6 = as.integer(df$is_baso),
      n_undiff_d6 = as.integer(df$is_undiff), n_mature_d6 = as.integer(df$is_mature)
    ),
    by = list(clone_id = df$clone_id),
    FUN = sum
  )
}

late_all <- summarise_late(late)
late_w1 <- summarise_late(late[late$Well == 1L, ])
names(late_w1)[-1] <- paste0(names(late_w1)[-1], "_well1")
late_w2 <- summarise_late(late[late$Well == 2L, ])
names(late_w2)[-1] <- paste0(names(late_w2)[-1], "_well2")

## clone covariates from day-2 C11 cells (mode)
mode_chr <- function(x) {
  x <- as_chr(x)
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

covar <- do.call(rbind, lapply(split(d2_c11[!is.na(d2_c11$clone_id), ], d2_c11$clone_id[!is.na(d2_c11$clone_id)]), function(dd) {
  data.frame(
    clone_id = dd$clone_id[1],
    library = mode_chr(dd$Library),
    starting_population = mode_chr(dd$Starting.population),
    n_libraries_d2 = length(unique(dd$Library)),
    stringsAsFactors = FALSE
  )
}))

clones <- data.frame(clone_id = c11_clones, stringsAsFactors = FALSE)
clones <- merge(clones, d2_c11_n, by = "clone_id", all.x = TRUE)
clones <- merge(clones, d2_all_n, by = "clone_id", all.x = TRUE)
clones <- merge(clones, late_all, by = "clone_id", all.x = TRUE)
clones <- merge(clones, late_w1, by = "clone_id", all.x = TRUE)
clones <- merge(clones, late_w2, by = "clone_id", all.x = TRUE)
clones <- merge(clones, covar, by = "clone_id", all.x = TRUE)

num_cols <- c(
  "n_d2_c11", "n_d2_any", "n_d6", "n_meg_d6", "n_baso_d6",
  "n_undiff_d6", "n_mature_d6",
  "n_d6_well1", "n_meg_d6_well1", "n_baso_d6_well1",
  "n_undiff_d6_well1", "n_mature_d6_well1",
  "n_d6_well2", "n_meg_d6_well2", "n_baso_d6_well2",
  "n_undiff_d6_well2", "n_mature_d6_well2"
)
for (cc in intersect(num_cols, names(clones))) {
  clones[[cc]][is.na(clones[[cc]])] <- 0L
  clones[[cc]] <- as.integer(clones[[cc]])
}
clones$frac_meg_d6 <- ifelse(clones$n_d6 > 0, clones$n_meg_d6 / clones$n_d6, NA_real_)
clones$frac_baso_d6 <- ifelse(clones$n_d6 > 0, clones$n_baso_d6 / clones$n_d6, NA_real_)
clones$frac_meg_mature_d6 <- ifelse(clones$n_mature_d6 > 0, clones$n_meg_d6 / clones$n_mature_d6, NA_real_)
clones$frac_baso_mature_d6 <- ifelse(clones$n_mature_d6 > 0, clones$n_baso_d6 / clones$n_mature_d6, NA_real_)
clones$frac_meg_d6_well1 <- ifelse(clones$n_d6_well1 > 0, clones$n_meg_d6_well1 / clones$n_d6_well1, NA_real_)
clones$frac_meg_d6_well2 <- ifelse(clones$n_d6_well2 > 0, clones$n_meg_d6_well2 / clones$n_d6_well2, NA_real_)
clones$in_both_later_wells <- clones$n_d6_well1 > 0 & clones$n_d6_well2 > 0
clones$has_d6 <- clones$n_d6 > 0
clones$has_meg <- clones$n_meg_d6 > 0
clones$has_baso <- clones$n_baso_d6 > 0

cell_map <- d2_c11[!is.na(d2_c11$clone_id), c("cell_id", "clone_id", "Library", "Cell.barcode",
                                             "Time.point", "Well", "Starting.population",
                                             "Cell.type.annotation", cfg$state_col, "cell_index0")]
names(cell_map)[names(cell_map) == cfg$state_col] <- "state"

saveRDS(list(
  clones = clones,
  d2_c11_cells = cell_map,
  work_meta = work[, c("cell_id", "clone_id", "cell_index0", "Time.point", "Well",
                       cfg$state_col, "Cell.type.annotation", "Library",
                       "n_clones", "multi_clone")]
), file.path(cfg$results_dir, "rds", "01_clone_dataset.rds"))
write_tsv(clones, file.path(cfg$results_dir, "tables", "01_clone_level_dataset.tsv"))
write_tsv(cell_map, file.path(cfg$results_dir, "tables", "01_d2_c11_cells.tsv"))

counts <- data.frame(
  item = c(
    "working_object_cells",
    "working_multi_clone_cells_dropped",
    "public_cells",
    "public_multi_clone_cells_dropped",
    "public_clones",
    "C11_day2_cells",
    "C11_day2_cells_multi_clone_dropped",
    "C11_day2_cells_with_unique_clone",
    "distinct_C11_clones",
    "C11_clones_with_annotated_day6_progeny",
    "C11_clones_with_at_least_one_Meg",
    "C11_clones_with_at_least_one_Baso",
    "C11_clones_in_both_later_wells"
  ),
  n = c(
    nrow(work),
    n_multi,
    nrow(pub),
    n_multi_pub,
    ncol(clone_mat),
    n_d2_c11,
    n_d2_c11_multi,
    n_d2_c11_cloned,
    length(c11_clones),
    sum(clones$has_d6),
    sum(clones$has_meg),
    sum(clones$has_baso),
    sum(clones$in_both_later_wells)
  ),
  stringsAsFactors = FALSE
)
write_tsv(counts, file.path(cfg$results_dir, "01_feasibility_counts.tsv"))
print(counts)

qs <- function(x) {
  x <- x[is.finite(x)]
  qs <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE)
  data.frame(min = qs[[1]], q25 = qs[[2]], median = qs[[3]], q75 = qs[[4]], q90 = qs[[5]], max = qs[[6]])
}
dist <- rbind(
  cbind(quantity = "n_d2_c11", qs(clones$n_d2_c11)),
  cbind(quantity = "n_d6", qs(clones$n_d6)),
  cbind(quantity = "n_d6_among_clones_with_d6", qs(clones$n_d6[clones$has_d6])),
  cbind(quantity = "n_mature_d6_among_clones_with_d6", qs(clones$n_mature_d6[clones$has_d6])),
  cbind(quantity = "n_meg_d6_among_meg_pos", qs(clones$n_meg_d6[clones$has_meg])),
  cbind(quantity = "n_d6_well1", qs(clones$n_d6_well1)),
  cbind(quantity = "n_d6_well2", qs(clones$n_d6_well2))
)
write_tsv(dist, file.path(cfg$results_dir, "tables", "01_size_distributions.tsv"))
print(dist)

grDevices::pdf(file.path(cfg$results_dir, "figures", "01_clone_size_histograms.pdf"), width = 10, height = 6)
old <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
hist(clones$n_d2_c11, breaks = 40, main = "Day-2 C11 cells / clone", xlab = "n")
hist(clones$n_d6, breaks = 40, main = "Day-6 progeny / clone", xlab = "n")
hist(clones$n_d6[clones$has_d6], breaks = 40, main = "Day-6 n | n_d6>0", xlab = "n")
hist(clones$n_meg_d6, breaks = 40, main = "Day-6 Meg / clone", xlab = "n")
hist(clones$n_d6_well1, breaks = 40, main = "Day-6 well 1 / clone", xlab = "n")
hist(clones$n_d6_well2, breaks = 40, main = "Day-6 well 2 / clone", xlab = "n")
par(old)
grDevices::dev.off()

message("[01] done. Inspect 01_feasibility_counts.tsv and size histograms before any prediction.")
message("[01] Choose min clone-size criteria from these distributions only — not from Meg-prediction performance.")

note <- c(
  "FEASIBILITY (written by 01, before any scores):",
  paste0("C11 day-2 cells=", n_d2_c11, "; with clone=", n_d2_c11_cloned, "; distinct C11 clones=", length(c11_clones), "."),
  paste0("With day-6 progeny=", sum(clones$has_d6), "; Meg-positive=", sum(clones$has_meg), "; Baso-positive=", sum(clones$has_baso), "; both later wells=", sum(clones$in_both_later_wells), "."),
  "Meg-positive C11 clones are uncommon. Keep n_d6>=1 and let quasibinomial weight by n_d6;",
  "min_n_d6=5 would drop most Meg-positive clones. Two-well tests will be underpowered.",
  "Baso-positive clones are more common than Meg-positive, so Baso held-out failure cannot be blamed on sample size alone."
)
writeLines(note, file.path(cfg$results_dir, "01_feasibility_NOTE.txt"))
message(paste(note, collapse = "\n"))
