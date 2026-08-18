# Held-out clone-fate prediction (LARRY / GSE140802)

This folder implements the design in `F:/projects/TIPS/doc/NatureComm/Larry_prediction.docx`.

Call this **held-out prediction of later clonal fate**, not prospective prediction. Day-2 sequenced cells were destroyed during profiling; they are barcode-matched clonal relatives of the day-6 cells, not literal ancestors.

```r
source("F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor/run_heldout_pipeline.R")
```

Default: run **Step 1 only** (clone dataset + feasibility counts), then stop. Inspect `01_feasibility_counts.tsv` and size histograms **before** looking at predictive performance. Then either:

```r
Sys.setenv(HELDOUT_RUN = "all")   # auto-lock size criteria from size-only quantiles, then predict
source(".../run_heldout_pipeline.R")
```

or edit `02_clone_size_criteria_PROPOSED.tsv`, copy it to `02_clone_size_criteria_LOCKED.tsv`, and run `HELDOUT_RUN=from_lock`.

Results: `F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/9_heldout_clone_fate/`

## Script map

| Script | Design step |
|--------|-------------|
| `01_construct_clone_dataset.R` | Clone-level dataset; counts; distributions |
| `02_lock_criteria_and_folds.R` | Size criteria (before performance); clone-grouped CV |
| `03_gene_sets.R` | Frozen comparator gene lists |
| `03b_build_tips_module_fold.R` | Per-fold TIPS: sources real `code_core_11_10vs17` 11.1–24.3 on training clones |
| `04_heldout_prediction.R` | Score held-out day-2 C11 cells; pool OOF |
| `05_stats_baselines_wells.R` | Quasibinomial OR, Spearman, bootstrap, permutation, AUROC/AUPRC, two-well |
| `06_branch_resolvability.R` | 8-row Meg vs Baso branch-resolvability table + Baso gene audit |
| `07_figures.R` | Main-figure panels a–g |
| `08_summary_figure.R` | Design three-panel + Meg held-out summary |
| `10_clone_bootstrap.R` | Paired clone-bootstrap CIs for TIPS − CTS / PPI / training Meg markers / node DE / MuTrans |
| `11_same_size_null.R` | Same-size CTS subsets and degree-matched PPI subnetworks vs fold TIPS modules |

## Branch-resolvability table (06)

Eight matched Meg vs Baso rows, computed before looking at clone-fate outcomes. Held-out prediction is joined last (row 8) and is not used to change Baso thresholds.

| Row | Diagnostic | How it is computed |
|-----|------------|--------------------|
| 1 | BioTIP empirical significance | Real C11 CTS summary: `MCI_P`, `IC_g_local_p`, `IC_s_local_p`, `localmax` (same state for both arms). |
| 2 | CTS–MuTrans overlap | Fisher's exact test of C11 CTS vs MuTrans TD gene lists; Holly `FET_venn3_*.tsv` appended when present. |
| 3 | Backbone coverage | Fraction of CTS ∪ TIPS-arm genes that are nodes **and** have degree ≥ 1 in the C11 STRING graph (`GSE140802_STRING_graph_perState_notsimplified.rds`, else the weighted HiG file). |
| 4 | Edge evidence | Observed mean delta from `edge_change_table(HiG_11, HiG_17/10)` restricted to dual-pull TIPS edges (dual-pull table if those edges are missing from the graph). **p-value:** Pearson Δ-correlation on those same TIPS edges after shuffling descendant labels among C11+arm cells. That null is not a full TIPS reweight per permutation. |
| 5 | Network integrity | Loads the **weighted** retained HiG (`PPI_weight/...combinedweighted.rds`); falls back to notsimplified. Reports `vcount`, `ecount`, components, and giant-component size for Meg `HiG_17` and Baso `HiG_10`. |
| 6 | Resampling stability | Bootstrap cells in C11 vs CF17/CM10; rank STRING pairs (else all pairs) by Pearson delta; mean gene/edge Jaccard across bootstraps and TIPS-gene selection frequency. |
| 7 | STRING–IID agreement | No IID arm at `leiden_r0_8`. Uses `jaccard_STRING_vs_IID_4_9vs11.tsv` and, if present, a metacell `code_core_4_9vs11` STRING vs IID vertex Jaccard. |
| 8 | External validation | Weinreb Table S3 FET on the frozen TIPS gene set, then 05 held-out rho / OR / AUROC joined last. |

Panel f in `07_figures.R` reads this 8-row table (`level`, `Meg`, `Baso`).


## Matching rule

Do **not** join on cell barcode alone. Public LARRY rows (130,887) are indexed by `cell_index0`. The working NM-trajectory object (92,527 cells) already stores `Library`, `Cell.barcode`, `Time.point`, `Well`, and `cell_index0`. Clone IDs are taken from `stateFate_inVitro_clone_matrix.mtx.gz` via that index.

Day-2 C11 **expression and state** come from the working object. Day-6 **fate counts** come from the **full public metadata**, so Meg/Baso progeny that left the NM-trajectory subset are still counted.

## Disagreements / locked choices

These are places I would not follow a naive reading of the design, or where the design left a gap:

1. **C11 is cell-level `leiden_r0_8` cluster 11, not MuTrans attractor 11.** Attractors 4 / 9 / 11 are metacell progenitor / Baso / Meg. Scoring day-2 clones on attractor 11 would score Meg cells. Primary TIPS module is the cell-level run `11_10vs17` (CP11 → CF17 Meg / CM10 Baso). Metacell TIPS `4_9vs11` is concordance only.

2. **Later fate uses Weinreb `Cell.type.annotation` (`Meg`, `Baso`), not leiden 17 / 10.** Cluster 17 is the TIPS CF construction cluster; using it as the endpoint would leak the clustering into the outcome.

3. **Default TIPS construction re-invokes `code_core_11_10vs17` (11.1–24.3) on training-clone cells.** C11/CTS stay frozen (design allows that). DEG, STRING reweighting, PageRank, and cisTarget dual-pull are refit per fold. `HELDOUT_TIPS_MODE=proxy` is a Pearson-delta smoke test, not the reported algorithm. `frozen_genes` is a leakage diagnostic.

4. **TIPS cell score = mean of training-fit gene-wise z-scored log-normalized expression.** The design did not specify a scoring rule. I will not invent PageRank-weighted scores as primary. Median is a sensitivity (`HELDOUT_SCORE_FUN=median`).

5. **AUROC/AUPRC stay secondary.** A Meg-positive clone is locked as `n_meg_d6 >= 1` among size-filtered clones, before any score is computed. Do not tune that cutoff on performance.

6. **Weinreb Table S3 fate genes are not a held-out comparator** (agreed). They appear only in the resolvability/enrichment audit.

7. **Do not change Baso TIPS thresholds after seeing clone-fate results** (agreed). `06` writes internal diagnostics first; held-out columns are joined last and are not inputs to any filter.

8. **With two arms this is branch-level resolvability, not a statistically calibrated confidence score** (agreed). Extending to all C11-to-mature arms is optional (`HELDOUT_ALL_ARMS=1`) and only after Meg/Baso works.

9. **Clone-size auto-lock uses n_d6 >= 1**, not a high size cutoff. Step 1 found only 25 Meg-positive C11 clones (8 if n_d6>=5). The design asked to keep numerator and denominator because clones differ in recovered progeny; quasibinomial already weights by n_d6. Dropping small clones would discard most of the Meg signal. Sensitivity still reports 3 / 5 / 10.

10. **25 Meg-positive clones is thin.** For `proxy` / `frozen_genes`, 02 auto-switches to repeated grouped 5-fold CV. Default `run_tips_pipeline` stays at 1 repeat (5 full TIPS runs). Two-well Meg tests will be underpowered. Baso has 61 positive clones, so a Baso held-out failure is not a sample-size excuse.

11. **Rows 4 (p-value) and 6 are Pearson Δ-correlation proxies** on existing TIPS/STRING edges. Observed row 4 also uses `edge_change_table()` on the real reweighted graphs. Neither permutes the full STRING-reweight + PageRank + cisTarget pipeline.

12. **Multi-clone cells are dropped**, not assigned to the first clone. Fold TIPS objects must have zero remaining test-clone cells by `cell_index0`. Empty TIPS modules keep every eligible test clone (score NA) and are counted as method failure, not missing data.

## What Step 1 found (already run)

| Item | n |
|------|---|
| C11 day-2 cells | 2331 |
| with clone assignment | 448 |
| distinct C11 clones | 318 |
| with day-6 progeny | 155 |
| Meg-positive | 25 |
| Baso-positive | 61 |
| in both later wells | 33 |

## Primary statistical test

Primary quasibinomial: `cbind(n_meg_d6, n_d6 - n_meg_d6) ~ scale(day2_TIPS) + library + starting_population`. Score-only OR is reported as a sensitivity. Denominator `n_d6` is **all recovered day-6 progeny** (including undifferentiated). Meg / mature progeny is a sensitivity (`n_mature_d6`). Secondary: out-of-fold Spearman, clone bootstrap CI, within-library fate permutation, AUROC/AUPRC. In-sample log-loss is not reported.

## Clone-bootstrap Δρ (10)

Paired clone bootstrap of **TIPS minus baseline** on the pooled OOF scores (gene lists are not refit). Same 155 clones as step 05. Default 5000 resamples (`HELDOUT_N_BOOT_DELTA`). Focal comparators: fixed C11 CTS, unweighted PPI, training Meg markers, node-level DE, MuTrans TD.

Primary statistic: Spearman vs later Meg fraction. Secondary: AUROC. Sensitivity: quasibinomial OR (`HELDOUT_BOOT_OR=0` to skip). Reports percentile and BCa 95% CIs, a Meg-positive–stratified CI, and a Williams/Steiger test on ranks. Inference is the **two-sided 95% CI**. A one-sided superiority test was not prespecified; do not treat one-sided bootstrap P(Δ ≤ 0) as a significance claim.

If the TIPS−CTS CI for Δρ includes 0, **do not claim predictive superiority**. The selling point remains branch-specific regulatory resolution and interpretability.

## Same-size network null (11)

Step 10 compares TIPS (k ≈ 5) to the **full** 71-gene CTS, which mixes module size with gene choice. Step 11 holds size fixed and asks whether the fold TIPS module is an unusually fate-informative compact set:

1. Random k-subsets of C11 CTS
2. Degree-matched k nodes from the C11 STRING graph
3. Connected k-node STRING subgraphs, then the candidate whose sorted log-degree vector is closest to TIPS

Each draw is scored with the same out-of-fold rule as TIPS. TIPS is **not** a CTS subset (`GATA2` is outside CTS); the CTS null is a competing compact CTS module. The existing HVG `random_matched` null in 04/05 is a different (and currently degenerate) test.

```r
Sys.setenv(HELDOUT_RUN = "11")
source("F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor/run_heldout_pipeline.R")
```

```r
Sys.setenv(HELDOUT_RUN = "10")
source("F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor/run_heldout_pipeline.R")
```

## Baselines (same held-out endpoint)

- TIPS Meg module (per-fold `code_core_11_10vs17` 11.1–24.3 on training clones)
- C11 CTS genes
- MuTrans Meg transition-driver genes (`td_genes_scores_4_to_11`, `|corr| > 0.7`)
- Training-derived mature Meg markers (day-6 Meg vs other, training clones only)
- Node-only progenitor-to-Meg DE (training C11 vs Meg)
- Unweighted PPI module (dual-pull nodes, no positive-delta filter)
- Expression- and size-matched random gene sets
- Frozen fate-blind TIPS genes (leakage diagnostic only)

## Environment

| Variable | Default | Purpose |
|----------|---------|---------|
| `HELDOUT_RUN` | `step1` | `step1`, `all`, `from_lock`, or a script number |
| `HELDOUT_TIPS_MODE` | `run_tips_pipeline` | `proxy` (Pearson smoke test) or `frozen_genes` |
| `HELDOUT_OVERWRITE_FOLDS` | unset | `HELDOUT_RUN=all` defaults to `1`. Set `0` to reuse identity-matched fold caches. |
| `HELDOUT_N_REPEATS` | `1` | extra CV repeats (each is a full 11.1–24.3 run per fold) |
| `HELDOUT_N_FOLDS` | `5` | clone-grouped folds |
| `HELDOUT_SEED` | `140802` | RNG |
| `HELDOUT_EDGE_N_PERM` | `200` | Row 4 descendant-label shuffles |
| `HELDOUT_STAB_N_BOOT` | `50` | Row 6 cell bootstraps |
| `HELDOUT_AUTOLOCK` | `0` | `1` to write LOCKED from PROPOSED without copying |
| `HELDOUT_N_BOOT_DELTA` | `5000` | clone resamples in step 10 |
| `HELDOUT_BOOT_OR` | `1` | `0` to skip quasibinomial OR in step 10 |
| `HELDOUT_N_NULL` | `500` | same-size null draws per family in step 11 |
