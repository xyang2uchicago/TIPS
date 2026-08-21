# TIPS

TIPS links transition-prone states in single-cell data to interpretable interaction subnetworks, prioritizing bridge regulators and “dual-pull” lineage programs during fate transitions. The rationale is critical-transition genes-derived progenitor networks were ordered yet selectively fragile. For more details, refer to [**regulated stochasticity tutorial**](https://xyang2uchicago.github.io/TIPS/tutorial/TIPS.html)

[**DOI**](https://doi.org/10.5281/zenodo.22031075)

### Repository Structure

```
TIPS/
├── R/                              # shared helper functions sourced by analysis scripts
├── imgs/                           # figures used in this README
├── tutorial/                       # step-by-step walkthrough (GSE87038)
└── examples/
    ├── config.R                    # here::here()-based path resolution, shared by all trees
    ├── Shared_Data/                 # cross-tree reference data (STRING/IID PPINs, cisTarget, gene lists)
    │
    ├── cardiac/                     # case study: cardiac development
    │   ├── data/                    # cardiac-only shared reference files
    │   ├── evaluation/              # cross-dataset evaluation figures/scripts
    │   ├── GSE175634/
    │   │   ├── data/
    │   │   ├── GSE175634_IID/
    │   │   │   ├── code/            # original analysis scripts
    │   │   │   ├── code_core/       # cleaned, portable pipeline
    │   │   │   └── results_core/    # networks, figures, tables
    │   │   └── GSE175634_STRING/    # same code / code_core / results_core layout
    │   ├── GSE87038/                 # same IID / STRING × code_core / results_core layout
    │   └── IbarraSoria2018/          # same layout
    │
    ├── hematoendothelial/           # case study: hematoendothelial development
    │   ├── GSE87038/                 # code_core_13 + results_core_13
    │   └── IbarraSoria2018/          # code_core_endothelial.b + results_core_endothelial.b
    │
    └── hematopoietic_LARRY/         # case study: hematopoiesis
        ├── data/                     # LARRY-wide shared reference data
        ├── metacell/
        │   ├── 7_data_MuTrans_TIPS_IID/
        │   │   ├── code_core_4_9vs11/
        │   │   └── results_core_4_9vs11/
        │   └── 7_data_MuTrans_TIPS_STRING/   # same layout
        └── single_cell/
            ├── code_core_11_10vs17/
            ├── results_core_11_10vs17/
            └── 9_Held_out_clone_fate_prediction_cursor/  # held-out clone-fate validation
```

### What is TIPS?
- Name: **Transcriptional Instability–guided Prediction of Subnetworks**.
- Goal: Identify **semi-stable transition states** obscured by program mixing (Fig. 1A).
- Output: Quantify **state-specific network rewiring** from uncommitted progenitors to descendant fates (Fig 1B).
- Method: BioTIP + **coexpression-weighted protein–protein interaction networks + robustness (targeted-attack) + chromatin integration** to resolve lineage-leaning arms (Fig. 1C).
<br>
<p align="center">
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_developmental_blackbox.jpg" width="33%" height="auto">
<br>
Figure 1A
<br>
<br>
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_pipeline.jpg" width="33%" height="auto">
<br>
Figure 1B
<br>
<br>
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_analysis_strategy.jpg" width="33%" height="auto">
<br>
Figure 1C
</p>
<br>

### Why TIPS?
- Vulnerable progenitor states are transient, heterogeneous, and difficult to capture using markers alone.
- Resolving how perturbations rewire developmental programs remains challenging.
- TIPS maps defect-vulnerable windows and prioritizes lineage-selective regulatory programs from mixed progenitor pools.

### Case Studies

We applied TIPS across three developmental settings (see [examples](https://github.com/xyang2uchicago/TIPS/tree/main/examples)), each testing a different axis of validity:

- **Hematopoiesis** — *We ask: does this score predict a real, ground-truth future outcome?*
  Dataset: [Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F.D. & Klein, A.M. Lineage tracing on transcriptional landscapes links state to fate during differentiation. *Science* 367 (2020).](https://pubmed.ncbi.nlm.nih.gov/31974159/)

- **Hematoendothelial development** — *We ask: does an independently-inferred module generalize to a completely different dataset?*
  Datasets: Dataset 1 — 12.7k genes in 11k E8.25 cells, with 16 predefined developing mesoderm subtypes ([Ibarra-Soria et al., 2018](https://www.nature.com/articles/s41556-017-0013-z)); Dataset 2 — 10.9k genes of 7,240 developing mesoderm cells collected at embryonic day (E) 8.25 when precursor cells of major organs have been formed ([Pijuan-Sala et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30787436/)).

- **Cardiac development** — *We ask: is TIPS prediction biologically and clinically relevant?*
  Datasets: Dataset 1 ([Ibarra-Soria et al., 2018](https://www.nature.com/articles/s41556-017-0013-z)); Dataset 2 ([Pijuan-Sala et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30787436/)); Dataset 3 — 38.9k genes in 230,786 human embryonic stem cells (hESCs) of 13 clusters ([Elorbany et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35061661/)). MuTrans evaluation is given at https://github.com/MohsenZand/ipsc_cardiomyocyte.

### Where to apply TIPS?
- TIPS connects single-cell instability signals to mechanistically grounded network architecture.
- Apply it to complex developmental trajectories with transient transition states (e.g., cardiogenesis).
- Apply it to disease settings where cells deviate from normal trajectories and occupy transitional or plastic states.


### How to use TIPS analysis? 

1. [**regulated stochasticity tutorial**](https://xyang2uchicago.github.io/TIPS/tutorial/TIPS.html) — "Ordered Yet Selectively Fragile Progenitor Networks": a detailed walkthrough of PPIN construction and topological robustness analysis on one of our key results (E8.25 Mouse Gastrulation, GSE87038, [Pijian-Sala 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038)). It covers the network-topology portion of TIPS (Fig. 3d,e of the manuscript); TF–target edge reweighting, CM/CF lineage-biased module resolution, CHD-program convergence and cross-dataset/disease transfer are demonstrated in `examples/` and described in the manuscript.

### How to install?
To use the newest TIPS package, clone/download this repository:
```
git clone https://github.com/xyang2uchicago/TIPS.git
```

### Acknowledgements
TIPS is made possible by contributions from the following authors: Xinan H Yang, Felix Yu, Horatio Ai, Tinjun Luo, and Mohsen Zand.
