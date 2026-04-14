# TIPS: Linking tipping-point instability to interaction networks for transition-state regulators 

TIPS links tipping-point instability in single-cell data to interpretable interaction subnetworks, prioritizing bridge regulators and “dual-pull” lineage programs during fate transitions. For more details, refer to [**TIPS tutorial**](https://xyang2uchicago.github.io/TIPS/tutorial/TIPS.html)

### What is TIPS?
- Name: **Transcriptional Instability Prediction of Subnetworks**.
- Goal: Identify **semi-stable transition states** obscured by program mixing (Fig. 1A).
- Output: Quantify **state-specific network rewiring** from uncommitted progenitors to descendant fates (Fig 1B).
- Method: BioTIP + **coexpression-weighted protein–protein interaction networks + robustness (targeted-attack) + chromatin integration** to resolve lineage-leaning arms (Fig. 1C).
<br>
<p align="center">
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_pipeline.jpg" width="33%" height="auto">
<br>
<br>
<br>
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_analysis_strategy.jpg" width="33%" height="auto">
</p>
<br>

### Why TIPS?
- Vulnerable progenitor states are transient, heterogeneous, and difficult to capture using markers alone.
- Pinpointing when—and how—perturbations rewire developmental programs remains challenging.
- TIPS maps defect-vulnerable windows and prioritizes lineage-selective regulatory programs from mixed progenitor pools.

### 3 Case Studies

We applied TIPS to three datasets regarding early mesoderm specifications (see [examples](https://github.com/xyang2uchicago/TIPS/tree/main/examples)). The three datasets are as follows:
1. #### Dataset 1 contained 12.7k genes in 11k E8.25 cells, with 16 predefined developing mesoderm subtypes ([Ibarra-Soria et al., 2018](https://www.nature.com/articles/s41556-017-0013-z)).
2. #### Dataset 2 consists of 10.9k genes of 7,240 developing mesoderm cells collected at embryonic day (E) 8.25 when precursor cells of major organs have been formed ([Pijuan-Sala et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30787436/)).
3. Dataset 3 consists of 38.9k genes in 230,786 human embryonic stem cells (hESCs) of 13 clusters  ([Elorbany et al., 2022]([https://pubmed.ncbi.nlm.nih.gov/28167799/](https://pubmed.ncbi.nlm.nih.gov/35061661/)). MuTrans evaluation is given at https://github.com/MohsenZand/ipsc_cardiomyocyte.

### Where to apply TIPS?
- TIPS connects single-cell instability signals to mechanistically grounded network architecture.
- Apply it to complex developmental trajectories with transient transition states (e.g., cardiogenesis).
- Apply it to disease settings where cells deviate from normal trajectories and occupy transitional or plastic states.


### How to use TIPS analysis? 

1. [**TIPS tutorial**](https://xyang2uchicago.github.io/TIPS/tutorial/TIPS.html): This is a detailed walkthrough of TIPS on one of our key results (E8.25 Mouse Gastrulation, GSE87038, [Pijian-Sala 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038)). 

### How to install?
To use the newest TIPS package, clone/download this repository:
```
git clone https://github.com/xyang2uchicago/TIPS.git
```

### Acknowledgements
TIPS is made possible by contributions from the following authors: Xinan H Yang and Felix Yu
