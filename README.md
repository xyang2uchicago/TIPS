# TIPS: an R-package for Characterization of Biological Tipping Points

For more details, refer to [**TIPS tutorial**](https://htmlpreview.github.io/?https://github.com/xyang2uchicago/TIPS/blob/main/tutorial/TIPS.html)

### What is TIPS?
- TIPS stands for **Transcriptional Instability Prediction of Subnetworks**.
- TIPS quantifies developmental transitions by linking transcriptional instability to interaction-network structure (Fig 1A).
- TIPS couples tipping-point detection (BioTIP)[https://github.com/xyang2uchicago/BioTIP] with state-specific, coexpression-weighted protein–protein interaction networks and robustness analyses (Fig 1B).
<br>
<p align="center">
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_pipeline.jpg" width="50%" height="auto">
<br>
<br>
<br>
<img src="https://github.com/xyang2uchicago/TIPS/blob/main/imgs/Fig3_TIPS_analysis_strategy.jpg" width="50%" height="auto">
</p>
<br>

### Why TIPS?
- Vulnerable progenitor states are transient, heterogeneous, and difficult to capture using markers alone.
- Pinpointing when—and how—perturbations rewire developmental programs remains challenging.
- TIPS identifies defect-vulnerable windows and prioritizes candidate regulators during fate
transitions.

### 3 Case Studies

We applied TIPS to three datasets (see [examples](https://github.com/xyang2uchicago/TIPS/tree/main/examples)). The three datasets are as follows:

1. Dataset 1 consists of 96 selected genes in 929 human embryonic stem cells (hESCs) of predefined 9 clusters ([Bargaje et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28167799/)). MuTrans evaluation is given at https://github.com/MohsenZand/ipsc_cardiomyocyte.
2. #### Dataset 2 consists of 10.9k genes of 7,240 developing mesoderm cells collected at embryonic day (E) 8.25 when precursor cells of major organs have been formed ([Pijuan-Sala et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30787436/)).
3. #### Dataset 3 contained 12.7k genes in 11k E8.25 cells, with 16 predefined developing mesoderm subtypes ([Ibarra-Soria et al., 2018](https://www.nature.com/articles/s41556-017-0013-z)).

### Where to apply TIPS?
- TIPS connects single-cell instability signals to mechanistically grounded network architecture.
- It is applicable to cardiac developmental and other organogenesis trajectories involving transitional cell states.

### How to use TIPS analysis? 

1. [**TIPS tutorial**](https://htmlpreview.github.io/?https://github.com/xyang2uchicago/TIPS/blob/main/tutorial/TIPS.html): This is a detailed walkthrough of TIPS on one of our key results (Mouse Gastrulation, GSE87038, [E8.25 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038)). 

### How to install?
To use the newest TIPS package, clone/download this repository:
```
git clone https://github.com/xyang2uchicago/TIPS.git
```

### Acknowledgements
TIPS is made possible by contributions from the following authors: Xinan H Yang and Felix Yu
