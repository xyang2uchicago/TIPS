# TIPS: an R-package for Characterization of Biological Tipping Points
### What is TIPS?
For Holly

<!-- <img src="https://github.com/xyang2uchicago/BioTIP/blob/master/results/Fig1_BioTIP_github.jpg">  -->

### Why TIPS?
For Holly

#### 3 Case Studies

We applied TIPS to three datasets (see [examples](https://github.com/xyang2uchicago/TIPS/tree/main/examples)). The three datasets are as follows:

1. Dataset 1 consists of 96 selected genes in 929 human embryonic stem cells (hESCs) of predefined 9 clusters ([Bargaje et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28167799/)). 
2. Dataset 2 consists of 10.9k genes of 7,240 developing mesoderm cells collected at embryonic day (E) 8.25 when precursor cells of major organs have been formed ([Pijuan-Sala et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30787436/)). 
3. Dataset 3 contained 12.7k genes in 11k E8.25 cells, with 16 predefined developing mesoderm subtypes ([Ibarra-Soria et al., 2018](https://www.nature.com/articles/s41556-017-0013-z)). 


### Where to apply TIPS?
For Holly

### How to use TIPS analysis? 


1. [TIPS tutorial](https://htmlpreview.github.io/?https://github.com/xyang2uchicago/TIPS/blob/main/tutorial/TIPS.html): This is a detailed walkthrough of TIPS on one of our key results (Mouse Gastrulation, GSE87038, [E8.25 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038)). 

<!-- 2. [Vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioTIP/inst/doc/BioTIP.html): This documented exampled case studies on bulk (GSE6136) and single-cell ([Nestorowa 2016](https://pubmed.ncbi.nlm.nih.gov/27365425/)) datasets. -->

### How to install?
To use the newest TIPS package, either clone/download this repository, or you can install TIPS with:

Will add when Bioconductor Package is created
<!-- ```r
library("devtools")
devtools::install_github("xyang2uchicago/TIPS")
```

You can install the released version of BioTIP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BioTIP")
```
or even better:
``` r
source('http://bioconductor.org/biocLite.R')
biocLite("BioTIP")
``` -->

### Variable Clarification

markers.up_all_ttest.rds – List containing upregulated genes per cluster determined by positive log fold change, with no minimum expression proportion filter applied.

markers.up_ttest_min.prop0.25.rds – List containing upregulated genes per cluster determined by positive log fold change, filtered to include only genes expressed in at least 25% of cells in the cluster.

DEG_perState_min.prop0.25_lfc0.6_FDFR0.05.rds – List containing significant upregulated differentially expressed genes per cluster (we call these HiG genes), filtered to genes expressed in ≥25% of cells, with logFC > 0.6 and FDR < 0.01.

### Acknowledgements
TIPS is made possible by contributions from the following authors: Xinan H Yang and Felix Yu
