
# Cell-type-specific expression and regulation in atypical Schinzel Giedion Syndrome (SGS)

__Jordan Whitlock, Tabea Soelter, Timothy Howton, Elizabeth Wilk, Vishal Oza, Brittany Lasseigne 2023__


__The University of Alabama at Birmingham (UAB), Heersink School of Medicine__
## Data Availability

All data for this project is publicly available:

* __GEO:__ [![GEO](https://img.shields.io/badge/GEO-GSE237816-pink)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816)

* __SRA:__ [![SRA](https://img.shields.io/badge/SRA-PRJNA996862-purple)](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA996862&o=acc_s%3Aa)

* __Repository:__ [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8190948.svg)](https://doi.org/10.5281/zenodo.8190948)

* __Docker:__ [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8190923.svg)](https://doi.org/10.5281/zenodo.8190923)

* __Data:__ [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8192482.svg)](https://doi.org/10.5281/zenodo.8192482)

* __Cell-type-specific PANDA Networks:__
   * __kidney__: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8199697.svg)](https://doi.org/10.5281/zenodo.8199697)
   * __cerebral cortex__: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8199844.svg)](https://doi.org/10.5281/zenodo.8199844)

### Citation
[![DOI](https://img.shields.io/badge/DOI-JCMM-blue)](http://doi.org/10.1111/jcmm.18001)

> Whitlock, Jordan H., Tabea M. Soelter, Timothy C. Howton, Elizabeth J. Wilk, Vishal H. Oza, and Brittany N. Lasseigne. 2023. “Cell-Type-Specific Gene Expression and Regulation in the Cerebral Cortex and Kidney of Atypical Setbp1S858R Schinzel Giedion Syndrome Mice.” Journal of Cellular and Molecular Medicine.



### Authors 
<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

[The Lasseigne Lab](https://www.lasseigne.org/)

- [@jordanwhitlock](https://github.com/jordanwhitlock)

- [@tsoelter](https://github.com/tsoelter)

- [@vishaloza](https://github.com/vishaloza)

- [@tchowton](https://github.com/tchowton)

- [@lizzyjoan](https://github.com/lizzyjoan)

- [@blasseigne](https://github.com/blasseigne) 

## Overview
![Copy of Fig1 (5)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/b120d83e-75f4-4bf4-a036-ad78e06828da)
   
## Approach

This repository provides a framework for investigating the cell-type-specific impact of genetic variants on gene expression and regulation.  

Here we provide code and data used to investigate SETBP1’s role as an epigenetic hub contributing to cell-type-specific differences in expression, TF activity, gene targeting, and regulatory rewiring: 

* Process 10X single-nuclei RNA-sequencing data using [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
* Perform Quality Control and process data with [Seurat](https://satijalab.org/seurat/)
* Assess TF Activity of SETBP1 and other TFs of interest ([decoupleR](https://saezlab.github.io/decoupleR/))
* Construct cell-type-specific bi-partite TF-gene regulatory networks using message passing algorithm [PANDA](https://netzoo.github.io/zooanimals/panda/)
* Carry out downstream network analyses 
    * Differential Community Detection ([ALPACA](https://netzoo.github.io/zooanimals/alpaca/))
    * Differential Gene Targeting
    * Network Rewiring
## Funding 

This work was supported in part by the UAB Lasseigne Lab funds, UAB Pilot Center for Precision Animal Modeling (C-PAM)(1U54OD030167), the UAB Predoctoral Training Grant in Cell, Molecular, and Developmental Biology (CMDB T32)(5T32GM008111-35
### License
[![License](https://img.shields.io/badge/LICENSE-MIT_License-yellow)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/LICENSE) 

This repository is licensed under the MIT License, see LICENSE documentation within this repository for more details.


