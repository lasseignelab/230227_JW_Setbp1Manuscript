
# Cell-type-specific expression and regulation in atypical Schinzel Giedion Syndrome (SGS)

__Jordan Whitlock, Tabea Soelter, Timothy Howton, Elizabeth Wilk, Vishal Oza, Brittany Lasseigne 2023__


__The University of Alabama at Birmingham (UAB), Heersink School of Medicine__
## Data Availability

[![Zenodo](https://img.shields.io/badge/Zenodo-add_zeondo_here!!!-green)](https://www.biorxiv.org/)
[![GEO](https://img.shields.io/badge/GEO-GSE237816-pink)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816)
[![SRA](https://img.shields.io/badge/SRA-PRJNA996862-purple)]([https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA996862&o=acc_s%3Aa))

All data for this project is publicly avialable on Zenodo or [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816)

### Citation
[![DOI](https://img.shields.io/badge/DOI-add_doi_here!!!-blue)](https://www.biorxiv.org/)

> This is a block quote placeholder for the citation



### Authors 
<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">
[__The Lasseigne Lab:__](https://www.lasseigne.org/)

- [@jordanwhitlock](https://github.com/jordanwhitlock)

- [@tsoelter](https://github.com/tsoelter)

- [@vishaloza](https://github.com/vishaloza)

- [@tchowton](https://github.com/tchowton)

- [@lizzyjoan](https://github.com/lizzyjoan)

- [@blasseigne](https://github.com/blasseigne) 

## Overview
![Copy of Fig1 (1)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/330fce42-bf88-4c88-9bc3-cc433bfd8272)
   
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


