
# Cell-type-specific expression and regulation in atypical Schinzel Giedion Syndrome (SGS)

__Jordan Whitlock, Tabea Soelter, Timothy Howton, Elizabeth Wilk, Vishal Oza, Brittany Lasseigne 2023__


__The University of Alabama at Birmingham (UAB), Heersink School of Medicine__
## Data Availability

[![Zenodo](https://img.shields.io/badge/Zenodo-add_zeondo_here!!!-green)](https://www.biorxiv.org/)
[![GEO](https://img.shields.io/badge/GEO-GSE237816-pink)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816)

All data for this project is publicly avialable on Zenodo or [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237816)

### Citation
[![DOI](https://img.shields.io/badge/DOI-add_doi_here!!!-blue)](https://www.biorxiv.org/)

> This is a block quote placeholder for the citation



### Authors 
<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

- [@jordanwhitlock](https://github.com/jordanwhitlock)

- [@tsoelter](https://github.com/tsoelter)

- [@vishaloza](https://github.com/vishaloza)

- [@tchowton](https://github.com/tchowton)

- [@lizzyjoan](https://github.com/lizzyjoan)

- [@blasseigne](https://github.com/blasseigne) [(Lasseigne Lab)](https://www.lasseigne.org/)

## Overview
{% figure caption:"Schematic overview of our approach (A) We generated single-nuclei RNA-seq (snRNA-seq) from male C57BL6/J (WT) and Setbp1S858R cerebral cortex and kidney tissues. (B) We processed and aggregated our data across samples to create cell-type-specific count matrices. (C) Next, we assessed cell-type-specific expression for Setbp1 and genes SETBP1 is known to target. (D) With decoupleR, we measured cell-type-specific TF activity. (E) We acquired protein-protein interaction data (PPI; green) from STRING and Transcription Factor-motif data (TF-motif; yellow) from CIS-BP enriched for SETBP1 targets. We then built cell-type-specific TF-gene regulatory networks with the data we generated by using the message-passing algorithm PANDA to identify regulatory relationships between TFs (circles), genes (rectangles), and proteins (triangles) for each S858R and WT cell type in both tissues and investigated differential communities using ALPACA, differential gene targeting, cooperativity analysis, and functional enrichment analysis." %}
    ![Fig1 (8)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/866508b6-aa6c-4380-83d2-be6142d40f27)
{% endfigure %}

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


