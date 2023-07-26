README
================
Jordan Whitlock
2023-07-26

# snRNA-seq quality control and pre-processing

## Purpose:

All scripts in `seurat_scripts`were used for Quality control and
processing of snRNA-seq data for S858R and WT cerebral cortex and kidney
tissues.

### Reproducibility:

Scripts were carried out in a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript/general)
using R version 4.1.3.

### Scripts:

![Copy of Fig1 (2)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/9931bd70-5094-48b2-90cd-ebdbc00cc51a)


Detailed information on the goal and purpose can be found within each
processing script.

    ## seurat_scripts
    ## +-- 01_Setbp1_QC.Rmd
    ## +-- 02_Setbp1_Integration.Rmd
    ## +-- 03_Setbp1_Clustering.Rmd
    ## +-- 04_Setbp1_MarkersCortex.Rmd
    ## +-- 04_Setbp1_MarkersKidney.Rmd
    ## +-- 05_Setbp1_CellTypesCortex.Rmd
    ## +-- 05_SoupX_Setbp1_CellTypesKidney.Rmd
    ## +-- 06_SoupX_Setbp1_AmbientRNAKidney.Rmd
    ## +-- 07_Setbp1_QC_postSoup.Rmd
    ## +-- 08_Setbp1_Integration_postSoup.Rmd
    ## +-- 09_Setbp1_Clustering_postSoup.Rmd
    ## +-- 10_Setbp1_MarkersKidney_postSoup.Rmd
    ## +-- 11_Setbp1_CellTypesKidney_postSoup.Rmd
    ## +-- 12_Setbp1_DGEcortex.Rmd
    ## +-- 12_Setbp1_DGEkidney.Rmd
    ## +-- 13_Setbp1_PathwayAnalysisCortexKidney.Rmd
    ## +-- 14_reactive_astrocytes.Rmd
    ## +-- 15_Set_Sox2_profiling.Rmd
    ## +-- 16_failed_repair_kidney_profiling.Rmd
    ## \-- Setbp1_target_list_construction.Rmd

# Network Construction and Analysis

## Purpose:

All scripts in `network_scripts` were used to build necessary inputs
required for network construction as well as downstream analyses
investigating TF activity, differential targeting, differential
cell-type-specific community structures, etc.

### Reproducibility:

**Network Construction and Tf Activity Analysis**

-   PANDA jobs were run using UAB Cheaha Supercomputer, and its SLURM
    scheduler. Jobs were run as an array within a
    [Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript_panda_1.0.1/general)
    container converted Singularity in order to execute.

-   Input construction and all other downstream analyses were carried
    out in a
    [Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript/general)
    as well.

### Scripts:

--------

### Part 1:
![Untitled (3)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/e2e73c4a-eb26-4b96-9847-0b81004bf1a7)


#### decoupleR Input Construction and TF Activity Analysis:

This directory contains all scripts needed for setting up required
inputs for [decoupleR](https://saezlab.github.io/decoupleR/). In order
to measure TF activity, a prior network containing regulatory
relationships between TF and genes for all cell types per condition for
both tissues were constructed using
[CollecTRI](https://github.com/saezlab/CollecTRI) (accessed on 230528,
prior is provided with this repository for both human and mouse). More
information
[here](https://www.biorxiv.org/content/10.1101/2023.03.30.534849v1).
Additionally, a script was included for users to construct a new prior
using CollecTRI. The prior network for each species was restructured to
be a 3-column matrix containing “source” (TFs), “target” (genes), and
“mor” (edge weight). In addition to a prior network, decoupleR also
needs an expression input in the gene x cell matrix format for each
tissue and condition. In addition, the *.err* and *.out* files for each
array job are included here to provide detailed information on the jobs.

    ## network_scripts/decoupleR_input_construction
    ## +-- 01_decoupleR_formatting_prior.Rmd
    ## \-- 02_MouseSetbp1_decoupleR_expression.Rmd

    ## network_scripts/decoupleR
    ## +-- 01_MouseSetbp1_decoupleR_inputs.Rmd
    ## +-- 02_decoupleR_analysis.R
    ## +-- 02_decoupleR_array_job.sh
    ## +-- 03_TF_activity_cortex.Rmd
    ## \-- 04_TF_activity_kidney.Rmd

1.  [01_MouseSetbp1_decoupleR_inputs.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/decoupleR/01_MouseSetbp1_decoupleR_inputs.Rmd)
2.  [02_decoupleR_analysis.R](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/decoupleR/02_decoupleR_analysis.R)
    and accompanying bash script
    [02_decoupleR_array_job.sh](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/decoupleR/02_decoupleR_array_job.sh)
3.  TF activity analysis scripts
    [03_TF_Activity_cortex.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/decoupleR/03_TF_activity_cortex.Rmd)
    and
    [03_TF_Activity_kidney.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/decoupleR/04_TF_activity_kidney.Rmd)

--------

### Part 2
![Untitled (2)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/c0c5f7e3-6284-4c82-be6b-512a97bc94a1)

#### PANDA Input Construction:

This directory contains all scripts needed to construct the
protein-protein interaction (PPI), transcription factor motif (TF-motif)
and cell-type specific expression matrices required for
[PANDA](https://netzoo.github.io/zooanimals/panda/) regulatory networks
for mouse.

    ## network_scripts/PANDA_input_construction
    ## +-- 01_Mouse_TFmotif.Rmd
    ## +-- 02_Mouse_TF_motif_enrichment.Rmd
    ## +-- 03_Mouse_ppi.Rmd
    ## \-- 04_MouseSetbp1_expression.Rmd

#### PANDA:

Gene regulatory networks were constructed using the scripts in this
directory. In addition, the *.err* and *.out* files for each array job
are included here to provide detailed information on the jobs. Before
building PANDA networks, the user must first construct the inputs needed
(ppi, TF-motif, and expression) and then build the *.txt*, a
list of expression input file paths needed for the array job.

    ## network_scripts/PANDA
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA.R
    ## \-- 02_PANDA_array.sh

#### Differential Targeting:

Gene targeting analyses were performed across and within cell types in
each respective condition as well as between conditions, which is
referred to as differential targeting. All scripts carrying out these
analysis in order to determine if there is a cell-type specific
enrichment in either tissue for gene targeting within the S858R,
Schinzel-Giedion Syndrome mice are within this directory.

    ## network_scripts/differential_targeting
    ## +-- 01_Setbp1_DiffTargetingWithin_Cortex.Rmd
    ## +-- 02_Setbp1_DiffTargetingWithin_Kidney.Rmd
    ## \-- 03_Setbp1_DiffTargeting_FEA.Rmd

1.  [01_Setbp1_DiffTargetingWithing_Cortex.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/differential_targeting/01_Setbp1_DiffTargetingWithin_Cortex.Rmd)
    and
    [01_Setbp1_DiffTargetingWithing_Kidney.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/differential_targeting/02_Setbp1_DiffTargetingWithin_Kidney.Rmd):
    code used to calculate differential gene targeting within cell types
    between conditions for both tissues

2.  Functional Enrichment Analysis (FEA) on differentially targeted
    genes
    [03_Setbp1_DiffTargeting_FEA.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/differential_targeting/03_Setbp1_DiffTargeting_FEA.Rmd)

#### Community Detection:

In order to understand differences in regulation across communities at
the cellular level, we detected communities for each cell-type-specific
network using [CONDOR](https://netzoo.github.io/zooanimals/condor/) and
found differential communities between S858R and WT mice using
[ALPACA](https://netzoo.github.io/zooanimals/alpaca/). Scripts for this
analysis is found here

    ## network_scripts/community_detection
    ## +-- 01_alpaca_array_job.sh
    ## +-- 01_alpaca_networks_array.R
    ## +-- 02_Setbp1_Community_Identification.Rmd
    ## +-- 03_alpaca_network_analysis.Rmd
    ## \-- file_pairs.txt

The ALPACA analysis consists of the following scripts:

1.  [file_pairs.txt](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/file_pairs.txt):
    contains a list of all paired comparisons for differential community
    detection analysis. One comparison per line containing the exact
    name of the PANDA objects you are wanting to use for ALPACA. See
    the example below of setting up the file pairs to identify cell
    type-specific community differences between S858R and WT samples.

astrocytes_controlcortexexpression_PANDA.Rdata
astrocytes_heterozygouscortexexpression_PANDA.Rdata

2.  [01_alpaca_array_job.sh](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/alpaca_array_job.sh):
    example bash script used for submitting ALPACA to UAB Cheaha
    Supercomputer SLURM scheduler. This script calls on the file pairs
    specified in \#1 and the accompanying R script and submits each as
    an array job with the specified partitions and memory.

3.  [01_alpaca_networks_array.R](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/alpaca_networks_array.R):
    code to run the actual analysis to identify differential communities
    using CONDOR and ALPACA. The `.Rdata` file for the communities
    constructed for every cell type is saved in a `membership/`
    directory and all other output files (`_ctrl_memb.txt`,
    `_final_memb.txt`, `_scores.txt` and `_DWBM.txt`) can be found in
    `main/`.

4.  [02_Community_Identification.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/02_Setbp1_Community_Identification.Rmd):
    code to identify the community number of Setbp1 in order to quantify
    overlap of community structures across cell types.

5.  [03_alpaca_network_analysis.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/03_alpaca_network_analysis.Rmd):
    code for calculating and visualizing Jaccard Similarity Index (JI)
    between Setbp1 differential Communities in S858R cell types.

#### Network Rewiring:

While little is known about the role of specific patient variants,
several mechanisms are hypothesized to contribute to altered
neurodevelopment in SGS, including chromatin remodeling, disrupted cell
cycle control, increased DNA damage, and modified PP2A complex activity.
To further evaluate the molecular impact of the S858R variant, we
investigated the magnitude and direction of the regulatory network and
changes in cooperativity network edge weights through quantifying
proteins working together (i.e., cooperating), where SETBP1 is acting as
a TF on its known target genes. Similar to the regulatory network, a
more positive edge weight between two proteins indicates a greater
likelihood they cooperate together. In contrast, a more negative edge
weight indicates greater confidence that the two proteins do not
cooperate. We used these regulatory or cooperativity edge weight
magnitude and direction changes to infer potential regulatory rewiring
due to the S858R variant.

    ## network_scripts/ppi_rewiring
    ## \-- 01_Setbp1_RegCoopNet.Rmd

1.  [01_Setbp1_RegCoopNet.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/ppi_rewiring/01_Setbp1_RegCoopNet.Rmd)
