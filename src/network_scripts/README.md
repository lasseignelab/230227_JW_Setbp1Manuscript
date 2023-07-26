README
================
Jordan Whitlock
2023-07-26

# Network Construction and Analysis

## Purpose:

All scripts in this directory are used to build necessary inputs
required for network construction as well as downstream analyses
investigating TF activity, differential targeting, cell-specific
community structures, etc.

### Reproducibility:

**Network Construction and Tf Activity Analysis** \* PANDA jobs were run
using UAB Cheaha Supercomputer and it’s SLURM scheduler. Jobs were run
as an array, within a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript_panda_1.0.1/general)
container converted Singularity in order to execute. \* Input
construction and all other downstream analyses were carried out in a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript/general)
as well.

### Scripts:

![Untitled (1)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/225dcd2f-f935-4656-ac13-816f7035b4ab)


#### decoupleR Input Construction and TF Activity Analysis:

This directory contains all scripts needed for setting up required
inputs for [decoupleR](https://saezlab.github.io/decoupleR/). In order
to measure TF activity, a prior network containing regulatory
relationships between TF and genes for all cell types per condition for
both tissues was constructed using
[CollecTRI](https://github.com/saezlab/CollecTRI) (accessed on 230528,
prior is provided with this repository for both human and mouse). More
information
[here](https://www.biorxiv.org/content/10.1101/2023.03.30.534849v1).
Additionally a script was included for users to constructed a new prior
using CollecTRI. The prior network for each species was restructured to
be a 3 column matrix containing “source” (TFs), “target” (genes), and
“mor” (edge weight). In addition to a prior network, decoupleR also
needs an expression input in the format of a gene x cell matrix for each
cell type and each condition. In addition, the *.err* and *.out* files
for each array job are included here to provide detailed information on
the jobs.

    ## decoupleR
    ## +-- 01_MouseSetbp1_decoupleR_inputs.Rmd
    ## +-- 02_decoupleR_analysis.R
    ## +-- 02_decoupleR_array_job.sh
    ## +-- 03_TF_activity_cortex.Rmd
    ## \-- 04_TF_activity_kidney.Rmd

![Untitled (2)](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/assets/62023125/b07d51a4-b5a0-4a29-9f50-6ecae8aff092)

#### PANDA Input Construction:

This directory contains all scripts needed to construct the
protein-protein interaction (ppi), transcription factor motif (TF-motif)
and cell-type specific expression matrices required for
[PANDA](https://netzoo.github.io/zooanimals/panda/) regulatory networks
for both human and mouse.

    ## PANDA_input_construction
    ## +-- 01_Mouse_TFmotif.Rmd
    ## +-- 02_Mouse_TF_motif_enrichment.Rmd
    ## +-- 03_Mouse_ppi.Rmd
    ## \-- 04_MouseSetbp1_expression.Rmd

#### PANDA:

Gene regulatory networks were constructed using the scripts in this
directory. In addition, the *.err* and *.out* files for each array job
are included here to provide detailed information on the jobs. Prior to
building PANDA networks, the user must first construct the inputs needed
(ppi, TF-motif, and expression) and then build the *.txt*, which is a
list of expression input file paths needed for the array job.

    ## PANDA
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

    ## differential_targeting
    ## +-- 01_Setbp1_DiffTargetingWithin_Cortex.Rmd
    ## +-- 02_Setbp1_DiffTargetingWithin_Kidney.Rmd
    ## \-- 03_Setbp1_DiffTargeting_FEA.Rmd

#### Community Detection:

In order to understand differences in regulation across communities at
the cellular level, we detected communities for each cell specific
network using [CONDOR](https://netzoo.github.io/zooanimals/condor/) and
found differential communities between S858R and WT mice using
[ALPACA](https://netzoo.github.io/zooanimals/alpaca/). Scripts for this
analysis are found here

    ## community_detection
    ## +-- 01_alpaca_array_job.sh
    ## +-- 01_alpaca_networks_array.R
    ## +-- 02_Setbp1_Community_Identification.Rmd
    ## +-- 03_alpaca_network_analysis.Rmd
    ## \-- file_pairs.txt

The ALPACA analysis consists of the following scripts:

1.  [file_pairs.txt](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/community_detection/file_pairs.txt):
    contains a list of all paired comparisons for differential community
    detections analysis. One comparison per line containing the exact
    name of the PANDA objects you are wanting to use for ALPACA. See
    below example of how to set up the file pairs to identify cell
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
    using CONDOR and ALPACA. The .Rdata file for the ccommunities
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

    ## ppi_rewiring
    ## \-- 01_Setbp1_RegCoopNet.Rmd

1.  [01_Setbp1_RegCoopNet.Rmd](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/ppi_rewiring/01_Setbp1_RegCoopNet.Rmd)
