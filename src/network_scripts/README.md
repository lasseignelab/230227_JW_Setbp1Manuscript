README
================
Jordan Whitlock
2023-05-16

# Network Construction and Analysis

## Purpose:

All scripts in this directory are used to build necessary inputs
required for network construction as well as downstream analyses
investigating TF activity, differential targeting, cell-specific
community structures, etc.

### Reproducibility:

**Network Construction** \* PANDA jobs were run using UAB Cheaha
Supercomputer and it’s SLURM scheduler. Jobs were run as an array,
within a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript_panda_1.0.1/general)
container converted Singularity in order to execute. \* Input
construction and all other downstream analyses were carried out in a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/setbp1_manuscript/general)
as well.

### Scripts:

#### PANDA Input Construction:

This directory contains all scripts needed to construct the
protein-protein interaction (ppi), transcription factor motif (TF-motif)
and cell-type specific expression matrices required for
[PANDA](https://netzoo.github.io/zooanimals/panda/) regulatory networks
for both human and mouse.

    ## PANDA_input_construction
    ## +-- 01_Mouse_TFmotif.Rmd
    ## +-- 02_Human_TFmotif.Rmd
    ## +-- 03_HumanMouse_TF_motif_enrichment.Rmd
    ## +-- 04_HumanMouse_ppi.Rmd
    ## \-- 05_MouseSetbp1_expression.Rmd

#### decoupleR Input Construction:

This directory contains all scripts needed for setting up required
inputs for [decoupleR](https://saezlab.github.io/decoupleR/). In order
to measure TF activity, a prior network containing regulatory
relationships between TF and genes for all cell types per condition for
both tissues was built using PANDA. These networks were restructured
here to be a 3 column matrix containing “source” (TFs), “target”
(genes), and “mor” (edge weight). In addition to a prior network,
decoupleR also needs an expression input in the format of a gene x cell
matrix for each cell type and each condition.

    ## decoupleR_input_construction
    ## +-- 01_decoupleR_formatting_prior.Rmd
    ## \-- 02_MouseSetbp1_decoupleR_expression.Rmd

#### decoupleR PANDA:

Gene regulatory networks were constructed using the scripts in this
directory. In addition, the *.err* and *.out* files for each array job
are included here to provide detailed information on the jobs. Prior to
building PANDA networks, the user must first construct the inputs needed
(ppi, TF-motif, and expression) and then build the *.txt*, which is a
list of expression input file paths needed for the array job. Here we
constructed a prior network for decoupleR containing all cell types for
each tissue and both conditions.

    ## decoupleR_PANDA
    ## +-- 01_MouseSetbp1_decoupleR_priorNetwork_expression.Rmd
    ## +-- 02_array_construction.R
    ## +-- 03_PANDA_array.sh
    ## +-- PANDA.R
    ## +-- dePANDA_20235539_0.err
    ## +-- dePANDA_20235539_0.out
    ## +-- dePANDA_20235539_1.err
    ## +-- dePANDA_20235539_1.out
    ## +-- dePANDA_20235539_2.err
    ## +-- dePANDA_20235539_2.out
    ## +-- dePANDA_20235539_3.err
    ## \-- dePANDA_20235539_3.out

#### PANDA:

Gene regulatory networks were constructed using the scripts in this
directory. In addition, the *.err* and *.out* files for each array job
are included here to provide detailed information on the jobs. Prior to
building PANDA networks, the user must first construct the inputs needed
(ppi, TF-motif, and expression) and then build the *.txt*, which is a
list of expression input file paths needed for the array job.

    ## PANDA
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA_array.sh
    ## \-- PANDA.R

#### Differential Targeting:

Gene targeting analyses were performed across and within cell types in
each respective condition as well as between conditions, which is
referred to as differential targeting. All scripts carrying out these
analysis in order to determine if there is a cell-type specific
enrichment in either tissue for gene targeting within the S858R,
Schinzel-Giedion Syndrome mice are within this directory.

    ## differential_targeting
    ## +-- 01_Setbp1_DiffTargetingWithin_Cortex.Rmd
    ## +-- 01_Setbp1_DiffTargetingWithin_Kidney.Rmd
    ## +-- 02_Setbp1_DiffTargetingAcross_Cortex.Rmd
    ## \-- 02_Setbp1_DiffTargetingAcross_Kidney.Rmd

#### Community Detection:

In order to understand differences in regulation across communities at
the cellular level, we detected communities for each cell specific
network using [CONDOR](https://netzoo.github.io/zooanimals/condor/) and
found differential communities between S858R and WT mice using
[ALPACA](https://netzoo.github.io/zooanimals/alpaca/). Scripts for this
analysis are found here

    ## community_detection
    ## \-- Setbp1_Communities_02.Rmd
