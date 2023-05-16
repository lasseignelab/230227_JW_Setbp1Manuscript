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
    ## +-- PANDA.R
    ## +-- PANDA_20202459_0.err
    ## +-- PANDA_20202459_0.out
    ## +-- PANDA_20202459_1.err
    ## +-- PANDA_20202459_1.out
    ## +-- PANDA_20202459_10.err
    ## +-- PANDA_20202459_10.out
    ## +-- PANDA_20202459_11.err
    ## +-- PANDA_20202459_11.out
    ## +-- PANDA_20202459_12.err
    ## +-- PANDA_20202459_12.out
    ## +-- PANDA_20202459_13.err
    ## +-- PANDA_20202459_13.out
    ## +-- PANDA_20202459_14.err
    ## +-- PANDA_20202459_14.out
    ## +-- PANDA_20202459_15.err
    ## +-- PANDA_20202459_15.out
    ## +-- PANDA_20202459_16.err
    ## +-- PANDA_20202459_16.out
    ## +-- PANDA_20202459_17.err
    ## +-- PANDA_20202459_17.out
    ## +-- PANDA_20202459_18.err
    ## +-- PANDA_20202459_18.out
    ## +-- PANDA_20202459_19.err
    ## +-- PANDA_20202459_19.out
    ## +-- PANDA_20202459_2.err
    ## +-- PANDA_20202459_2.out
    ## +-- PANDA_20202459_20.err
    ## +-- PANDA_20202459_20.out
    ## +-- PANDA_20202459_21.err
    ## +-- PANDA_20202459_21.out
    ## +-- PANDA_20202459_22.err
    ## +-- PANDA_20202459_22.out
    ## +-- PANDA_20202459_23.err
    ## +-- PANDA_20202459_23.out
    ## +-- PANDA_20202459_24.err
    ## +-- PANDA_20202459_24.out
    ## +-- PANDA_20202459_25.err
    ## +-- PANDA_20202459_25.out
    ## +-- PANDA_20202459_26.err
    ## +-- PANDA_20202459_26.out
    ## +-- PANDA_20202459_27.err
    ## +-- PANDA_20202459_27.out
    ## +-- PANDA_20202459_28.err
    ## +-- PANDA_20202459_28.out
    ## +-- PANDA_20202459_29.err
    ## +-- PANDA_20202459_29.out
    ## +-- PANDA_20202459_3.err
    ## +-- PANDA_20202459_3.out
    ## +-- PANDA_20202459_30.err
    ## +-- PANDA_20202459_30.out
    ## +-- PANDA_20202459_31.err
    ## +-- PANDA_20202459_31.out
    ## +-- PANDA_20202459_32.err
    ## +-- PANDA_20202459_32.out
    ## +-- PANDA_20202459_33.err
    ## +-- PANDA_20202459_33.out
    ## +-- PANDA_20202459_34.err
    ## +-- PANDA_20202459_34.out
    ## +-- PANDA_20202459_35.err
    ## +-- PANDA_20202459_35.out
    ## +-- PANDA_20202459_36.err
    ## +-- PANDA_20202459_36.out
    ## +-- PANDA_20202459_37.err
    ## +-- PANDA_20202459_37.out
    ## +-- PANDA_20202459_38.err
    ## +-- PANDA_20202459_38.out
    ## +-- PANDA_20202459_39.err
    ## +-- PANDA_20202459_39.out
    ## +-- PANDA_20202459_4.err
    ## +-- PANDA_20202459_4.out
    ## +-- PANDA_20202459_40.err
    ## +-- PANDA_20202459_40.out
    ## +-- PANDA_20202459_41.err
    ## +-- PANDA_20202459_41.out
    ## +-- PANDA_20202459_42.err
    ## +-- PANDA_20202459_42.out
    ## +-- PANDA_20202459_43.err
    ## +-- PANDA_20202459_43.out
    ## +-- PANDA_20202459_44.err
    ## +-- PANDA_20202459_44.out
    ## +-- PANDA_20202459_45.err
    ## +-- PANDA_20202459_45.out
    ## +-- PANDA_20202459_46.err
    ## +-- PANDA_20202459_46.out
    ## +-- PANDA_20202459_47.err
    ## +-- PANDA_20202459_47.out
    ## +-- PANDA_20202459_48.err
    ## +-- PANDA_20202459_48.out
    ## +-- PANDA_20202459_49.err
    ## +-- PANDA_20202459_49.out
    ## +-- PANDA_20202459_5.err
    ## +-- PANDA_20202459_5.out
    ## +-- PANDA_20202459_6.err
    ## +-- PANDA_20202459_6.out
    ## +-- PANDA_20202459_7.err
    ## +-- PANDA_20202459_7.out
    ## +-- PANDA_20202459_8.err
    ## +-- PANDA_20202459_8.out
    ## +-- PANDA_20202459_9.err
    ## \-- PANDA_20202459_9.out

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
