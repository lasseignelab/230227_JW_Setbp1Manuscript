FROM bioconductor/bioconductor_docker:RELEASE_3_14

# location for mounting
#RUN mkdir /home/rstudio/data


##----------------------------------------------------------------------------##
##  R PACKAGES
##----------------------------------------------------------------------------##

#### General CRAN Packages ####

RUN R -e 'install.packages(c("tidyverse", "remotes",  "BiocManager", "devtools", "knitr", "lintr", "styler", "here"))'

#### Project Specific CRAN Packages ####

RUN R -e 'install.packages(c("NMF", "SoupX", "fgsea", "gprofiler2", "circlize", "msigdbr", "enrichR", "lattice", "grr", "rafalib", "markdown", "reticulate", "harmony", "fields", "KernSmooth", "ROCR", "paralell", "pheatmap", "Seurat", "ggplot2", "Matrix.utils", "readr", "dplyr", "clustree", "sctransform", "igraph", "reshape2", "visNetwork", "RColorBrewer", "plotly", "data.table", "cowplot", "magrittr", "Matrix", "purrr", "tibble", "png", "hdf5r"))'

#### Project Specific GitHub Packages ####

RUN R -e 'devtools::install_version("Rcpp", version = "1.0.8", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_github("saeyslab/nichenetr")'
RUN R -e 'devtools::install_github("saezlab/liana")'
RUN R -e 'remotes::install_github( repo = "chris-mcginnis-ucsf/DoubletFinder", ref = "master")'
RUN R -e 'remotes::install_github(repo = "theislab/kBET", ref = "master")'
RUN R -e 'remotes::install_github(repo = "netZoo/netZooR", ref = "master")'

#### Project Specific Bioconductor Packages ####

# RUN R -e 'BiocManager::install(c("GEOquery", "rrvgo", "AnnotationHub", "org.Mm.eg.db", "SingleCellExperiment", "ComplexHeatmap", "scran", "limma", "Biobase", "scater", "edgeR", "S4vectors", "apeglm", "DESeq2"))'

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  python-dev \
  python3-pip \
  python3-dev \
  python3-pip \
  python3-tk


RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev libxt-dev \
  libxml2-dev \
  libudunits2-dev \
  libhdf5-dev \
  libv8-dev \
  libgdal-dev \
  xorg libx11-dev libglu1-mesa-dev \
  samtools libboost-all-dev \
  libgsl-dev

##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##

RUN pip3 install leidenalg


