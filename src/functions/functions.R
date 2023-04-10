# Function-scatterplot
scatterPlot <- function(
    simMatrix, reducedTerms, size = "score", addLabel = TRUE,
    labelSize = 3) {
  if (!all(sapply(c("ggplot2", "ggrepel"), requireNamespace,
                  quietly = TRUE
  ))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ",
         "Consider installing them before using this function.",
         call. = FALSE
    )
  }
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)),
                eig = TRUE,
                k = 2
  )
  df <- cbind(as.data.frame(x$points), reducedTerms[match(
    rownames(x$points),
    reducedTerms$go
  ), c("term", "parent", "parentTerm", "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = "V1", y = "V2", color = "parentTerm")) +
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) +
    ggplot2::scale_color_discrete(guide = "none") +
    ggplot2::scale_size_continuous(
      guide = "none",
      range = c(0, 25)
    ) +
    ggplot2::scale_x_continuous(name = "") +
    ggplot2::scale_y_continuous(name = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = "parentTerm"),
                                  data = subset(df, parent == rownames(df)), box.padding = grid::unit(
                                    0.5,
                                    "lines"
                                  ), size = labelSize, max.overlaps = 20
    )
  } else {
    p
  }
}

# function-parallel_paranSweep_v3_noapprox; used in paramSweep_v3_noapprox, setting RunPCA approx = FALSE
parallel_paramSweep_v3_noapprox <- function (n, n.real.cells, real.cells, pK, pN, data, orig.commands, 
                                             PCs, sct) 
{
  sweep.res.list = list()
  list.ind = 0
  print(paste("Creating artificial doublets for pN = ", pN[n] * 
                100, "%", sep = ""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  if (sct == FALSE) {
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                               model.use = orig.commands$ScaleData.RNA$model.use, 
                               do.scale = orig.commands$ScaleData.RNA$do.scale, 
                               do.center = orig.commands$ScaleData.RNA$do.center, 
                               scale.max = orig.commands$ScaleData.RNA$scale.max, 
                               block.size = orig.commands$ScaleData.RNA$block.size, 
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                            npcs = length(PCs), approx = FALSE, rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                            verbose = FALSE)
  }
  if (sct == TRUE) {
    require(sctransform)
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs), approx = FALSE)
  }
  print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                            PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[, 1:n.real.cells]
  print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[, i] <- order(dist.mat[, i])
  }
  ind <- round(nCells * max(pK)) + 5
  dist.mat <- dist.mat[1:ind, ]
  print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, 
                                 ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1
    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1), i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }
    sweep.res.list[[list.ind]] <- pANN
  }
  return(sweep.res.list)
}

#function-paramSweep-v3-noapprox; setting RunpCA approx to false, used in defin_pK function
paramSweep_v3_noapprox <- function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3_noapprox, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3_noapprox, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

# Function-define-pK-noapprox; used in doubletfinder script, setting approx for RunPCA = FALSE
define_pK_noapprox <- function(sample){
  sweep_res_list <- paramSweep_v3_noapprox(sample, PCs = 1:30, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  pK <- as.numeric(as.character(bcmvn$pK))
  BCmetric <- bcmvn$BCmetric
  pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
  
  par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)
  plot(x = pK, y = BCmetric, pch = 16, type = "b",
       col = "blue", lty = 1)
  abline(v = pK_choose, lwd = 2, col = "red", lty = 2)
  title("The BCmvn distributions")
  text(pK_choose, max(BCmetric), as.character(pK_choose), pos = 4, col = "red")
}

# function doubletFinder_v3_noapprox ; used in doublet_removal, setting RunPCA approx = FALSE
doubletFinder_v3_noapprox <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                       sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), approx = FALSE, rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs), approx = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbor_types[i, ] <- table(doublet_types1[neighbors - 
                                                        n_real.cells]) + table(doublet_types2[neighbors - 
                                                                                                n_real.cells])
          neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, 
          ])
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}

# function doublet_removal used in doubletFinder script; setting RunPCA approx = FALSE
doublet_removal_noapprox <- function(sample, expected, pK){
  # define the expected number of doublets in nuclei.
  nExp <- round(ncol(sample) * expected)  # expect 10.2% doublets
  sample <- doubletFinder_v3_noapprox(sample, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:30)
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name <- colnames(sample@meta.data)[grepl("DF.classification", colnames(sample@meta.data))]
  plot1 <- cowplot::plot_grid(ncol = 2, DimPlot(sample, group.by = "orig.ident") + NoAxes(),
                              DimPlot(sample, group.by = DF.name) + NoAxes())
  plot(plot1)
  VlnPlot(sample, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  #REMOVE DOUBLETS:
  sample <- sample[, sample@meta.data[, DF.name] == "Singlet"]
  assign("no_doublets", sample, envir = globalenv())
  #revisualize:
  plot2 <- cowplot::plot_grid(ncol = 2, DimPlot(sample, group.by = "orig.ident") + NoAxes(),
                              DimPlot(sample, group.by = DF.name) + NoAxes())
  plot(plot2)
}

# aggregate matrix function
##Matrix.utils function not loading into Docker, so grabbed source code (https://rdrr.io/cran/Matrix.utils/src/R/Matrix.utils.R)
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

# dmCast function
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
      if(nlevels(x)!=length(na.omit(unique(x))))
        x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })
  
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

# Functions for PANDA network job for loop
#make function for loading .Rdata to reassign to variable in loop 
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#run 'makePanda' function to make PANDA networks
makePanda <- function(motif, ppi, expression){
  
  #remove expression values = 0 
  #expression <- expression[rowSums(expression[])>0,] #removing values with zeroes, this is commented out because removing them here will cause different sized matrices and result in downstream network comparisons not being possible
  
  #run PANDA
  panda(expr = expression, ppi = ppi, motif = motif, progress = TRUE, mode = "intersection") #running in default mode 
  
}

# Functions for Setbp1_AllCortex_PANDAComparison_positive_01

## function-for-functional-enrichment-pathway-analysis
fea <- function(genes, organism, max_term = 1000, min_term = 5){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  fea_result_filt <- fea_result$result %>% dplyr::filter(., term_size < max_term & term_size > min_term) #std for max_term is 1000 and min_term is 10
  # # select the top ___ pathways for plotting --- suggest 50, could do more, but difficult to see visually
  # fea_result_filt <- fea_result %>% top_n(n = pathway_number)
  return(fea_result_filt)
}

fea_no_sig <- function(genes, organism, max_term = 1000, min_term = 5){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       FALSE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  #fea_result_filt <- fea_result
  # remove arbitrary pathways --- do not want pathways too "generic"
  #\fea_result_filt <- fea_result$result %>% dplyr::filter(., term_size < max_term & term_size > min_term) #std for max_term is 1000 and min_term is 10
  # # select the top ___ pathways for plotting --- suggest 50, could do more, but difficult to see visually
  # fea_result_filt <- fea_result %>% top_n(n = pathway_number)
  return(fea_result_filt)
}

## function-for-plotting-functional-enrichment-pathway-analysis
bubbleplot <- function(fea_result_filt){
  plot <- ggplot(fea_result_filt, aes(x = intersection_size, y = reorder(term_name, -p_value), size = recall, fill =
                                        p_value)) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_size(range = c(2, 10), name = "# Genes Matched to Term") + 
    scale_fill_distiller(palette = "Purples") + 
    labs(x = "Intersection Size", y = "Functional Enrichment Terms")
  return(plot)
}


## function-GO-positive-negative-combined
combined_fea <- function(genes, organism, max_term = 1000, min_term = 5){
  posgenes <- as.list(genes %>% filter(value == "positive"))
  pos_fea_filt <- fea(genes= posgenes, organism = organism, max_term = max_term, min_term = min_term) %>% mutate(direction = "positive")
  neggenes <- as.list(genes %>% filter(value == "negative"))
  neg_fea_filt <- fea(genes = neggenes, organism = organism, max_term = max_term, min_term = min_term) %>% mutate(direction = "negative")
  combined_filt <- rbind(pos_fea_filt, neg_fea_filt)
  combined_fea_result <- list(control_enriched = neg_fea_filt, condition_enriched = pos_fea_filt, combined = combined_filt)
  return(combined_fea_result)
}

#function DGE fea
fea_DGE <- function(genes, organism){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  fea_result <- fea_result$result %>% filter(term_size < 1000 | term_size > 10)
  # select the top 50 pathways for plotting --- could do more, but difficult to see visually
  fea_result_filt <- fea_result %>% top_n(n = 50)
  return(fea_result_filt)
}

combined_fea_DGE <- function(genes, organism){
  upgenes <- as.list(genes %>% filter(diffexpressed == "UP"))
  up_fea_filt <- fea(genes = upgenes, organism = organism) %>% mutate(direction = "upregulated")
  downgenes <- as.list(genes %>% filter(diffexpressed == "DOWN"))
  down_fea_filt <- fea(genes = downgenes, organism = organism) %>% mutate(direction = "downregulated")
  combined_filt <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = combined_filt)
  return(fea_result)
}

#setbp1 gene set enrichment functions side by side fea
fea_set <- function(genes, organism){
  # create gprofiler2 query ---
  fea_result <- gost(query = genes, organism = organism, ordered_query = FALSE, multi_query = FALSE, significant = 
                       FALSE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 
                       0.05, correction_method = "bonferroni", domain_scope = "annotated", numeric_ns =
                       "", sources = NULL, as_short_link = FALSE) 
  # remove arbitrary pathways --- do not want pathways too "generic"
  #fea_result <- fea_result$result %>% filter(term_size < 1000 | term_size > 10) $removed because gene set is small
  # select the top 50 pathways for plotting --- could do more, but difficult to see visually
  fea_result_filt <- fea_result$result %>% top_n(n = 50)
  return(fea_result_filt)
}

combined_fea_set <- function(genes, organism){
  upgenes <- as.list(genes %>% filter(diffexpressed == "UP"))
  up_fea_filt <- fea_set(genes = upgenes, organism = organism) %>% mutate(direction = "upregulated")
  downgenes <- as.list(genes %>% filter(diffexpressed == "DOWN"))
  down_fea_filt <- fea_set(genes = downgenes, organism = organism) %>% mutate(direction = "downregulated")
  combined_filt <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = combined_filt)
  return(fea_result)
}


## function-aggregate-matrix
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

## DmCast Function 
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
      if(nlevels(x)!=length(na.omit(unique(x))))
        x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })
  
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

# split violin plot function
library(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             } else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

# function targeting-Calc on panda regNet for gene and TF
targetingCalc <- function(regNetmatrix, variable_name, edge_weight_name, condition){
  #rearrange dataframe 
  regNetmatrix <- melt(regNetmatrix, varnames = c("TF", "gene"), value.name = "edge_weight_name")#melting dataframe
  print("datframe melted")
  regNetmatrix$edge_weight_name_pos <- ifelse(regNetmatrix$edge_weight_name < 0, 0, regNetmatrix$edge_weight_name) #replacing all negatives as a 0 and storing in new column
  print("subsetting only positive edge weights")
  regNetmatrix <- regNetmatrix[,c(1,2,4)]
  regNetmatrix
  
  #calculate gene targeting
  print("calculating gene targeting")
  Gene.targeting <- aggregate(.~gene, regNetmatrix[-1], sum) #removing TF column and calculating targeting for all edge weights and when edge weight is only positive
  #set column names based on condition
  print("renaming columns")
  if(condition == "het"){
    colnames(Gene.targeting) <- c("gene", "het_edge_weight_pos")
    print("plotting gene targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_GeneTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
      hist(Gene.targeting$het_edge_weight_pos)
      dev.off()
  } else {
    colnames(Gene.targeting) <- c("gene", "ctrl_edge_weight_pos")
    print("plotting gene targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_GeneTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(Gene.targeting$ctrl_edge_weight_pos)
    dev.off()
  }
  #reassign variable 
  print("assigning variable name to object")
  variable_name <- as.character(variable_name)
  assign(paste0(variable_name, "_gene_targeting_", condition), Gene.targeting, envir = .GlobalEnv)
  print("gene targeting calculation complete")
  
  #calculate TF targeting
  print("calculating TF targeting")
  TF.targeting <- aggregate(.~TF, regNetmatrix[-2], sum) #same as above but for TF instead of gene
  #set column names based on condition
  print("renaming columns")
  if(condition == "het"){
    colnames(TF.targeting) <- c("TF", "het_edge_weight_pos")
    print("plotting TF targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_TFTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(TF.targeting$het_edge_weight_pos)
    dev.off()
  } else {
    colnames(TF.targeting) <- c("TF", "ctrl_edge_weight_pos")
    print("plotting TF targeting score distribution")
    png(file = paste0(here("results/diff_targeting/"), variable_name, condition, "_TFTargetingScoresDist.png"),
        width = 1000,
        height = 1000)
    hist(TF.targeting$ctrl_edge_weight_pos)
    dev.off()
  }
  #reassign variable
  print("assigning variable name to object")
  variable_name <- as.character(variable_name)
  assign(paste0(variable_name, "_TF_targeting_", condition), TF.targeting, envir = .GlobalEnv)
  print("TF targeting calculation complete")
}

# function-targeting_heatmap; used in targeting
targeting_heatmap <- function(annotation_colors, data, meta_colname, plot_path, rowtitle, plot_title){
  #plotting all 
  ##grabbing metadata and annotations
  meta <- as.data.frame(colnames(data))
  colnames(meta) <- meta_colname
  rownames(meta) <- meta[,1]
  
  ##set heatmap annotations
  heat.anno = HeatmapAnnotation(df = meta, show_annotation_name = TRUE, col = annotation_colors)
  
  ##ensure column order matches annotation table
  data <- data[,rownames(meta), drop = FALSE]
  
  ##convert data to matrix
  mat <- as.matrix(data)
  
  ##plot heatmap 
  png(filename = plot_path,
      width = 1000,
      height = 1000)
  print(Heatmap(mat,
                col = colorRampPalette(brewer.pal(8,"Blues")) (25),
                heatmap_legend_param = list(title = "targeting score"),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_order = NULL,
                show_row_dend = TRUE,
                show_column_dend = TRUE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                use_raster = TRUE,
                raster_device = c("png"),
                bottom_annotation = NULL,
                top_annotation = heat.anno,
                column_title = plot_title, row_title = rowtitle, row_title_side = "right"))
  dev.off()
}