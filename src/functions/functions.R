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

# Function-define-pK
define_pK <- function(sample){
  sweep_res_list <- paramSweep_v3(sample, PCs = 1:30, sct = FALSE)
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

# Function-remove-doublets
doublet_removal <- function(sample, expected, pK){
  # define the expected number of doublets in nuclei.
  nExp <- round(ncol(sample) * expected)  # expect 10.2% doublets
  sample <- doubletFinder_v3(sample, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:30)
  
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
