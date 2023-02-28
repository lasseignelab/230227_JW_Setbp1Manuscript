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