## clusterProfiler simplify methods modified to calculate kappa score similarity 
## between terms and to use data.table for faster table manipulations
##' simplify output from enrichGO by removing redundancy of enriched GO terms
##'
##'
##' @name simplify
##' @docType methods
##' @rdname simplify-methods
##' @title simplify method
##' @param x output of enrichGO
##' @param cutoff similarity cutoff
##' @param by feature to select representative term, selected by 'select_fun' function
##' @param select_fun function to select feature passed by 'by' parameter
##' @param measure method to measure similarity
##' @param semData GOSemSimDATA object
##' @return updated enrichResult object
##' @exportMethod simplify
##' @references issue #28
##' \url{https://github.com/GuangchuangYu/clusterProfiler/issues/28}
##' @aliases simplify,enrichResult-method
##' @author Guangchuang Yu
setMethod("simplify", signature(x="enrichResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL, use_data_table = F, use_bioc_annotationdbi = T) {
            library(magrittr)
            if (!x@ontology %in% c("BP", "MF", "CC"))
              stop("simplify only applied to output from enrichGO...")
            
            
            x@result %<>% simplify_internal(., cutoff, by, select_fun,
                                            measure, x@ontology, semData, use_data_table, use_bioc_annotationdbi, x@keytype)
            
            return(x)
          }
)

##' @author Guangchuang Yu & Vitalii Kleshchevnikov
setMethod("simplify", signature(x="gseaResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL, use_data_table = F, use_bioc_annotationdbi = T) {
            library(magrittr)
            if (!x@setType %in% c("BP", "MF", "CC"))
              stop("simplify only applied to output from enrichGO...")
            
            
            x@result %<>% simplify_internal(., cutoff, by, select_fun,
                                            measure, x@setType, semData, use_data_table, use_bioc_annotationdbi, x@keytype)
            
            return(x)
          }
)

##' @importFrom GOSemSim mgoSim
##' @importFrom GOSemSim godata
##' @importFrom tidyr gather
simplify_internal <- function(res = x@result, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel", ontology = x@ontology, semData, use_data_table, use_bioc_annotationdbi, keytype = x@keytype) {
  if (missing(semData) || is.null(semData)) {
    if (measure == "Wang") {
      semData <- GOSemSim::godata(ont = ontology)
    } else {
      if(measure != "kappa") stop("godata should be provided for IC-based methods...")
    }
  } else {
    if (ontology != semData@ont) {
      msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
      stop(msg)
    }
  }
  
  if (measure == "Wang"){
    sim <- GOSemSim::mgoSim(res$ID, res$ID,
                            semData = semData,
                            measure=measure,
                            combine=NULL) 
  }
  
  if(measure == "kappa"){
    library(data.table)
    library(org.Hs.eg.db)
    # find GO term to protein mapping
    semData <- clusterProfiler:::get_GO_data(OrgDb = 'org.Hs.eg.db', ont = ontology, keytype = keytype)
    # measure distance between selected terms (important to get all GO to protein mappings)
    semData = data.table(names(unlist(semData$EXTID2PATHID)),unlist(semData$EXTID2PATHID))
    colnames(semData) = c(keytype,"GO")
    semData = unique(semData)
    sim = categ_dist(mapping_table = semData, terms_to_compare = res$ID)$similarity_matrix
  }
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  if(use_data_table == F){
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- tidyr::gather(sim.df, go2, similarity, -go1) 
    sim.df <- sim.df[!is.na(sim.df$similarity),]
    
    ## feature 'by' is attached to 'go1'
    sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
    sim.df$go2 <- as.character(sim.df$go2)
    
    ID <- res$ID
    
    GO_to_remove <- character()
    for (i in seq_along(ID)) {
      ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
      ## if length(ii) == 1, then go1 == go2
      if (length(ii) < 2)
        next
      
      sim_subset <- sim.df[ii,]
      
      jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
      
      ## sim.df <- sim.df[-ii[-jj]]
      GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }
  }
  
  if(use_data_table == T){
    library(data.table)
    sim.df <- as.data.table(sim, keep.rownames = "go1")
    sim.df <- melt(sim.df,id.vars = "go1",variable.name = "go2",
                   value.name = "similarity", variable.factor =F) 
    sim.df <- sim.df[!is.na(similarity),]
    
    ## feature 'by' is attached to 'go1'
    res2 = as.data.table(res)
    sim.df <- merge(sim.df, res2[, c("ID",by), with = F], by.x="go1", by.y="ID")
    sim.df$go2 <- as.character(sim.df$go2)
    
    ID <- res2$ID
    
    GO_to_remove <- character()
    for (i in seq_along(ID)) {
      ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
      ## if length(ii) == 1, then go1 == go2
      if (length(ii) < 2)
        next
      
      sim_subset <- sim.df[ii,]
      
      jj <- which(sim_subset[, by, with =F] == select_fun(sim_subset[, by, with =F]))
      
      ## sim.df <- sim.df[-ii[-jj]]
      GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }
  }
  
  res[!res$ID %in% GO_to_remove, ]
}
