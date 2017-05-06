###########################################################################
# dotplot modified to order and color by any column of gseaResult and enrichResult
# order high to low, top to bottom on the graph
##' dotplot for gseaResult
##'
##'
##' @rdname dotplot-methods
##' @aliases dotplot,gseaResult,ANY-method
##' @exportMethod dotplot
##' @author Guangchuang Yu
setMethod("dotplot", signature(object="gseaResult"),
          function(object, x="geneRatio", colorBy="p.adjust", orderBy = "GeneRatio", showCategory=10, split=NULL, font.size=12, title="") {
            dotplot_internal(object = object, x, colorBy, orderBy, showCategory, split, font.size, title)
          }
)

##' dotplot for enrichResult
##'
##'
##' @rdname dotplot-methods
##' @aliases dotplot,enrichResult,ANY-method
##' @param object an instance of enrichResult
##' @param x variable for x axis
##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of category
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @exportMethod dotplot
##' @author Guangchuang Yu
setMethod("dotplot", signature(object="enrichResult"),
          function(object, x="geneRatio", colorBy="p.adjust", orderBy = "GeneRatio", showCategory=10, split=NULL, font.size=12, title="") {
            dotplot_internal(object, x, colorBy, orderBy, showCategory, split, font.size, title)
          }
)

##' @importFrom ggplot2 fortify
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
dotplot_internal <- function(object, x="geneRatio", colorBy="p.adjust", orderBy = "GeneRatio", showCategory=10, split=NULL, font.size=12, title="") {
  colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue", "enrichmentScore"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    size <- "GeneRatio"
  } else {
    stop("x should be geneRatio or count...")
  }
  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  if(colorBy == "enrichmentScore") {lows = "blue"; highs = "red"}
  if(colorBy != "enrichmentScore") {lows = "red"; highs = "blue"}
  
  idx <- order(df[,orderBy], decreasing = FALSE)
  df$Description <- factor(df$Description, levels=df$Description[idx])
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() + scale_color_gradient(low=lows, high=highs) +
    ggtitle(title) + theme_dose(font.size) + 
    ylab(ifelse(orderBy == "p.adjust",">> adjusted p-value increasing >>", ""))+
    xlab(ifelse(x == "GeneRatio" & class(object) == "gseaResult",
                "the fraction of proteins from a gene set which are over- or underrepresented", 
                ifelse(x == "GeneRatio" & class(object) == "enrichResult",
                       "the fraction of proteins from a gene set in the analysed set", 
                       x)))
}

##' convert enrichResult object for ggplot2
##'
##'
##' @title fortify
##' @param model enrichResult object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param by one of Count and GeneRatio
##' @param order logical
##' @param drop logical
##' @param split separate result by 'split' variable
##' @param ... additional parameter
##' @importFrom ggplot2 fortify
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
  fortify.internal(model, data, showCategory, by, order, drop, split, ...) 
}

##' @method fortify enrichResult
##' @export
fortify.gseaResult <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
  fortify.internal(model, data, showCategory, by, order, drop, split, ...) 
}


fortify.internal <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
  res <- as.data.frame(model)
  if (inherits(model, "gseaResult")) {
    res$Count <- str_count(res$core_enrichment, "/") + 1
    res$.sign <- "activated"
    res$.sign[res$NES < 0] <- "suppressed"
  }
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  if (inherits(model, "gseaResult")) {
    res$GeneRatio <- res$Count / res$setSize
  } else if (inherits(model, "enrichResult")) {
    res$GeneRatio <- parse_ratio(res$GeneRatio)
  }
  
  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }
  
  topN <- function(res, showCategory) {
    if ( is.numeric(showCategory) ) {
      if ( showCategory <= nrow(res) ) {
        res <- res[1:showCategory,]
      }
    } else { ## selected categories
      res <- res[res$ID %in% showCategory,]
    }
    return(res)
  }
  
  if (is.null(split)) {
    res <- topN(res, showCategory)
  } else {
    lres <- split(res, as.character(res[, split]))
    lres <- lapply(lres, topN, showCategory = showCategory)
    res <- do.call('rbind', lres)
  }
  
  res$Description <- factor(res$Description,
                            levels=rev(res$Description))
  
  return(res)
}

str_count <- function(string, pattern="") {
  sapply(string, str_count_item, pattern=pattern)
}

str_count_item <- function(string, pattern = "") {
  length(gregexpr(pattern, string)[[1]])
}

parse_ratio <- function(ratio) {
  gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
  gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
  return(gsize/gcsize)
}
###########################################################################
