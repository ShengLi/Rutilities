# Rutilities R package
# hCluster.R function

#---------------------------------------------------------------------------------------
# regular R functions to be used in S3 functions

# correlation function return dist object
# x matrix, each row is a sample.
# method only support "pearson"
# abs logical. If use absolute value for correlation to distance.
# diag logical. Inlcude diagnol value or not
# upper output upper panel or lower panel.
distCor=function(x, method="pearson", abs=TRUE, diag=FALSE, upper=FALSE)
{
  if (!is.na(pmatch(method, "pearson"))) 
    method <- "pearson"
  METHODS <- c("pearson")
  method <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid correlation method")
  if (method == -1) 
    stop("ambiguous correlation method")
  
  xcor = cor(t(x), method=METHODS[method])
  if(abs)
    xcor = 1-abs(xcor)
  else
    xcor = 1-xcor
  if(upper)
    d <- xcor[upper.tri(xcor,diag=diag)]
  else
    d <- xcor[lower.tri(xcor,diag=diag)]
  attr(d, "Size") <- nrow(x)
  attr(d, "Labels") <- dimnames(x)[[1L]]
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- METHODS[method]
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
} 

# end of regular functions to be used in S3 functions
#---------------------------------------------------------------------------------------

#' Hierarchical Cluster Analysis
#' 
#' @param x a \code{matrix} object, each column is a sample, each row is a feature
#' @param dist.method the distance measure to be used. This must be one of "\code{correlation}", "\code{euclidean}", "\code{maximum}", "\code{manhattan}", "\code{canberra}", "\code{binary}" or "\code{minkowski}". Any unambiguous substring can be given. (default:"\code{correlation}")
#' @param hclust.method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "\code{ward}", "\code{single}", "\code{complete}", "\code{average}", "\code{mcquitty}", "\code{median}" or "\code{centroid}". (default:"\code{ward}")
#' @param plot a logical value indicating whether to plot hierarchical clustering. (default:TRUE) 
#'
#' @export
#' @docType methods
#' @rdname hCluster-methods
#' @aliases hCluster,ANY-method
hCluster=function(x, dist.method="correlation", hclust.method="ward", plot=TRUE){
  DIST.METHODS <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski")
  dist.method <- pmatch(dist.method, DIST.METHODS)

  HCLUST.METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
        "median", "centroid")
  hclust.method <- pmatch(hclust.method, HCLUST.METHODS)
  if (is.na(hclust.method)) 
    stop("invalid clustering method")
  if (hclust.method == -1) 
    stop("ambiguous clustering method")

  if(DIST.METHODS[dist.method] == "correlation")
    d = distCor(t(x))
  else
    d=dist(scale(t(x)), method=DIST.METHODS[dist.method]);
  
  hc=hclust(d, HCLUST.METHODS[hclust.method]);
  
  if(plot){
    plclust(hc,hang=-1, main=paste("CpG dinucleotide methylation clustering\nDistance: ",
                                   DIST.METHODS[dist.method],sep=""), xlab = "Samples");
  }
  return(hc)
  }

#' @return a \code{tree} object of a hierarchical cluster analysis using a set of dissimilarities for the n objects being clustered.

