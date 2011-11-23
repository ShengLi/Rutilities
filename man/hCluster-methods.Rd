\docType{methods}
\name{hCluster}
\alias{hCluster}
\alias{hCluster,ANY-method}
\title{Hierarchical Cluster Analysis}
\usage{
  hCluster(x, dist.method = "correlation",
    hclust.method = "ward", plot = TRUE)
}
\arguments{
  \item{x}{a \code{matrix} object, each column is a sample,
  each row is a feature}

  \item{dist.method}{the distance measure to be used. This
  must be one of "\code{correlation}", "\code{euclidean}",
  "\code{maximum}", "\code{manhattan}", "\code{canberra}",
  "\code{binary}" or "\code{minkowski}". Any unambiguous
  substring can be given. (default:"\code{correlation}")}

  \item{hclust.method}{the agglomeration method to be used.
  This should be (an unambiguous abbreviation of) one of
  "\code{ward}", "\code{single}", "\code{complete}",
  "\code{average}", "\code{mcquitty}", "\code{median}" or
  "\code{centroid}". (default:"\code{ward}")}

  \item{plot}{a logical value indicating whether to plot
  hierarchical clustering. (default:TRUE)}
}
\description{
  Hierarchical Cluster Analysis
}

