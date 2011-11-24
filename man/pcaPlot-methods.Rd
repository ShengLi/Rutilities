\docType{methods}
\name{pcaPlot}
\alias{"pcaPlot,ANY-methods"}
\alias{pcaPlot}
\title{Principal Components Analysis for scatter plot and screeplot}
\usage{
  pcaPlot(x, cor = TRUE, screeplot = FALSE)
}
\arguments{
  \item{x}{a \code{matrix} object, each column is a sample,
  each row is a feature}

  \item{cor}{a logical value indicating whether the
  calculation should use the correlation matrix or the
  covariance matrix. (default: TRUE)}

  \item{screeplot}{a logical value indicating whether to
  plot the variances against the number of the principal
  component. (default: FALSE)}
}
\value{
  The form of the value returned by \code{pcaPlot} is the
  summary of principal component analysis by
  \code{princomp}.
}
\description{
  Principal Components Analysis for scatter plot and
  screeplot
}
\examples{
x = rnorm(100, 100, 5)
y = cbind(x, 2*x+rnorm(100), 1:100, 1:100 + rnorm(100,0.5,2))
colnames(y) = c('test1','test2','ctrl1','ctrl2')
pcaPlot(y,screeplot=TRUE)
pcaPlot(y)
}

