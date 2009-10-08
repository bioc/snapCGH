\name{imputeMissingValues}
\alias{imputeMissingValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Imputing log2 ratios}
\description{
  Imputing log2 ratios
}
\usage{
imputeMissingValues(seg, chrominfo = chrominfo.Mb, maxChrom =
23, smooth = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{Object of class \code{\link[snapCGH:SegList]{SegList}}}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{maxChrom}{Highest chromosome to impute}
  \item{smooth}{smoothing parameter for the lowess procedure}
}
\details{
  There are two main reasons to impute data.
  One is that given that imputation is reasonable,
  one can increase the analytical power and improve results.
  Another, more practical, is that at the moment many widely used
  functions in R do not support missing values.
  While procedures such as kNN imputations is widely used for gene
  expression data, it is more powerful to take advantage of the
  genomic structure of the array CGH data and use a smoother.
  Note that we perform only one pass of smoothing.
  If there still remain missing values, they are imputed by the
  median on the chromosome or chromosomal arm where applicable.
}
\value{
  Computes and returns the imputed log2 ratio matrix of the aCGH object.
}
\seealso{
  \code{\link[snapCGH:SegList]{SegList}}
}


\keyword{models}% at least one, from doc/KEYWORDS
