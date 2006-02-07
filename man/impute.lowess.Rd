\name{imputeMissingValues}
\alias{imputeMissingValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Imputing log2 ratios}
\description{
  Imputing log2 ratios
}
\usage{
imputeMissingValues(MA, chrominfo = chrominfo.Mb, maxChrom =
23, smooth = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MA}{Object of class \code{\link[limma]{MAList}}}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{maxChrom}{Highest chromosome to impute}
  \item{smooth}{smoothing parameter for the lowess procedure}
}
\details{
  There are two main reasons to impute data.
  One is that given that imputation is reasonable,
  one can increase the analytical power and improve results.
  Another, more practical, is that at the moment many widely used
  fuctions in R do not support missing values.
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
  \code{\link[limma]{MAList}}
}
\examples{

%MA1 <- process.MAList(MA1,chrom.below.threshold=1,chrom.remove.threshold=24,prop.missing=1,method.of.averaging = mean)

%## Imputing the log2 ratios 

%MA2 <- imputeMissingValues(MA1,maxChrom=24)

}
\keyword{models}% at least one, from doc/KEYWORDS
