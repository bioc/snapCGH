\name{findAmplifications}
\alias{findAmplifications}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to determine high level amplifications}
\description{
  This function identifies high level amplifications by considering the
  height, the width of an amplicon relative to the urrounding
  clones. Only narrow peaks much higher than its neigbors are considered
  as high level amplifications.
}
\usage{
findAmplifications(segList, absValSingle = 1, absValRegion = 1.5, diffVal1 = 1,
diffVal2 = 0.5, maxSize = 10000, translen.matr = res3$translen.matrix,
trans.matr = res3$trans.matr, aber = res2$aber, outliers = res1$outlier,
pred = res1$pred.out, pred.obs = res1$pred.obs.out)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{segList}{An object of class \code{\link[snapCGH:SegList]{SegList}}}
  \item{absValSingle}{A clone is declared to be an amplification if it
    is a focal aberration or an outlier and its value exceeds
    absValSingle}
  \item{absValRegion}{A clone is an amplification if if a clone belong
    to a region with width less than maxSize and  observed value for a
    clones is greater than absValRegion}
  \item{diffVal1}{Clone is an amplification if it is an aberration and
    greater by diffVal1 than max of the two surrounding stretches}
  \item{diffVal2}{Clone is an amplification if it is an outlier, its
    observed values is greater by diffVal2 than max of the two
    surrounding stretches}
  \item{maxSize}{The clones  may not be declared as amplifications if
    they belong to the states with spanning more than maxSize}
  \item{translen.matr}{State length matrix. The output of the
    \code{\link{findTranslocations}}}
  \item{trans.matr}{Transition matrix. The output of the
    \code{\link{findTranslocations}}}
  \item{aber}{Aberration matrix. The output of the
    \code{\link{findAberrations}}}
  \item{outliers}{Outliers matrix. The output of the
    \code{\link{findOutliers}}}
  \item{pred}{Predicted values matrix. The output of the
    \code{\link{findOutliers}}}
  \item{pred.obs}{Predicted values matrix with observed values assigned
    to the outliers. The output of the \code{\link{findOutliers}}}
  }
\details{
  Note that all the distances are in Megabases and all the heights are
  on log2ratio scale.
}
\value{
  \item{amplif.matrix}{Binary matrix with a row for each clone and
    column for each sample. "1" indicates amplification}
  ...
}
\seealso{
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH  data, Fridlyand et.al., \emph{JMVA}, 2004}
\author{Jane Fridlyand}
\keyword{models}% at least one, from doc/KEYWORDS
