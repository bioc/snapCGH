\name{runPicard}
\alias{runPicard}
\title{Results of segmenting an MAList data object using the Picard et
  al algorithm found in the tilingArray library}
\description{
  Wrapper calling the Picard et al segmentation algorithm on an MAList
  object.  This function requires the library \code{DNAcopy} to be loaded.
}
\usage{
  runPicard(input, maxSeg = 5, maxk = 200, criteria = "BIC")
}
\arguments{
  \item{input}{An object of class \code{\link[limma]{MAList}} or \code{\link[snapCGH:SegList]{SegList}}}
  \item{maxSeg}{integer of length 1, maximum number of segments (= 1 +
    maximum number of change points)}
  \item{maxk}{integer of length 1, maximum length of a single segment}
  \item{criteria}{Criteria for model selection.  Options are "none",
    "AIC" and "BIC" (default)}
}
\value{
  The function returns an object of class \code{\link[snapCGH:SegList]{SegList}}
}
\seealso{
  \code{\link[tilingArray]{segment}}
  \code{\link[limma]{MAList}}
  \code{\link{runHomHMM}}
  \code{\link{runGLAD}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\author{Mike Smith}
\keyword{methods}
