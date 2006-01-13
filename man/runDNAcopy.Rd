\name{runDNAcopy}
\alias{runDNAcopy}
\title{Results of segmenting an MAList data object using the DNAcopy library}
\description{
  The results of segmenting data from copy number array experiments from
  programs such as circular binary segmentation (CBS).  This function requires
  the library \code{DNAcopy} to be loaded.
}
\usage{
runDNAcopy(MA)
}
\arguments{
  \item{MA}{An object of class \code{\link[limma]{MAList}}}
}
\value{
  The function returns an object of class \code{\link[snapCGH:SegList]{SegList}}
}
\seealso{
  \code{\link[limma]{MAList}}
  \code{\link{runHomHMM}}
  \code{\link{runGLAD}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\keyword{methods}