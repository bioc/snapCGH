\name{cbind}
\alias{cbind.SegList}
\title{Combine SegList Objects}
\description{
Combine a series of \code{SegList} objects. 
}
\usage{
\method{cbind}{SegList}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{SegList} objects}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
\code{cbind} combines data objects assuming the same gene lists but different arrays.

For \code{cbind}, the matrices of expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
There is no \code{rbind}
function for the \code{\link[snapCGH:SegList]{SegList}} object.  This is
  because it would require \code{SegLists} with mutually exclusive
  chromosomes or would result in combining multiple different segmentations
  for the same chromosome, which is pointless.  If \code{rbind} is
  required perform it on an \code{\link[limma:malist]{MAList}} and then
    segment it.
}
\value{
An \code{\link[snapCGH:SegList]{SegList}} object holding data from all the arrays and all genes from the individual objects.
}
\author{Gordon Smyth, modified by Mike Smith for SegList object}
\seealso{
  \code{\link[base]{cbind}} in the base package.
 }
\keyword{manip}
