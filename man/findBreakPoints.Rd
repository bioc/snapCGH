\name{findBreakPoints}
\alias{findBreakPoints}
\title{Returns the start and end of segments.}
\description{
  Function to returns the start and end of segments when given a SegList
  and an array.  Currently only used within the
  \code{\link{plotSegmentedGenome}} function.
}

\usage{
findBreakPoints(seg, array)
}

\arguments{
  \item{seg}{An object of class \code{\link[snapCGH:SegList]{SegList}}.}
  \item{array}{Numeric value corresponding to a column in seg.}
}
 
\author{Mike Smith}

\keyword{methods}
