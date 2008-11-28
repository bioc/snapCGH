\name{compareSegmentations}
\alias{compareSegmentations}
\alias{compareBreakPoints}

\title{Function for comparing segmentation methods to a known truth}
\description{
  This function takes a SegList and compares the breakpoints indicated
  in other SegLists with this original one.
}
\usage{
compareSegmentations(TrueSeg,offset = 0,...)
}

\arguments{
  \item{TrueSeg}{An object of class
    \code{\link[snapCGH:SegList]{SegList}} which is scored against.
    Normally the output from \code{\link{simulateData}}.}
  \item{offset}{Integer value between 0 and 2 specifying how close (in
    number of clones) to a true breakpoint the segmentation method must be before it is scored.}
  \item{...}{One or more objects of class
    \code{\link[snapCGH:SegList]{SegList}}.  These are compared to TrueSeg.}
}
\value{
  The method returns a list containing two 
  matrices.  The first of these, \$TPR, contains the true positive rate, 
  whilst the second, \$FDR, holds the false discovery rate.  Both of these
  matrices are arranged such that a row represents a segmentation method and 
  each column is an array.
}

\author{John Marioni and Mike Smith}

\keyword{manip}
