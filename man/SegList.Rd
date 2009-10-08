\name{SegList-class}
\docType{class}
\alias{SegList-class}
\alias{print.SegList}
\alias{show,SegList-method}
\alias{[.SegList}
\title{Segmentation States - class}

\description{
A list based class for storing the results of a segmentation algorithm.
They are generally created by running one of the following functions
\code{\link{runHomHMM}}, \code{\link{runGLAD}} or \code{\link{runDNAcopy}}
on an \code{\link[limma:malist]{MAList} object}.
}

\section{Slots/List Components}{
Objects should contain the following list components:
\tabular{ll}{
  \code{pred}:\tab Predicted value of the state.\cr
  \code{disp}:\tab Dispersion.\cr
  \code{obs}:\tab Observed value.\cr
  \code{state}:\tab Numeric value.\cr
  \code{nstates.hmm}: \tab The number of states per chromosome.  Each row represents a chromosome and each column is an array.\cr
  \code{genes}:\tab data.frame that contains the chromosome and position on the chromosome for each clone.  Used for plotting functions.\cr
}
Optional:\cr
\tabular{ll}{
  \code{rpred}:\tab Smoothed value for the clone.\cr
  \code{prob}:\tab Probability of the clone being in the assigned state.\cr
}
}

\section{Methods}{
\code{SegLists} can be subsetted and combined.  They also return dimensions
so functions such as \code{\link[base]{dim}}, \code{\link[base]{nrow}} and
\code{\link[base:nrow]{ncol}} are also defined.  \code{SegList} inherits the
\code{\link[methods]{show}} method from the Limma class
\code{\link[limma]{LargeDataObject}}.  This means that the \code{SegList}
will print in a relatively compact way.
}
\author{Mike Smith}

\keyword{classes}
\keyword{data}
