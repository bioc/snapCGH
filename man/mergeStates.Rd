\name{mergeStates}
\alias{mergeStates}
\alias{MergeLevels.new}
\alias{MergeLevels.old}
\alias{combine.func}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Funtion to merge states based on their state means}
\description{
  \code{mergeStates} takes the output of a segementation algorithm in
  the form of a \code{\link[snapCGH:SegList]{SegList}} and iteratively merges the
  states with means closer than a supplied threshold.
}
\usage{
mergeStates(segList, MergeType = 1, pv.thres=0.0001, ansari.sign=0.01, minDiff = 0.25)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{segList}{
    Object of class \code{\link[snapCGH:SegList]{SegList}}.
  }
  \item{MergeType}{
    Select either 1 or 2.  1 uses a new merging algorithm developed by
    Hanni Willenbrock and Jane Fridlyand.
  }
  \item{pv.thres}{Significance threshold for Wilcoxon test for
    level merging.  Used when MergeType = 1.}
  \item{ansari.sign}{Significance threshold for Ansari-Bradley test. Used when MergeType = 1.}
  \item{minDiff}{The states whose predicted values are less than minDiff
    apart are merged into one state and all the predicited values are
    recomputed.  Used when MergeType = 2.   
  }
}
\details{
  This function is intended to reduce effect of the possible small
  magnitude technological artifacts on the structure determination.
}
\value{
  A \code{\link[snapCGH:SegList]{SegList}} object is returned with the merged states
  stored in the \code{pred} list.
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand}
\seealso{\code{\link[snapCGH:SegList]{SegList}},
  \code{\link{fitHMM}},
  \code{\link{run.GLAD}},
  \code{\link{run.DNAcopy}}
  }

\keyword{models}% at least one, from doc/KEYWORDS
