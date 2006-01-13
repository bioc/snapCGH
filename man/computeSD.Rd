\name{computeSD}
\alias{computeSD}
\alias{}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to estimate experimental variability of a sample}
\description{
  This functions estimate experimental variability of a given
  sample. This value can be used to rank samples in terms of the quality
  as well as to derive thresholds for declaring gained and lost clones.
}
\usage{
computeSD(segList, maxmadUse = 0.2, maxmedUse = 0.2,
               maxState = 3, maxStateChange = 10, minClone = 5,
               maxChrom = 22)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{segList}{
    Object of class \code{\link[snapCGH:SegList]{SegList}}.
  }
  \item{maxmadUse}{Maximum median absolute deviation allowed to
    controbute to the overall variability calculation.
  }
  \item{maxmedUse}{Maximum median value for a state allowed to
    contribute to the calculation.
  }
  \item{maxState}{Maximum number of the states on a given chromosome for
    the states from that chromosome to be allowed to enter noise
    variability calculation.
  }
  \item{maxStateChange}{Maximum number of changes from state to state on
    a given chromosome for that chromosome to enter  noise variability
    calculation.
  }
  \item{minClone}{Minimum number of clones in a state for clones in that
    sate to enter variability calculation.
  }
  \item{maxChrom}{Maxiumum chromosomal index (generally only autosomes
    are used for this calculation.
  }
}
\details{
  Median absolute deviation is estimated in all the states passing the
  criteria defined by the parameters of the function. Then median of all
  MADs on individual chromosomes as well as across all chromosomes is
  taken to estimate chromosomal  experimental variability and sample
  experimental variability.
}
\value{
  \item{madChrom}{Returns a matrix containing estimated variability for
    each chromosome for each sample.}
  \item{madGenome}{Returns a vector with estimate of experimental
    varibility for each sample.}
}
\seealso{
  \code{\link[snapCGH:SegList]{SegList}}
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Mike Smith}
\keyword{models}% at least one, from doc/KEYWORDS
