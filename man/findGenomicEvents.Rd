\name{findGenomicEvents}
\alias{findGenomicEvents}
\title{Finds the genomic events associated with each of the array CGH samples}
\description{
  Finds the genomic events associated with each of the array CGH
  samples. Events include whole chromosomal gains and losses,
  aberrations, transitions, amplifications and their respective counts
  and sizes. The hmm states has to be computed before using this
  function.
}
\usage{
findGenomicEvents(segList, maxChrom = 22, factor = 5, maxClones = 1,
                    maxLen = 1000, absValSingle = 1, absValRegion = 1,
                    diffVal1 = 1, diffVal2 = 0.5, maxSize = 10000,
                    pChrom.min = 0.9, medChrom.min = 0.1,
                    maxmadUse = 0.2, maxmedUse = 0.2, 
                    maxState = 3, maxStateChange = 10, minClone = 5)
}
\arguments{
  \item{segList}{
    Object of class \code{\link[snapCGH:SegList]{SegList}}.
  }
  \item{maxChrom}{Highest chromosomal number to find events.}
  \item{factor}{Determines outliers. See \code{\link{findOutliers}}.}
  \item{maxClones}{Determines aberrations. See \code{\link{findAberrations}}.}
  \item{maxLen}{Determines aberrations. See \code{\link{findAberrations}}.}
  \item{absValSingle}{Determines amplifications. See \code{\link{findAmplifications}}.}
  \item{absValRegion}{Determines amplifications. See \code{\link{findAmplifications}}.}
  \item{diffVal1}{Determines amplifications. See \code{\link{findAmplifications}}.}
  \item{diffVal2}{Determines amplifications. See \code{\link{findAmplifications}}.}
  \item{maxSize}{Determines amplifications. See \code{\link{findAmplifications}}.}
  \item{pChrom.min}{
    Determines whole chromosomal gains and losses. Chromosome should
    contain no transitions, have its absolute median equal or greater
    than \code{medChrom.min} and at least \code{medChrom.min} has to be
    greater or less than 0.
  }
  \item{medChrom.min}{
    Determines whole chromosomal gains and losses. Chromosome should
    contain no transitions, have its absolute median equal or greater
    than \code{medChrom.min} and at least \code{medChrom.min} has to be
    greater or less than 0.
  }
  \item{maxmadUse}{See: \code{\link{computeSD}}}
  \item{maxmedUse}{See: \code{\link{computeSD}}}
  \item{maxState}{See: \code{\link{computeSD}}}
  \item{maxStateChange}{See: \code{\link{computeSD}}}
  \item{minClone}{See: \code{\link{computeSD}}}
}
\details{
  The default parameters generally work. Threshold for merging may be
  changed depending on the expected normal cell contamination and/or
  expected magnitude of the changes. AIC model generally works,
  however, may need to be readjusted depending on how liberal or
  conservative one wants to be in finding genomic events. We recommend
  BIC criterion with delta = 1 for noisier data.
}

\value{
  \item{num.transitions}{
    matrix of dimensions maxChrom by number of samples. It contains
    number of transitions that were recorded on a given chromosome for a
    given sample.
  }
  \item{num.amplifications}{
    matrix of dimensions maxChrom by number of samples It contains
    number of amplifications that were recorded on a given chromosome
    for a given sample.
  }
  \item{num.aberrations}{
    matrix of dimensions maxChrom by number of samples. It contains
    number of focal aberrations that were recorded on a given chromosome
    for a given sample.
  }
  \item{num.outliers}{
    matrix of dimensions maxChrom by number of samples.  It contains
    number of outliers that were recorded on a given chromosome for a
    given sample.
  }
  \item{num.transitions.binary}{
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more transitions were recorded on a
    given chromosome for a given sample.
  }
  \item{num.amplifications.binary}{
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more amplifications were recorded on a
    given chromosome for a given sample.
  }
  \item{num.aberrations.binary}{
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more focal aberrations were recorded on
    a given chromosome for a given sample.
  }
  \item{num.outliers.binary}{
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more outliers were recorded on a given
    chromosome for a given sample.
  }
  \item{whole.chrom.gain.loss}{
    matrix of dimensions maxChrom by number of samples. Positive entry
    indicates that a given chromosome was gained in a given sample,
    negative entry indicates that a given chromosome was lost in a given
    sample, 0 entry is normal chromosome and NA marks chromosomes with
    one or more transition.
  }
  \item{size.amplicons}{
    matrix of dimensions maxChrom by number of samples. Reports size of
    a given chromosome that is amplified (kb units) in a given sample.
  }
  \item{num.amplicons}{
    matrix of dimensions maxChrom by number of samples. Reports number
    of disjoint amplicons on a given chromosome for a given sample.
  }
  \item{outliers}{
    list containing 3 matrices of dimensions number of clones by number
    of samples. See \code{\link{findOutliers}}.
  }
  \item{aberrations}{
    list containing a matrix of dimensions number of clones by number
    of samples. See \code{\link{findAberrations}}.
  }
  \item{transitions}{
    list containing 2 matrices of dimensions number of clones by number
    of samples. See \code{\link{findTranslocations}}.
  }
  \item{amplifications}{
    list containing a matrix of dimensions number of clones by number
    of samples. See \code{\link{findAmplifications}}.
  }
}
\seealso{
  \code{\link{runHomHMM}}
  \code{\link{mergeStates}}
  \code{\link{findAberrations}}
  \code{\link{findTranslocations}}
  \code{\link{findAmplifications}}
  \code{\link{findOutliers}}
}
\keyword{models}
