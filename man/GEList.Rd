\name{GEList-class}
\docType{class}
\alias{GEList-class}
\title{Genomic Events - class}

\description{
  A list based class for storing the location of points of interest
  within the data set following the segmentation step.  Events include whole chromosomal gains and losses,
  aberrations, transitions, amplifications and their respective counts
  and sizes.  They are
  generally created by running the function \code{\link{findGenomicEvents}}.
}

\section{Slots/List Components}{
  Objects should contain the following list components:
  \tabular{ll}{
    \code{num.transitions}:\tab matrix of dimensions maxChrom by number of samples. It contains number of transitions that were recorded on a given chromosome for a given sample.\cr
     \code{num.amplifications}:\tab
    matrix of dimensions maxChrom by number of samples It contains
    number of amplifications that were recorded on a given chromosome
    for a given sample.\cr
  \code{num.aberrations}:\tab
    matrix of dimensions maxChrom by number of samples. It contains
    number of focal aberrations that were recorded on a given chromosome
    for a given sample.\cr
\code{num.outliers}:\tab
    matrix of dimensions maxChrom by number of samples.  It contains
    number of outliers that were recorded on a given chromosome for a
    given sample.\cr
  \code{num.transitions.binary}:\tab
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more transitions were recorded on a
    given chromosome for a given sample.\cr
  \code{num.amplifications.binary}:\tab
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more amplifications were recorded on a
    given chromosome for a given sample.\cr
  \code{num.aberrations.binary}:\tab
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more focal aberrations were recorded on
    a given chromosome for a given sample.\cr
  \code{num.outliers.binary}:\tab
    binary matrix of dimensions maxChrom by number of samples. Non-zero
    entry indicates whether 1 or more outliers were recorded on a given
    chromosome for a given sample.\cr
  \code{whole.chrom.gain.loss}:\tab
    matrix of dimensions maxChrom by number of samples. Positive entry
    indicates that a given chromosome was gained in a given sample,
    negative entry indicates that a given chromosome was lost in a given
    sample, 0 entry is normal chromosome and NA marks chromosomes with
    one or more transition.\cr
  \code{size.amplicons}:\tab
    matrix of dimensions maxChrom by number of samples. Reports size of
    a given chromosome that is amplified (kb units) in a given sample.\cr
  \code{num.amplicons}:\tab
    matrix of dimensions maxChrom by number of samples. Reports number
    of disjoint amplicons on a given chromosome for a given sample.\cr
  \code{outliers}:\tab
    list containing 3 matrices of dimensions number of clones by number
    of samples. See \code{\link{findOutliers.func}}.\cr
  \code{aberrations}:\tab
    list containing a matrix of dimensions number of clones by number
    of samples. See \code{\link{findAber.func}}.\cr
  \code{transitions}:\tab
    list containing 2 matrices of dimensions number of clones by number
    of samples. See \code{\link{findTrans.func}}.\cr
  \code{amplifications}:\tab
    list containing a matrix of dimensions number of clones by number
    of samples. See \code{\link{findAmplif.func}}.\cr
  }
}

\seealso{
  \code{\link[snapCGH:SegList]{SegList}}
  \code{\link{plotSegmentationStates}}
}

\author{Mike Smith}

 \keyword{classes}
 \keyword{data}