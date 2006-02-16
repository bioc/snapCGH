\name{plotGainLoss}
\alias{plotGainLoss}
\alias{changeProp.func}
\alias{changeProp.overall.func}
\alias{table.bac.func}
\alias{lengthGain.na}
\alias{propGain.na}
\alias{lengthLoss.na}
\alias{propLoss.na}
\alias{prop.na}
\alias{create.resT}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{frequency plots and significance analysis}
\description{
  The main application of this function is to plot the frequency of
  changes.
}
\usage{
plotGainLoss(segList, resT = NULL, pheno = rep(1,
                 ncol(log2ratios(segList))), chrominfo =
                 chrominfo.Mb, X = TRUE, Y = FALSE, rsp.uniq =
                 unique(pheno), all = length(rsp.uniq) == 1 &&
                 is.null(resT), titles = if (all) "All Samples" else
                 rsp.uniq, cutplot = 0, thres = 0.25, factor = 2.5, ylm
                 = c(-1, 1), p.thres = c(0.01, 0.05, 0.1), numaut = 22,
                 onepage = TRUE, colored = TRUE)
}
\arguments{
  \item{segList}{Object of class \code{SegList}}
  \item{resT}{Data frame having the same structure as the result of
    applying \code{\link[aCGH:mt.maxT]{mt.maxT}} or \code{\link[aCGH:mt.maxT]{mt.minP}} functions
    from Bioconductor's \code{multtest} package for multiple testing.
    The result is a data frame including the following 4 components:
    'index', 'teststat', 'rawp' and 'adjp'.
  }
  \item{pheno}{phenotype to compare.}
  \item{chrominfo}{
    Chromosomal information. Defaults to
    \code{\link{chrominfo.Mb}}
  }
  \item{X}{Include X chromosome? Defaults to yes.}
  \item{Y}{Include Y chromosome? Defaults to no.}
  \item{rsp.uniq}{\code{rsp.uniq} specified the codes for the groups of
    interest. Default is the unique levels of the phenotype. Not used
    when \code{all} is T.
  }
  \item{all}{
    \code{all} specifies whether samples should be analyzed by subgroups
    (T) or together (F).
  }
   \item{titles}{
    \code{titles} names of the groups to be used. Default is the unique
    levels of the \code{pheno}.
  }
  \item{cutplot}{only clones with at least \code{cutplot} frequency of
    gain and loss are plotted.
  }
  \item{thres}{\code{thres} is either a vector providing unique
    threshold for each sample or a vector of the same length as number
    of samples (columns in \code{data}) providing sample-specific
    threshold. Used for plotting
    the gain/loss frequency data as well as for clone screening and for
    significance analysis when \code{threshold} is TRUE.Defaults to 0.25
  }
  \item{factor}{\code{factor} specifies the number by which experimental
    variability should be multiplied.  Defaults to 2.5}
  \item{ylm}{\code{ylm} vertical limits for the plot}
  \item{p.thres}{
    \code{p.thres} vector of p-value ciut-off to be plotted. computed
    conservatively as the threshold corresponding to a given adjusted
    p-value.
  }
  \item{numaut}{\code{numaut} number of the autosomes}
  \item{onepage}{
    \code{onepage} whether all plots are to be plotted on one page or
    different pages. When more than 2 groups are compared, we recommend multiple pages.
  }
  \item{colored}{Is plotting in color or not? Default is TRUE.}
}
\examples{

}
\keyword{htest}
\keyword{hplot}
