\name{fit.model}
\alias{fit.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting a heterogneous HMM to the log2 ratios on a particular chromosome.}
\description{
  This function (which is a workhorse of process.data) fits five heterogeneous HMMs to the log2 ratios on a particular chromosome. It then uses either the
  AIC or BIC to determine which of the five models is optimal before using a scaled version of the Viterbi algorithm to assign clones to states with the
  same underlying copy number.
}
\usage{
fit.model(sample, chrom, dat, datainfo = clones.info, covariates, vr =
0.01, iterlim = iterlim, criteria = 1, delta = 1, nlists = 1,var.fixed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{ If there are multiple samples, the number of the sample to be segmented }
  \item{chrom}{ The chromosome on which the segmentation is to be carried out on }
  \item{dat}{ The log2 ratios obtained from the clones located on that chromosome }
  \item{datainfo}{ A dataframe containing information about the clones on that chromosome (name, chromosome and location (in Mbs)) }
  \item{covariates}{ A matrix contianing the covariate information for the clones located on the chromosome to be segmented. It should have length one less
                     than the number of clones as covariate information is not used when segmenting the first clone on the chromosome. }
  \item{vr}{ SUPERFLUOUS - TO BE REMOVED }
  \item{iterlim}{ The maximum number of times the iterative procedure should be carried out }
  \item{criteria}{Select 1 if you want to use the AIC to distinguish
    between states with different numbers of underlying states or 2 if
    you wish to use the BIC. }
  \item{delta}{ A parameter to be set if you want to use the BIC }
  \item{nlists}{ SUPERFLUOUS - TO BE REMOVED }
  \item{var.fixed}{ Logical variable - TRUE if you want to tie the variance to be the same across all states. Defaults to FALSE }
}

\value{
 The output is in the same format as that obtained when the nlm function is applied.
}
\author{ John Marioni }

\keyword{ models}

