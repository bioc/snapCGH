\name{fit.model}
\alias{fit.model}
\alias{run.nelder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting a heterogneous HMM to the log2 ratios on a particular chromosome.}
\description{
  This function fits five homogeneous HMMs to the log2 ratios on a particular chromosome. It then uses either the
  AIC or BIC to determine which of the five models is optimal before using a scaled version of the Viterbi algorithm to assign clones to states with the
  same underlying copy number.
}
\usage{
fit.model(sample, chrom, dat, datainfo = clones.info, useCloneDists = TRUE, covariates,
aic = TRUE, bic = FALSE, delta = 1, var.fixed=FALSE, epsilon = 1e-06,
numiter = 30000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{ If there are multiple samples, the number of the sample to be segmented }
  \item{chrom}{ The chromosome on which the segmentation is to be carried out on }
  \item{dat}{ The log2 ratios obtained from the clones located on that chromosome }
  \item{datainfo}{ A dataframe containing information about the clones
    on that chromosome (name, chromosome and location (in Mbs)) }
  \item{useCloneDists}{Boolean stating whether the distance between
    clones should be incorportated into the HMM.  If false then the HMM
    become homogeneous.}
  \item{covariates}{ A matrix contianing the covariate information for the clones located on the chromosome to be segmented. It should have length one less
                     than the number of clones as covariate information is not used when segmenting the first clone on the chromosome. }
  \item{aic}{Set to true if you want to use the aic. This is the
    default.  Only one of aic and bic should be set to true.}
  \item{bic}{Set to true if you want to use the bic.}
  \item{delta}{ A parameter to be set if you want to use the BIC }
  \item{var.fixed}{ Logical variable - TRUE if you want to tie the
    variance to be the same across all states. Defaults to FALSE }
  \item{epsilon}{.}	
  \item{numiter}{Number of iterations to be used in the
    optimization alogithm.}
}

\value{
 The output is in the same format as that obtained when the nlm function is applied.
}
\author{John Marioni and Mike Smith}

\keyword{ models}

