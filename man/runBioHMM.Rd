\name{runBioHMM}
\alias{runBioHMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ This function implements the BioHMM}
\description{
  This function reads in a dataset of log2 ratios and the corresponding clone and covariate information. It calculates a heterogeneous HMM when there are 
  1,2,3,4 or 5 underlying states and chooses between them using either the AIC or BIC. It then assigns clones using a modified version of the Viterbi
  algorithm.
}
\usage{
runBioHMM(input, useCloneDists = TRUE, covariates, criteria="AIC", delta=NA
,var.fixed=FALSE, epsilon = 1e-06, numiter = 30000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{An object of class \code{\link[limma:malist]{MAList}}or \code{\link[snapCGH:SegList]{SegList}}}
  \item{useCloneDists}{Boolean stating whether the distance between
    clones should be incorportated into the HMM.  If false then the HMM
    becomes homogeneous.}
  \item{covariates}{ This is a dataframe containing information about covariate factors. The first two columns should be Chrom (giving the 
                     chromosome on which a clone is located) and Mb (giving the position of the chromosome along a particular chromosome in Megabases).
                     The order should be the same as that described above with the following crucial difference. No covariate information about the first
                     clone is used in the segmentation. Hence, for each chromosome, there should be one less row in the covariate dataframe than in the 
                     datainfo dataframe corresponding to this missing chromosome. This is important if the transition matrix is to be calculated 
                     correctly. }
  \item{criteria}{Options are AIC or BIC depending upon which we want to use to distinguish between the number of states}
  \item{delta}{ A variable to be assigned if the BIC is used. }
  \item{var.fixed}{ Logical variable - TRUE if you want to tie the
    variance to be the same across all states. Defaults to FALSE}
  \item{epsilon}{Stopping criterion for the optimization algorithm.}
  \item{numiter}{Number of iterations to be used in the
    optimization algorithm.}
}
\value{
  The model returns an object of class \code{\link[snapCGH:SegList]{SegList}}.
}
\references{ Marioni, J.C., Thorne, N.P., Tavar\'e, S., BioHMM: a heterogeneous Hidden Markov Model for segmenting array CGH data, submitted}
\author{ John Marioni and Mike Smith}

\keyword{ models}
