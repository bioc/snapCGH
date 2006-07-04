\name{runHomHMM}
\alias{runHomHMM}
\alias{mu1.func}
\alias{states.hmm.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{A function to fit unsupervised Hidden Markov model}
\description{
  This function fits an unsupervised Hidden Markov model to a given \code{\link[limma]{MAList}}
}
\usage{
runHomHMM(MA, vr = 0.01,
                maxiter = 100, criteria = "AIC", delta = NA,
                full.output = FALSE, eps = 0.01)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{MA}{an object of class \code{\link[limma]{MAList}}}
  \item{vr}{Gets passed to the function
    \code{\link[pkg:repeated]{hidden}} as the \code{pshape} arguement.}
  \item{maxiter}{Gets passed to the function
    \code{\link[pkg:repeated]{hidden}} as the \code{iterlim}
    arguement. }
  \item{criteria}{Choice of which selection criteria should be used in
    the algorithm.  The choices are either AIC or BIC}.
  \item{delta}{Delta value used of the BIC is selected.  If no value is
    entered it defaults to 1.}
  \item{full.output}{if true the SegList output includes a probability
    that a clone is in its assigned state and a smoothed value for the
    clone.}
  \item{eps}{parameter controlling the convergence of the EM
    algorithm. }
  }

\seealso{
  \code{\link{runDNAcopy}}
  \code{\link{runGLAD}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\keyword{models}% at least one, from doc/KEYWORDS

