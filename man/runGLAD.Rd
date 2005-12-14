\name{runGLAD}
\alias{runGLAD}
\title{Results of segmenting an aCGHList data object using the GLAD library}
\description{
  This function allows the detection of breakpoints in genomic profiles
  obtained by array CGH technology and affects a status (gain, normal
  or lost) to each clone.  It requires that the library \code{GLAD} is loaded.
}
\usage{
runGLAD(MA, smoothfunc="aws", base=FALSE, sigma, bandwidth=10,
round=2, lambdabreak=8, lambdacluster=8, lambdaclusterGen=40,
type="tricubic", param=c(d=6), alpha=0.001, method="centroid",
nmax=8, verbose=FALSE, ...)
}
\arguments{
  \item{MA}{An object of class \code{\link[limma]{MAList}}}
  
  \item{smoothfunc}{Type of algorithm used to smooth \code{LogRatio} by a
    piecewise constant function. Choose either \code{lawsglad}, \code{\link[aws]{aws}} or
    \code{\link[aws]{laws}}.}

  \item{base}{If \code{TRUE}, the position of clone is the physical position onto
    the chromosome, otherwise the rank position is used.}

  \item{sigma}{Value to be passed to either argument \code{sigma2}   
    of\code{ \link[aws]{aws}} function or \code{shape} of
    \code{\link[aws]{laws}}. If \code{NULL}, sigma is calculated from
    the data.}
  
  \item{bandwidth}{Set the maximal bandwidth \code{hmax} in the
    \code{\link[aws]{aws}} or  \code{\link[aws]{laws}} function. For
    example, if \code{bandwidth=10} then the \code{hmax} value is set
    to 10*\eqn{X_N} where \eqn{X_N} is the position of the last clone.}

  \item{round}{The smoothing results are rounded or not depending on
    the \code{round} argument. The \code{round} value is passed to the
    argument \code{digits} of the \code{\link[base]{round}} function.}

  \item{lambdabreak}{Penalty term (\eqn{\lambda'}) used during the 
      \bold{Optimization of the number of breakpoints} step.}

  \item{lambdacluster}{Penalty term (\eqn{\lambda*}) used during the \bold{MSHR
      clustering by chromosome} step.}
  
  \item{lambdaclusterGen}{Penalty term (\eqn{\lambda*}) used during the \bold{HCSR
      clustering throughout the genome} step.}

  \item{type}{Type of kernel function used in the penalty term during the \bold{Optimization of the
      number of breakpoints} step, the \bold{MSHR
      clustering by chromosome} step and the \bold{HCSR
      clustering throughout the genome} step.}
  
  \item{param}{Parameter of kernel used in the penalty term.}

  \item{alpha}{Parameter for ?????.}
  
  \item{method}{The agglomeration method to be used during the \bold{MSHR
      clustering by chromosome} and the \bold{HCSR
      clustering throughout the genome} clustering steps.}
  
  \item{nmax}{Maximum number of clusters (N*max) allowed during
    the the \bold{MSHR
      clustering by chromosome} and the \bold{HCSR
      clustering throughout the genome} clustering steps.}

  \item{verbose}{If \code{TRUE} some information are printed}  
  
  \item{...}{}
}
\details{
  For a detailed explanation of the GLAD algorithm please see the
  relevant section of the GLAD manual: \code{\link[GLAD]{glad}}
  }
\value{
  The function returns an object of class \code{\link[snapCGH:SegList]{SegList}}
}
\seealso{
  \code{\link[GLAD]{glad}}
  \code{\link[limma]{MAList}}
  \code{\link{fitHMM}}
  \code{\link{runDNAcopy}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\keyword{methods}
