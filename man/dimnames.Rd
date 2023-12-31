\name{dimnames}
\alias{dimnames}
\alias{dimnames,SegList-method}
\alias{dimnames.SegList}
\title{Retrieve the Dimension Names of an RGList, MAList or SegList Object}
\description{
Retrieve the dimension names of a microarray data object.
}
\usage{
\method{dimnames}{SegList}(x)
}
\arguments{
  \item{x}{An object of class \code{SegList}}
}
\details{
The dimension names of an microarray object are the same as those of the most important matrix component of that object.

A consequence is that \code{rownames} and \code{colnames} will work as expected.
}
\value{
Either \code{NULL} or a list of length 2.
If a list, its components are either 'NULL' or a character vector the length of the appropriate dimension of \code{x}.
}
\author{Gordon Smyth, edited by Mike Smith}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
}
\keyword{array}
