\name{dim}
\alias{dim}
\alias{length}
\alias{dim.RGList}
\alias{dim.MAList}
\alias{dim.SegList}
\alias{dim,SegList-method}
\alias{length.RGList}
\alias{length.MAList}
\alias{length.SegList}
\alias{length,SegList-method}
\title{Retrieve the Dimensions of an RGList, MAList or SegList Object}
\description{
Retrieve the number of rows (genes) and columns (arrays) for an RGList, MAList or SegList object.
}
\usage{
\method{dim}{SegList}(x)
\method{length}{SegList}(x)
}
\arguments{
  \item{x}{an object of class \code{RGList}, \code{MAList} or \code{SegList}}
}
\details{
Microarray data objects share many analogies with ordinary matrices in which the rows correspond to spots or genes and the columns to arrays.
These methods allow one to extract the size of microarray data objects in the same way that one would do for ordinary matrices.

A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} and so on also work.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{Gordon Smyth, modified by Mike Smith for SegList object}
\seealso{
  \code{\link[base]{dim}} in the base package.
  
%  \link{2.Classes} gives an overview of data classes used in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA <- new("MAList",list(M=M,A=A))
dim(M)
ncol(M)
nrow(M)
length(M)
}
\keyword{array}
