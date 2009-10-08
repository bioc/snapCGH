\name{read.clonesinfo}
\alias{read.clonesinfo}
\title{Reading chromsome and positional information about each clone.}
\description{
  Function to read the chromosomal position information of each clone
  and incorporate it into the genes data.frame of the relevant object.
}

\usage{
read.clonesinfo(file, RG, path = NULL, sep="\t", quote="\"")
}

\arguments{
  \item{file}{Name of the file containing the chromosomal information.}
  \item{RG}{Object containing a \$genes data.frame that the information
    should be incorporated into.}
  \item{path}{Path to the chromosomal information file.}
  \item{sep}{Identifying the column seperator in the designated file.}
  \item{quote}{Identifying the quotation character used in the
    designated file.}
}

\author{Mike Smith}

\keyword{methods}
