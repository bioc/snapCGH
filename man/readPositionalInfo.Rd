\name{readPositionalInfo}
\alias{readPositionalInfo}

\title{readPositionalInfo}
\description{
  This function automatically inserts information about the chromosomal
  postional of a clone into the \$genes matrix of an RGList or MAList.
  This information is used in all the available segmentation methods as
  well as many of the plotting functions available in snapCGH.
}
\usage{
readPositionalInfo(input, source, path = NULL)
}
\arguments{
  \item{input}{An object of class \code{\link[limma]{RGList}} or
    \code{\link[limma]{MAList}}}
  \item{source}{Defines which platform or technology this data was
    produced on.  Currently supported options are: "aglient",
    "bluefuse", "nimblegen".  This list will be expanded in time.}
  \item{path}{Optional parameter to specify where the original data is
    stored. Defaults to the current working directory.  This option is
    only required for reading "bluefuse" data at the moment, as
    chromosome information isn't read by limma as default.}
}

\keyword{manip}

