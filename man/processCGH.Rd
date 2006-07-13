\name{processCGH}
\alias{processCGH}
\title{Process data in an MAList}
\description{
  This function takes object of class \code{MAList} and it re-orderes
  and filters clones based on their mapping information and proportion missing. It also average duplicated clones and imputes missing values for clones
  that are still NA after the filetering step.
}
\usage{
processCGH(input, maxChromThreshold = 22, minChromThreshold
                 = 1, method.of.averaging = NULL, ID = "ID")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{Object of class \code{\link[limma]{MAList}} or \code{\link[limma]{RGList}}}
  \item{maxChromThreshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 repsectivelly, in for Mus musculus 20
    and 21, etc. Remove chromosomes from segmentation analysis which are
                 greater than this value.}
  \item{minChromThreshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 repsectivelly, in for Mus musculus 20
    and 21, etc. Remove chromosomes from segmentation analysis which re
                 lower than this value.}
   \item{method.of.averaging}{If left as the default no combining of replicate spots takes place.  Other wise this should specify a function which takes a vector of
                 duplicates and combines them into a single value.}
   \item{ID}{Name of column in \code{RG\$genes} corresponding to the clone names. For most
                 data the default will work, however for affy data the
                 value for ID should be \code{"CloneName"}}
	       }
\value{
  Object of class \code{\link[snapCGH]{SegList}}
}
\author{Jane Fridlyand, Peter Dimitrov, John Marioni and Mike Smith}
\seealso{
  \code{\link[limma]{MAList}}
}

}
\keyword{file}% at least one, from doc/KEYWORDS
