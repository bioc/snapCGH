\name{processCGH}
\alias{processCGH}
\alias{prop.na}
\title{Process data in an MAList}
\description{
  This function takes object of class \code{MAList} and it re-orders
  and filters clones based on their mapping information and proportion missing.
  It also average duplicated clones and imputes missing values for clones
  that are still NA after the filtering step. Note that imputation will
  only take place if duplicated clones are removed.
}
\usage{
processCGH(input, maxChromThreshold = 22, minChromThreshold = 1,
           method.of.averaging = NULL, ID = "ID", prop.missing = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{Object of class \code{\link[limma:malist]{MAList}} or
    \code{\link[limma:rglist]{RGList}}}
  \item{maxChromThreshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 respectively, in for Mus musculus 20
    and 21, etc. Remove chromosomes from segmentation analysis which are
    greater than this value.}
  \item{minChromThreshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 respectively, in for Mus musculus 20
    and 21, etc. Remove chromosomes from segmentation analysis which are
    lower than this value.}
  \item{method.of.averaging}{If left as the default no combining of replicate
    spots takes place. Otherwise this should specify a function which takes a
    vector of duplicates and combines them into a single value.}
  \item{ID}{Name of column in \code{RG\$genes} corresponding to the clone names.
    For most data the default will work, however for affy data the value for ID
    should be \code{"CloneName"}}
  \item{prop.missing}{For each probe the proportion of NA's is
     calculated, and the probe is kept for further analysis if the
     proportion of NA's is less than missing.prop}
}
\value{
  Object of class \code{\link[snapCGH]{SegList}}
}
\author{Jane Fridlyand, Peter Dimitrov, John Marioni and Mike Smith}
\seealso{
  \code{\link[limma:malist]{MAList}}
}

\keyword{file}% at least one, from doc/KEYWORDS
