\name{plotSummaryProfile}
\alias{plotSummaryProfile}
\alias{length.num.func}
\alias{prop.num.func}
\title{plotSummaryProfile}
\description{
  This function display the genomic events and tests for the differences
  between groups if they are specified.
}
\usage{
plotSummaryProfile(input, geList, response = as.factor(rep("All",
                 ncol(input))), titles =
                 unique(response[!is.na(response)]), X = TRUE, Y =
                 FALSE, maxChrom = 23, chrominfo = chrominfo.basepair,
                 num.plots.per.page = length(titles), factor = 2.5,
                 thresAbs = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{an object of either \code{\link[snapCGH:SegList]{SegList}} or
    \code{\link[limma]{MAList}} class}
  \item{geList}{an object of \code{\link[snapCGH:GEList]{GEList}} class}
  \item{response}{
    phenotype to compare. defaults to all the samples being analyzed
    together.
  }
  \item{titles}{
    titles for the groups, defaults to the name of the \code{response}
  }
  \item{X}{
    logical indicating whether X needs to be shown
  }
  \item{Y}{
    logical indicating whether Y needs to be shown
  }  
  \item{maxChrom}{this parameter controls how many chromosomes will be
    plotted, from 1 to maxChrom}
  \item{chrominfo}{
    a chromosomal information associated with the mapping of the data
  }
  \item{num.plots.per.page}{
    number of frequency plots per page. Default is the number of
    groups
  }
  \item{factor}{\code{factor} specifies the number by which experimental variability should be multiples. Used only when tumor specific variability in \code{\link[snapCGH:SegList]{SegList}} is not NULL. Defaults to 2.5}
  \item{thresAbs}{Threshold for gain. Set very high for homozygous deletion}
}
\details{
  This function utilizes output of the \code{\link{find.genomic.events}}
  by plotting it and testing between groups. The test are performed
  using kruskal-wallis rank test.
}
% \value{
%   plot containing a number of events
% }
\seealso{
  \code{\link[limma]{MAList}}
  \code{\link{find.genomic.events}}
}

}
\keyword{hplot}% at least one, from doc/KEYWORDS
