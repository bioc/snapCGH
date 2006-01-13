\name{threshold}
\alias{threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to indicate gain or loss for each clone for each sample}
\description{
 This function outputs a matrix containing gain/loss indicator for each clone and sample. 
}
\usage{
threshold(dat, thresAbs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{log2ratios of the relevant array CGH object}
  \item{thresAbs}{Global or sample-specific threshold for gain}
  }
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
     Returns a matrix with a row for each clone and column for each sample. "1" indicates gain and "-1" indicates loss.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jane Fridlyand, Ritu Roydasgupta}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{\code{\link{plotGainLoss}}}
\examples{

%data(colorectal)

%factor <- 3
%tbl <- threshold(log2ratios(colorectal),posThres=factor*(sd.samples(colorectal)$madGenome))
%rownames(tbl) <- clone.names(colorectal)
%colnames(tbl) <- sample.names(colorectal)
%tbl[1:5,1:5]

}
\keyword{htest}

