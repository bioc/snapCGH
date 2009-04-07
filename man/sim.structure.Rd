\name{simulateData}
\alias{simulateData}
\alias{generate.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function for simulating aCGH data and the corresponding clone layout}
\description{
  This simulation scheme operates in two stages. Initially, we simulate the layout of clones before using a modified version of the scheme developed by Willenbrock et al., 2005 to generate the \eqn{\log_2} ratios. For each simulated clone layout we generate 20 sets of simulated \eqn{\log_2} ratios from one of five templates. Additionally, we also take account of the cellularity of the test sample in our simulation.
}
\usage{
simulateData(nArrays = 20, chrominfo = NULL, prb.short.tiled = 0.5,
                 prb.long.tiled = 0.5, non.tiled.lower.res = 0.9,
                 non.tiled.upper.res = 1.1, length.clone.lower = 0.05,
                 length.clone.upper = 0.2, tiled.lower.res = -0.05,
                 tiled.upper.res = 0, sd = NULL, output = FALSE,
                 prb.proportion.tiled = c(0.2, 0.2, 0.2, 0.2, 0.2),
                 zerolengthnontiled = NULL, zerolengthtiled = NULL,
                 nonzerolengthnontiled = NULL, nonzerolengthtiled =
                 NULL, seed = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nArrays}{ The number of arrays we want to simulate }
  \item{chrominfo}{ The information about chromosome length/centromere location to be used. Defaults to the information provided in aCGH package of Jane Fridlyand and Peter Dimitrov. }
  \item{prb.short.tiled}{ The probability of a tiled region on the short arm of the simulated chromosome (defaults to 0.5).}
  \item{prb.long.tiled}{ The probability of a tiled region on the long arm of the simulated chromosome (defaults to 0.5).}
  \item{non.tiled.lower.res}{ The lower limit for the distance (in Mbs) between adjacent clones in non-tiled regions of the genome (defaults to 0.9Mb).}
  \item{non.tiled.upper.res}{ The upper limit for the distance (in Mbs) between adjacent clones in non-tiled regions of the genome (defaults to 1.1Mb). }
  \item{length.clone.lower}{ The lower limit for the length (in Mbs) of a clone (this defaults to 0.05Mb). }
  \item{length.clone.upper}{ The upper limit for the length (in Mbs) of a clone (this defaults to 0.2Mb). }
  \item{tiled.lower.res}{ The lower limit for the distance (in Mbs) between adjacent clones in tiled regions of the genome (defaults to -0.05Mb). }
  \item{tiled.upper.res}{ The upper limit for the distance (in Mbs) between adjacent clones in tiled regions of the genome (defaults to 0Mb). }
  \item{sd}{ The standard deviation of the simulated data in each of the states. Defaults to being randomly sampled between 0.1 and 0.2. }
  \item{output}{ A logical variable which is TRUE if you want the output to be written to txt files in the present working directory. Defaults to FALSE. }
  \item{prb.proportion.tiled}{ Given that an arm of a chromosome contains a tiled region this variable (which is a vector of length 5) gives the probability that 20,30,40,50 or 100\% of the chromosome is tiled. It defaults to (0.2,0.2,0.2,0.2,0.2) }
  \item{zerolengthnontiled}{ The empirical distribution for regions of the genome which are non-tiled and contain no copy number gains or losses. Defaults to zero.length.distr.non.tiled }
  \item{zerolengthtiled}{ The empirical distribution for regions of the genome which are tiled and contain no copy number gains or losses. Defaults to zero.length.distr.tiled }
  \item{nonzerolengthnontiled}{ The empirical distribution for regions of the genome which are non-tiled and contain no copy number gains or losses. Defaults to non.zero.length.distr.non.tiled }
  \item{nonzerolengthtiled}{ The empiricial distribution for regions of the genome which are tiled and contain copy number gains or losses. Defaults to non.zero.length.distr.tiled }
  \item{seed}{Seed value allowing simulation to be reproduced if the same seed value is set.}
}
\details{
  For more details see the article by Marioni and Thorne published in Bioinformatics.
}
\value{
  The function returns a list containing the following elements.
  \item{clones}{Gives the start, end and midpoint of the simulated clones.}
  \item{class.output }{A list of the true underlying state clones are assigned to for each of the twenty simulations associated with each clone layout.}
  \item{class.matrix}{Defines the true underlying state clones are assigned to in each of the five classes}
  \item{classes}{Which of the five class outputs has been used to simulate the \eqn{\log_2} ratios}
  \item{datamatrix}{A matrix containing twenty columns each of which contains the simulated \eqn{\log_2} ratios associated with each of the simulations for a particular clone layout.}
  \item{samples}{Gives information about the cellularity associated with each of the samples.}
}
\references{ See the relevant article in Bioinformatics or the following website: www.damtp.cam.ac.uk/user/jcm68 }
\author{ John Marioni and Natalie Thorne }

\examples{
}
\keyword{ datasets }% at least one, from doc/KEYWORDS

