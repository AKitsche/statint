\name{potato}
\alias{potato}
\docType{data}
\title{
Potato data set
}
\description{
Data set from an experiment to investigate a potential differing response among nine varieties of starch potato (Solanum tuberosum) on water deficiency stress.
}
\usage{potato}
\format{
  A data frame with 216 observations on the following 5 variables.
  \describe{
    \item{\code{genotype}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G} \code{H} \code{I}, identifying the variety  of starch potato (Solanum tuberosum)}
    \item{\code{treatment}}{a factor with levels \code{0.1} \code{0.2} \code{0.3} \code{control}, specifying the water deficiency condition}
    \item{\code{length}}{a numeric vector,  the length of the plants}
    \item{\code{shoot}}{a numeric vector, the shoot fresh weight}
    \item{\code{root}}{a numeric vector, the root fresh weight}
  }
}
\details{
A total number of nine varieties of starch potato was selected for cultivation. Each variety was cultivated under four different water conditions: standard control medium, 0.1M, 0.2M and 0.3M sorbitol medium. An increasing sorbitol concentration was used to induce an increasing water deficiency stress on the plants. For each treatment by variety combination 6 replicates were investigated. The trial was planned in a completely cross-classified treatment structure, in a completely randomized design. For each experimental unit the length, the fresh weight of the shoot and the fresh weight of the root was measured. 
}
\source{
This data set was part of an experiment conducted at the Woody Plant and Propagation Physiology Section in the Institute of Horticultural Production Systems at the Leibniz Universitaet Hannover in 2013 (from M.Sc. Christin Hadden, Prof. Dr. Traud Winkelmann). 
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(potato)
str(potato) 

}
\keyword{datasets}
