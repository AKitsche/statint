\name{lettuce}
\alias{lettuce}
\docType{data}
\title{
Lettuce data set
}
\description{
Simulated data, inspired by the example presented in Bradu and Gabriel.  
}
\usage{lettuce}
\format{
  A data frame with 48 observations on the following 3 variables.
  \describe{
    \item{\code{Soil}}{a factor with levels \code{S1} \code{S2} \code{S3}, indicating the soil type}
    \item{\code{Fertilizer}}{a factor with levels \code{cal.phos} \code{control} \code{magn.phos} \code{pot.metaphos}, indicating the type of phosphate fertilizer}
    \item{\code{Weight}}{a numeric vector, indicating dry matter in grams per plot}
  }
}
\details{
The experiment was conducted to analyse the effects of soil type and
phosphate fertilizers on lettuce crops. The primary response variable was dry matter
measured in grams per plot. Three different soil types, namely S1, S2 and S3, and four different levels of phosphate fertilization (including an untreated control) were
investigated in a balanced, completely cross-classified treatment structure, laid out
as completely randomized design with four replications per treatment combination.
}
\source{
Bradu, D. and Gabriel, K. R. (1974). Simultaneous Statistical Inference on
Interactions in Two-Way Analysis of Variance. Journal of the American
Statistical Association, 69(346):428-436.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(lettuce)
str(lettuce) 
}
\keyword{datasets}
