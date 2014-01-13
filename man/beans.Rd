\name{beans}
\alias{beans}
\docType{data}
\title{
Bush beans data set
}
\description{
This experiment was conducted to investigate the effect of row spacing on the yield of different varieties of bush beans. Due to the different growth habits of the considered
varieties it was assumed that the spacing effect differs between the varieties. The
selected four varieties differ such that "New Era" and "Big Green" form low, bushy
plants and the two varieties "Little Gem" and "Red Lake" form erect plants with few
branches. The chosen row spacings were of 20, 40 and 60 cm between rows. A
randomized complete block design with four blocks and 12 plots per block was
used. The yield of dried beans in kilograms per plot was determined after harvest
time.
}
\usage{data(beans)}
\format{
  A data frame with 48 observations on the following 4 variables.
  \describe{
    \item{\code{Variety}}{a factor with levels \code{BigGreen} \code{LittleGem} \code{NewEra} \code{RedLake}, indicating the varieties of bush beans}
    \item{\code{Spacing}}{a numeric vector, indicating the row spacing}
    \item{\code{Block}}{a factor with levels \code{I}, \code{II}, \code{III}, \code{IV}, indicator variable for the 4 blocks}
    \item{\code{Yield}}{a numeric vector, indicating the yield of dried beans in kilograms per plot}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Petersen, R. G. (1985). Design and analysis of experiments, volume 66 of
Statistics. Dekker, New York NY u.a.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(beans)
str(beans) 
}
\keyword{datasets}
