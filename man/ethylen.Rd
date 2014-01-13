\name{ethylen}
\alias{ethylen}
\docType{data}
\title{
Ethylene data set
}
\description{
Ethylen data set from an experiment to analyse ethylene receptor blockers for improvement of postharvest characteristics of Pelargonium zonale hybrids.
}
\usage{ethylen}
\format{
  A data frame with 90 observations on the following 9 variables.
  \describe{
    \item{\code{Cult}}{a numeric vector}
    \item{\code{Eth}}{factor with 2 levels, "no" and "yes" specifying the presence or absence of ethylen treatment}
    \item{\code{C1MCP}}{factor with 5 levels, "control", "prod1",   "prod2",   "prod3",   "water"}
    \item{\code{Hue4}}{blabla content after 4 days}
    \item{\code{Hue8}}{blabla content after 8 days}
    \item{\code{Chroma4}}{blabla content after 4 days}
    \item{\code{Chroma8}}{blabla content after 8 days}
    \item{\code{Chloro4}}{Chlorophyll content after 4 days}
    \item{\code{Chloro8}}{Chlorophyll content after 8 days}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
MSc thesis from Syariful Mubarok "Use of ethylene receptor blockers for improvement of postharvest characteristics of Pelargonium zonale hybrids" at the Section of Floriculture within the Institute of Floriculture and Woody Plant Science in the Leibniz Universitaet Hannover.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ethylen)
str(ethylen)
ethylen$Cult  <- factor(ethylen$Cult)
library(ggplot2)
p <- ggplot(ethylen, aes(Treatment, Hue4))
p + geom_boxplot(aes(fill= Cult)) + facet_grid(. ~ Eth)

}
\keyword{datasets}
