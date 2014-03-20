\name{plotcontrast}
\alias{plotcontrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot contrast matrix in a tile plane with rectangles
}
\description{
This function plots a matrix of contrast in a tile plane with rectangles. Is is very useful to disply more complicated contrast matrices like product type interaction contrast matrices.
}
\usage{
plotcontrast(cmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmat}{
a contrast matrix with rownames a the names comparisons of interest and colnames as the levels of a factor
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Frank Schaarschmidt
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fa<-factor(rep(rep(LETTERS[1:5], rep(4,5)), each=3))
fb<-factor(rep(rep(c(1,2,3,4), rep(3,4)), times=5))

cma<-contrMat(table(fa), type="Dunnett")
cmb<-contrMat(table(fb), type="Tukey")

InteractionMatrix<-iacontrast(fa=fa, fb=fb, cma=cma, cmb=cmb, method="c")

plotcontrast(InteractionMatrix$cmab)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
