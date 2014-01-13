\name{iacontrast}
\alias{iacontrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating labelled product type interaction matrix
}
\description{
Function to create a product type interaction contrast matrix from two one-way contrast matrices according to Gabriel et al.  The two one-way contrast matrices are constructed from a prespecified set of contrast or from user defined contrast matrices. 
}
\usage{
iacontrast(fa, fb, typea = "Dunnett", typeb = "Dunnett", 
           cma = NULL, cmb = NULL, droplevels = TRUE, 
           abbrevnames = NULL, orderby = NULL, 
           sep = " - ", cw = ",", cb = ":", method = "at")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa}{
factor variable for the first factor (A) with at least two levels
}
  \item{fb}{
factor variable for the second factor (B) with at least two levels
}
  \item{typea}{
character string, one of those as described for type in contrMat (multcomp), ignored if cma is given. If provided, typea is ignored.
}
  \item{typeb}{
character string, one of those as described for type in contrMat (multcomp), ignored if cmb is given. If provided, typeb is ignored.
}
  \item{cma}{
contrast matrix, with as many columns as there are levels in fa (before or after dropping unused levels, see droplevels below).
}
  \item{cmb}{
contrast matrix, with as many columns as there are levels in fb (before or after dropping unused levels, see droplevels below).
}
  \item{droplevels}{
logical, indicating whether unused levels of the input factors fa, fb should be dropped. If TRUE, the number of columns in cma, cmb, should fit the number of used levels. If FALSE, the number of columns in cma, cmb shoudl fit the number of level in the leve?l attribute of the factors.
}
  \item{abbrevnames}{
a named list of arguments to be passed to abbreviate the factor level names for the construction of names for the interaction contrasts
}
  \item{orderby}{
CURRENTLY NOT USED a,b, whether the output contrasts should be additionally ordered according to factor a, b 
}
  \item{sep}{
character string to be used to separate different levels of the cellmeans in comparisons in the contrasts names
}
  \item{cw}{
character string, to be used when collapsing original factor levels within a factor (pooling contrasts)
}
  \item{cb}{
character string, to be used when collapsing original factor levels to cellmean levels in the contrasts names
}
  \item{method}{
a character string specifying the lebelling of contrasts, must be one of "at" (default) or "cellmeans" .
}
}
\details{
Generates a labelled product type interaction contrast matrix from two one way contrast matrices. The user can either use select from a set of prespecified contrast in the function contrMat from the package multcomp. 
}
\value{
a list with two elements:
%%  If it is a LIST, use
  \item{fab }{the cellmeans factor (with one level for each level combination of fa, fb)}
  \item{cmab }{interaction contrast matrix, with columns and naming of contrast fitting the level order in fab.}
%% ...
}
\references{
K. R. Gabriel, J. Putter, and Y. Wax. Simultaneous Confidence Intervals for product-type Interaction Contrasts. Journal
of the Royal Statistical Society Series B - Statistical Methodology, 35(2): 234-244, 1973.

A. Kitsche, F. Schaarschmidt. Analysis of statistical interactions in factorial experiments. Journal of Agronomy and Crop Science, 2014
}
\author{
Frank Schaarschmidt, Andreas Kitsche
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fa<-factor(rep(rep(LETTERS[1:5], rep(3,5)), each=3))
fb<-factor(rep(rep(c(1,2,3), rep(3,3)), times=5))

iacontrast(fa=fa, fb=fb, typea="Dunnett", typeb="Tukey")
iacontrast(fa=fa, fb=fb, typea="Dunnett", typeb="Tukey", method="c")


cmbp<-rbind(c(1,0,0),
            c(0,0.5,0.5))
            
iacontrast(fa=fa, fb=fb, typea="Dunnett", cmb=cmbp)
iacontrast(fa=fa, fb=fb, typea="Dunnett", cmb=cmbp, method="c")

cmap<-rbind(c(1,1,1,0,0)/3,
            c(0,0,0,1,1)/2)
            
iacontrast(fa=fa, fb=fb, cma=cmap, typeb="Tukey")
iacontrast(fa=fa, fb=fb, cma=cmap, typeb="Tukey", method="c")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
