\name{chi3scen1}
\alias{chi3scen1}
\title{
Pearson's index for three-way contingency tables under Scenario 1 (prescribed probabilities)
}
\description{
It provides the Pearson's index, e.g. chi-square index, partitioning under the Scenario 1 when probabilities are homogeneous.
}
\usage{
chi3scen1(X, pi=rep(1/dim(X)[[1]],dim(X)[[1]]), pj=rep(1/dim(X)[[2]],dim(X)[[2]]), 
pk=rep(1/dim(X)[[3]],dim(X)[[3]]), digits = 3)
}
\arguments{
  \item{X}{The three-way contingency table.
}
\item{pi}{The input parameter for specifying the theoretical probabilities of rows categories. 
When \code{scen = 1}, they can be prescribed by the analyst. By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is 
\code{pi = rep(1/dim(X)[[1]],dim(X)[[1]])}.}
\item{pj}{The input parameter for specifying the theoretical probabilities of  columns categories. 
When \code{scen = 1}, they can be prescribed by the analyst. By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is
\code{pj = rep(1/dim(X)[[2]],dim(X)[[2]])}.}
\item{pk}{The input parameter for specifying the theoretical probabilities of  tube categories. 
When \code{scen = 1}, they can be prescribed by the analyst. By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is
\code{pk = rep(1/dim(X)[[3]],dim(X)[[3]])}.}
  \item{digits}{
The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.
}
}
\value{
Description of the output returned
\item{z}{The chi-square index partition under Scenario 1, we get seven terms of the chi-square partition, three main terms,
two bivariate terms and a trivariate term.
The output is in a matrix, the four rows of this matrix indicate
the index, the  percentage of the explained inertia, the degree of freedom, the p-value, respectively.}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373. \cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R and Takane Y}
\note{
This function belongs to the class \code{chi3class}. 
}
\examples{
##---- Should be DIRECTLY executable !! ----
data(olive)
chi3scen1(olive)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{chi-square}% use one of  RShowDoc("KEYWORDS")
\keyword{index}% __ONLY ONE__ keyword per line
