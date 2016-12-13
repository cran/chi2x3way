\name{chi3scen2}
\alias{chi3scen2}
\title{
Pearson's index for three-way contingency tables under Scenario 2
}
\description{
It provides the Pearson's index, e.g. chi-square index, partitioning under the Scenario 2 when probabilities are set equal to the observed margin frequencies}
\usage{
chi3scen2(X, digits = 3)
}
\arguments{
  \item{X}{The three-way contingency table.
}
  \item{digits}{
The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.
}
}
\value{
 \item{z }{The chi-square index partition under Scenario 2, we get five terms of the chi-square partition, three bivariate terms and a trivariate one. 
The output is in a matrix, the four rows of this matrix indicate
the index, the percentage of the explained inertia, the degree of freedom, the p-value, respectively.}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373.\cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R, Takane Y and Beh EJ}
\note{
This function belongs to the class \code{chi3class}. 
}

\examples{
##---- Should be DIRECTLY executable !! ----
data(olive)
chi3scen2(olive)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{chi-square}% use one of  RShowDoc("KEYWORDS")
\keyword{index}% __ONLY ONE__ keyword per line
