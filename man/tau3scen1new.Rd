\name{tau3scen1new}
\alias{tau3scen1new}
\title{
Marcotorchino's index for three-way contingency tables under Scenario 1. Revised formulation.
}
\description{
It provides the partition of the Marcotorchino's index and of the $C_M$-statistic, revised formula,  under the Scenario 1 when probabilities are set homogeneous.
}
\usage{
tau3scen1new(X, pi=rep(1/dim(X)[[1]],dim(X)[[1]]), pj=rep(1/dim(X)[[2]],dim(X)[[2]]), 
pk=rep(1/dim(X)[[3]],dim(X)[[3]]), digits = 3)
}
\arguments{
\item{X}{The three-way contingency table.}
\item{pi}{The input parameter for specifying the theoretical probabilities of rows categories. When \code{scen = 1}, they can be prescribed by the analyst. \cr
By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is \code{pi = rep(1/dim(X)[[1]],dim(X)[[1]])}.
}
\item{pj}{The input parameter for specifying the theoretical probabilities of  columns categories. When \code{scen = 1}, they can be prescribed by the analyst. \cr
By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is \code{pj = rep(1/dim(X)[[2]],dim(X)[[2]])}.
}
\item{pk}{The input parameter for specifying the theoretical probabilities of  tube categories. When \code{scen = 1}, they can be prescribed by the analyst. \cr
By default, they are set equal among the categories, homogeneous margins (uniform probabilities), that is \code{pk = rep(1/dim(X)[[3]],dim(X)[[3]])}.
}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. By default, \code{digits = 3}.
}
}
\value{
Description of  of the output returned
\item{z}{The Marcotorchino's index partition under Scenario 1, we get seven terms partitioning  the Marcotorchino's index and  the revised $C_M$-statistic, 
three main terms, two bivariate terms and a trivariate term.
The output is in a matrix, the six rows of this matrix indicate the tau index numerator,  the tau index, the percentage of explained inertia, the revised $C_M$-statistic,
the degree of freedom, the p-value, respectively.}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2015) Partitions of Pearson's chi-square statistic for frequency tables: A comprehensive account.  Computational Statistics, 31, 1429-1452.\cr
Lombardo R Carlier A D'Ambra L (1996) Nonsymmetric correspondence analysis for three-way contingency tables. Methodologica,  4, 59-80.\cr
Marcotorchino F (1985) Utilisation des comparaisons par paires en statistique des contingencies: Partie III. Etude du Centre Scientifique, IBM, France. No F 081 
}
\author{
Lombardo R, Takane Y and Beh EJ}
\note{
This function belongs to the class \code{chi3class}. 
}
\examples{
data(olive)
tau3scen1new(olive)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marcotorchino index}% use one of  RShowDoc("KEYWORDS")
\keyword{three-way index}