\name{tau3scen2}
\alias{tau3scen2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marcortchino's index for three-way contingency tables under Scenario 2
}
\description{
It provides the Marcotorchino' index partitioning under  Scenario 2 when probabilities are set equal to the observed marginal frequencies.
}
\usage{
tau3scen2(X, digits = 3)
}
\arguments{
\item{X}{The three-way contingency table.}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
}
\value{
\item{z}{The Marcotorchino's index partition under Scenario 2, we get five terms of the tau index partition, three bivariate terms and a trivariate one. 
The output is in a matrix, , the six rows of this matrix indicate the tau index numerator,  the tau index, the percentage of explained inertia,
 the $C_M$-statistic, the degree of freedom, the p-value, respectively.}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Lombardo R Carlier A D'Ambra L (1996) Nonsymmetric correspondence analysis for three-way contingency tables. Methodologica,  4, 59-80.\cr
Loisel S and Takane Y (2015) Partitions of Pearson's chi-square statistic for frequency tables: A comprehensive account.  Computational Statistics, 31, 1429-1452. \cr
Marcotorchino F (1985) Utilisation des comparaisons par paires en statistique des contingencies: Partie III. Etude du Centre Scientifique, IBM, France. No F 081 
}
\author{
Lombardo R, Takane Y and Beh EJ
}
\note{
This function belongs to the class \code{chi3class}. 
}
\examples{
##---- Should be DIRECTLY executable !! ----
data(olive)
tau3scen2(olive)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marcotorchino index}% use one of  RShowDoc("KEYWORDS")
\keyword{three-way index}% __ONLY ONE__ keyword per line
\keyword{partition}