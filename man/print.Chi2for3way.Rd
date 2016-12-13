\name{print.Chi2for3way}
\alias{print.Chi2for3way}
\title{
Printing function of the results of the three-way index partition 
}
\description{
This function prints the results of the three-way index partition. \cr
The input parameter is the name of the output of the main function  \code{Chi2for3way}. 
}
\usage{
\S3method{print}{Chi2for3way}(x, digits = 3, ...) 
}
\arguments{
\item{x}{The name of the output object, for example say \code{res}, used with the main function \code{Chi2for3way}.}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{The value of output returned depends on the kind of three-way-index and scenario considered
item{X}{The three-way contingency table.}
\item{indexparts}{The three-way index partition \code{indexparts}. When \code{indextype = "chi2"} this output gives the chi-squared 
partition while \code{indextype = "tauM"} returns the partition of Marcotorchino's index, of its numerator and its related $C_M$-statistics. 
Further, it also returns the percentage of explained inertia, the degrees of freedom and the p-value of each term of the partition.}
\item{simulaout}{When the input parameter \code{simulation = TRUE}, the output object \code{simulaout} includes the randomly generated 
three-way tables and their hypothesized probabilities \code{pi}, \code{pj}, \code{pk}. The output also includes the observed marginal frequencies 
\code{margI}, \code{margJ} and \code{margK} and the index partition of $C_M$, $X^2$ and $C^{S}_M$ for each of the \code{nboot} tables.
    
Further, the output also gives the mean of the theoretical and observed distributions of the classical and revised $C_M$-statistics, the variance 
of the theoretical and observed distributions of the classical and revised $C_M$-statistics.
When \code{simulation = FALSE}, then
 \code{simulaout = NULL}.
}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373. \cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
Lombardo R Carlier A D'Ambra L (1996)  Nonsymmetric correspondence analysis for three-way contingency tables. Methodologica,  4, 59-80. \cr
}
\author{
Lombardo R, Takane Y and Beh EJ}
\examples{
data(olive)
resolive<-Chi2for3way(olive,scen=2,simulation=FALSE) 
print.Chi2for3way(resolive)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
