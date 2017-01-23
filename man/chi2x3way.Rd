\name{chi2x3way}
\alias{chi2x3way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chi-square and Marcotorchino's index for three-way contingency tables
}
\description{
It performs  \cr
1) the computation of the Pearson's index and its partitioning for three-way contingency tables under two Scenarios. When the input parameter \code{scen==1} then the theoretical probabilities are 
prescribed by the analyst (by default they are set homogeneous). When   the input parameter \code{scen==2} then the theoretical probabilities are  estimated from the data.\cr
2) the computation of the Marcotorchino's index and its partitioning for three-way contingency tables under the  two Scenarios. When the input parameter \code{scen==1} then the theoretical probabilities are 
prescribed by the analyst (by default they are set homogeneous). When   the input parameter \code{scen==2} then the theoretical probabilities are  
estimated from the data. In order to check the distribution of the Marcotorchino's index under the two Scenarios, it is possible to look at the results of a simulation study setting the input parameter 
\code{simulation=TRUE}.
}
\usage{
chi2x3way(X,  indextype = "chi2", scen = 2, simulation = FALSE,
nboots = 1000, nran = 1000,
pi = rep(1/dim(X)[[1]],dim(X)[[1]]),
pj = rep(1/dim(X)[[2]],dim(X)[[2]]),pk = rep(1/dim(X)[[3]],dim(X)[[3]]), digits = 3)
}
\arguments{
\item{X}{
The three-way contingency table.}
\item{indextype}{The input parameter  for specifying what index should be considered. 
By default, the partition of the classical three-way 
Pearson index \code{indextype = "chi2"} is selected. The analyst can also partition Marcotorchino's index by defining the input parameter 
\code{indextype = "tauM"}.
}
\item{scen}{The input parameter for specifying what scenario should be considered. By default, \code{scen = 1}, so that the probabilities are defined 
as being fixed and homogeneous among the categories (i.e. {Scenario 1}). When \code{scen = 2}, the expected 
frequencies are set to be equal to the observed marginal frequencies (i.e. {Scenario 2}).}
\item{simulation}{
A flag parameter, \code{simulation}, is included for specifying whether simulations are included as part of the analysis. 
When \code{simulation = TRUE}, three-way contingency tables are randomly generated under the  different scenarios specified by \code{scen}. 
Note that for investigating the index distributions, a randomly generated contingency table which consists of at least one cell frequency that is less than five is automatically discarded. 
When \code{simulation = TRUE}, the distribution of the terms from the partition of the classic $C_M$-statistic, associated with Marcotorchino's index, the revised $C^{S}_{M}$-statistic 
and Pearson's chi-squared index are graphically depicted and compared using QQ-plots. By default, \code{simulation = FALSE}.}
\item{nboots}{The input parameter for specifying the number of random  three-way contingency  tables to be generated when \code{simulation = TRUE}. By default, \code{nboots = 1000}.}
\item{nran}{The input parameter for specifying the total number of samples of each randomly generated contingency table when \code{simulation = TRUE}. By default, \code{nran = 1000}.}
\item{pi}{The input parameter \code{pi} specifies the probabilities assigned to the row categories.
When \code{scen = 1}, they can be arbitrarily defined by the analyst. By default, the parameter 
is set to reflect homogeneous marginal (uniform) probabilities so that \code{pi = rep(1/dim(X)[[1]], dim(X)[[1]])}. 
When \code{scen = 2} the hypothesized probabilities cannot be prescribed by the analyst and are set equal to the observed row margins of the 
three-way table.}
\item{pj}{The input parameter \code{pj} specifies the probabilities assigned to the column categories.
When \code{scen = 1}, they can be arbitrarily 
defined by the analyst. By default, the parameter is set to reflect homogeneous marginal (uniform) probabilities so that 
\code{pi = rep(1/dim(X)[[2]], dim(X)[[2]])}. When \code{scen = 2} the hypothesized probabilities 
cannot be prescribed by the analyst and are set equal to the observed column margins of the three-way table.}
\item{pk}{The input parameter \code{pk} specifies the probabilities assigned to the tube categories.
When \code{scen = 1}, 
they can be arbitrarily defined by the analyst. By default, the parameter is set to reflect homogeneous marginal (uniform) probabilities 
so that \code{pi = rep(1/dim(X)[[3]], dim(X)[[3]])}. When \code{scen = 2} the hypothesized probabilities  are set equal to the observed tube
 margins of the three-way table.}
\item{digits}{The minimum number of decimal places used for displaying the numerical summaries of the analysis is set by the parameter \code{digits}.

By default, \code{digits = 3}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{X}{The three-way contingency table of dimension IxJxK.}
 \item{indexparts}{The three-way index partition \code{indexparts}. When \code{indextype = "chi2"},
this output gives the chi-squared partition, while \code{indextype = "tauM"} returns the partition of Marcotorchino's index, 
$tau_M$ and its related $C_M$-statistics. 
Further, it also returns the percentage of explained inertia, the degrees of freedom and the p-value of each term of the partition.

}
\item{simulaout}{ When the input parameter \code{simulation = TRUE}, the output includes the object \code{simulaout} which returns \code{nboot} 
number of randomly generated three-way contingency tables. The output also includes the row, column and tube hypothesized probabilities 
\code{pi}, \code{pj}, and \code{pk}, and their observed marginal frequencies defined by the object name \code{margI}, \code{margJ} and \code{margK}, respectively.
Furthermore, the output includes the empirical distribution of each term of the partition of 
the $C_M$-statistic, chi-squared statistic and $C^{S}_{M}$-statistic based on the \code{nboots} randomly generated contingency tables.

    
When \code{simulation = FALSE}, then  \code{simulaout = NULL}.

}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373. \cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
Lombardo R Carlier A D'Ambra L (1996)  Nonsymmetric correspondence analysis for three-way contingency tables. Methodologica,  4, 59-80.
}
\author{
Lombardo R, Takane Y and Beh EJ}
\note{
This function recalls internally many other functions, depending on the setting of the input parameter \code{indexype}. 
It recalls one of the four functions which does a partition under two different Scenarios.
These two Scenarios depend on the theoretical probabilities: 1) the theoretical  probabilities can be prescribed by the analysy. By default, when \code{scen = 1}, 
they are set all equal (homogeneity margins); 2) when \code{scen = 2}, the theoretical probabilities are estimated from the data.
After performing a partition, it gives the output object necessary for printing the results. The print  function is  \code{print.Chi2for3way} . 
This function belongs to the class \code{chi3class}. 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
## The function is currently defined as
data(olive)
chi2x3way(olive, scen = 2, indextype = "tauM", simulation = FALSE, nboots = 100, nran = 1000,
pi = rep(1/dim(olive)[[1]],dim(olive)[[1]]), pj = rep(1/dim(olive)[[2]],dim(olive)[[2]]),
pk = rep(1/dim(olive)[[3]],dim(olive)[[3]]), digits = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{partition}% use one of  RShowDoc("KEYWORDS")
\keyword{chi-squared index}% __ONLY ONE__ keyword per line
\keyword{Marcotorchino index}