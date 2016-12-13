\name{r3dtable}
\alias{r3dtable}
\title{
Simulations for generating three-way contingency tables 
}
\description{
It allows 
1) the generation of \code{nboots=1000}  randomly tables where the 
row, column, tube probabilities can be prescribed by the analyst. 
By default, they are uniform. 
}
\usage{
r3dtable(I = 3, J = 3, K = 3, pi=rep(1/I,I), pj=rep(1/J,J), pk=rep(1/K,K), nboots = 1000, 
nran = 10000, digits = 3)
}
\arguments{
\item{I}{The number   \code{I} is set equal to the rows of the input table \code{X}.
}
\item{J}{The number   \code{J} is set equal to the columns of the input table \code{X}.
}
\item{K}{The number   \code{K} is set equal to the tubes of the input three-way table \code{X}.
}
\item{pi}{
The prescribed row probabilities. By default, they are homogeneous.
}
\item{pj}{The prescribed column probabilities.  By default, they are homogeneous.
}
\item{pk}{The prescribed tube probabilities.  By default, they are homogeneous.
}
\item{nboots}{The number of the random three-way tables  that you want to generate.
}
\item{nran}{The total number of individuals of each generated three-way table.
}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
}
\value{
\item{XB}{The \code{nboots=1000} randomly generated three-way tables.}
\item{XB[[i]]$pi }{The row, prescribed probabilities of the i.th randomly generated three-way table.}
 \item{XB[[i]]$pj }{The column, prescribed probabilities of the i.th randomly generated three-way table.}
 \item{XB[[i]]$pk}{The tube, prescribed probabilities of the i.th randomly generated three-way table.}
\item{margI }{The row observed margins of the randomly generated three-way table.}
 \item{margJ }{The column observed margins of the randomly generated three-way table.}
 \item{margK }{The tube observed margins of the randomly generated three-way table.}
 }
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R, Takane Y, Beh EJ}
\note{
This function allows the generation of random tables under the complete independence with different  theoretical probabilities. 
}
\examples{
r3dtable(I = 3, J = 3, K = 3, pi=rep(1/3,3), pj=rep(1/3,3), pk=rep(1/3,3), 
nboots = 10, nran = 1000, digits = 3) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}% use one of  RShowDoc("KEYWORDS")
\keyword{three-way table}% __ONLY ONE__ keyword per line
\keyword{random tables}