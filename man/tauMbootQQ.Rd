\name{tauMbootQQ}
\alias{tauMbootQQ}
\title{
Simulations for observing the Marcotorchino's index distribution 
}
\description{
It allows 
1) the generation of \code{nboots=1000}  randomly tables where the 
row, column, tube probabilities are set equal to the observed margins of the three-way table considered for the partition under Scenario 2. 
While under Scenario 1, the row, column, tube probabilities are prescribed by the analyst. By default, they are homogeneous. 
}
\usage{
tauMbootQQ(rows = 3, cols = 3, tubs = 3, nboots = 1000, nran = 10000, digits = 3, 
scen=2, pi, pj, pk)
}
\arguments{
\item{rows}{The number  of \code{rows} is set equal to the rows of the input table \code{X}.
}
\item{cols}{The number  of \code{cols} is set equal to the columns of the input table \code{X}.
}
\item{tubs}{The number  of \code{tubs} is set equal to the tubes of the input three-way table \code{X}.
}
\item{nboots}{The number of three-way tables randomly generated.
}
\item{nran}{The total number of individuals of each generated three-way table.
}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
\item{scen}{The input parameter for specifying the Scenario under which the theoretical probabilities are computed. Under Scenario 1 the probabilities 
are prescribed by the analyst, by default they are set homogeneous.}
\item{pi}{
The prescribed row probabilities. By default they are equal to the row margins of the input three-way table \code{X}.
}
\item{pj}{The prescribed column probabilities. By default they are equal to the column margins of the input three-way table \code{X}.
}
\item{pk}{The prescribed tube probabilities. By default they are equal to the tube margins of the input three-way table \code{X}.
}
}
\value{
\item{XG}{The \code{nboots=1000} randomly generated three-way tables.}
 \item{margI }{The row observed margins of the randomly generated three-way table.}
 \item{margJ }{The column observed margins of the randomly generated three-way table.}
 \item{margK }{The tube observed margins of the randomly generated three-way table.}
 \item{ytau}{The table of the terms of the Marcotorchino's index  and of the $C_M$-statistic partition, associated to each of the randomly generated three-way table.}
 \item{ytauNew}{The table of the new expression of the terms of the $C_M$-statistic partition, associated to each of the randomly generated three-way table.}
 \item{ychi}{The table of the terms of the chi-square index partition associated to each of the randomly generated three-way table.}
\item{chidf}{The table of the degree of freedom related to each terms of the chi-square index partition of the randomly generated three-way table.}
\item{cont}{The number of  the randomly generated three-way table whose margin products is less than 5.}
}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373. \cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R, Takane Y and Beh EJ}
\note{
This function allows the generation of random tables under different  theoretical probabilities. 
It  allow to depict graphically the three-way index distribution. From calling the function \code{QQplot}, a  QQ-plot is produced for each term of the partition of
three indices: the Pearson's index, the classic $C_M$-statistic and of the revised $C_M$-statistic. 
}
\examples{
tauMbootQQ(rows = 3, cols = 3, tubs = 3, nboots = 10, nran = 1000, digits = 3, 
pi=rep(1/3,3), pj=rep(1/3,3), pk=rep(1/3,3)) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}% use one of  RShowDoc("KEYWORDS")
\keyword{three-way table}% __ONLY ONE__ keyword per line
\keyword{random tables}