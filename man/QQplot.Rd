\name{QQplot}
\alias{QQplot}
\title{
QQ-plots of the simulated distribution of the Marcotorchino's index partition  
}
\description{
It allows \cr
the graphical representation of each term of  the Marcotorchino's index partition under scenario 2 or 1. 
}
\usage{
QQplot(nsample=100, yobs, nameC, taudf)
}
\arguments{
  \item{nsample}{
The number  of random tables to generate. For each table, the terms of  index partition are computed.
}
  \item{yobs}{
The term of the index partition, it represents the observed distribution of the C-statistic associated to the term of the partition.
}
  \item{nameC}{
The label of the term of the index partition.
}
  \item{taudf}{
The number of degree of freedom associated to the term of the index partition.
}
}
%\value{
%}
\references{
Beh EJ and Lombardo R (2014) Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons. \cr
Carlier A Kroonenberg PM (1996) Biplots and decompositions in two-way and three-way correspondence analysis. Psychometrika, 61, 355-373. \cr
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R, Takane Y and Beh EJ}
\note{
This function is called from the function tauMbootQQ and  allows to depict graphically the three-way index distribution. 
A  QQ-plot is produced for each term of the index partition. 
}
%\examples{
%QQplot(nsample = 1000, 
%    digits = 3, pi=rep(1/rows,rows), pj=rep(1/cols,cols), pk=rep(1/tubs,tubs)) 
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}% use one of  RShowDoc("KEYWORDS")
\keyword{ three-way table}% __ONLY ONE__ keyword per line
