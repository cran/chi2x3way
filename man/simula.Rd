\name{simula}
\alias{simula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulations of three-way contingency tables
}
\description{
Allows the generation of random contingency tables with prescribed theoretical probabilities
}
\usage{
simula(I, J, K, nran = 1000, pi, pj, pk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{I}{
The input parameter for specifying the number of rows.}
  \item{J}{
The input parameter for specifying the number of columns.}
  \item{K}{
The input parameter for specifying the number of tubes.}
  \item{nran}{
The input parameter for specifying the number of individuals in each table.}
  \item{pi}{
The input parameter for specifying the theoretical probability  of rows.}
  \item{pj}{
The input parameter for specifying the theoretical probability  of columns.}
  \item{pk}{
The input parameter for specifying the theoretical probability  of tubes.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Lancaster H O (1951) Complex contingency tables treated by the partition of the chi-square. Journal of Royal Statistical Society, Series B, 13, 242-249. \cr
Loisel S and Takane Y (2016) Partitions of Pearson's chi-square ststistic for frequency tables:  A comprehensive account. Computational Statistics, 31, 1429-1452.\cr
}
\author{
Lombardo R, Takane Y and Beh EJ
}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
simula(4,3,2,pi=rep(1/4,4),pj=rep(1/3,3),pk=rep(1/2,2))
## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{contingency table}% use one of  RShowDoc("KEYWORDS")
\keyword{simulation}% __ONLY ONE__ keyword per line
