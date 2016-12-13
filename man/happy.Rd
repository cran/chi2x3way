\name{happy}
\alias{happy}
\docType{data}
\title{
Three-way contingency table by Davis (1977)
}
\description{
The data consists of 3 rows, 5 columns and 4 tubes. The rows represent
the levels of Happiness, the columns are concerned
with the number of siblings and the tubes with the  Education.
}
\usage{data(olive)}
\format{
  The format is:
 row names [1:6] "H1" "H2" "H3"  \cr
 col names [1:3] "S1" "S2" "S3" "S4" "S5"
 tube names [1:2]  "E1" "E2" "E3" "E4"
}
\references{
Davis JA (1977) Codebook for the 1977 General Social Survey. Chicago: National Opinion Research Centre. 
}
\examples{
data(happy)
dim(happy)
dimnames(happy)
}
\keyword{datasets}
