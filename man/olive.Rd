\name{olive}
\alias{olive}
\docType{data}
\title{
Agresti's data (1990), three-way contingency table
}
\description{
The data consists of 6 rows, 3 columns and 2 tubes. The rows represent
the preference for black olives of Armed Forces personnel, the columns are concerned
with the location and the tubes with the  urbanization.
}
\usage{data(olive)}
\format{
  The format is:
 row names [1:6] "A" "B" "C" "D" "E" "F" \cr
 col names [1:3] "NW" "NE" "SW"
 tube names [1:2]  "urban" "rural"
}
\references{
Agresti A (1990) Categorical Data Analysis  (pg. 335). John Wiley & Sons.
}
\examples{
data(olive)
dim(olive)
dimnames(olive)
}
\keyword{datasets}
