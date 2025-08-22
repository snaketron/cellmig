\name{d}
\alias{d}
\alias{d_mini}
\docType{data}
\title{
Example dataset \code{d}
}
\description{
This dataset \code{d} contains simulated cell migration speed data from an 
imaginery experiment. It includes migration speed values (column \code{v}) 
for different samples (column \code{sample}) treated with chemical compounds 
(column \code{compound}) administered at a dost (column \code{dose}) on
experimental plate with identifier (column \code{plate}). The dataset 
\code{d_mini} contains a subset of \code{d} with only four of the compounds 
(C1 to C4) and three doses (D3, D5, and D5).
}
\usage{
data("d", package = "cellmig")
}
\format{
  A data frame with 6,300 (cells) observations with the following features.
  \describe{
    \item{\code{v}}{a numeric vector with cell migration speeds}
    \item{\code{well}}{a character vector}
    \item{\code{compound}}{a character vector}
    \item{\code{dose}}{a character vector}
    \item{\code{plate}}{a character vector}
    \item{\code{offset}}{indicator: 1 = control, 0 = non-control. Used for
    batch correction}
  }
}
\details{
The code used to generate this data is in  inst/script/sim_s.R
}
\source{
The code used to generate this data is in  inst/script/sim_s.R
}
\references{
The code used to generate this data is in  inst/script/sim_s.R
}
\examples{
data(d, package = "cellmig")
data(d_mini, package = "cellmig")
str(d)
str(d_mini)
}
\keyword{datasets}
