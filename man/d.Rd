\name{d}
\alias{d}
\alias{d_mini}
\docType{data}
\title{
Example dataset \code{d}
}
\description{
The dataset \code{d} contains simulated cell migration speed data from an
imaginary experiment. It includes cell migration speed measurements
(column \code{v}) for individual cells measured in different wells,
treated with chemical compounds (\code{compound}) at varying doses
(\code{dose}) and measured on experimental plates (\code{plate}).

The dataset \code{d_mini} is a reduced subset of \code{d}, containing data
for four compounds (C1--C4) and three doses.
}
\usage{
data("d", package = "cellmig")
data(d_mini, package = "cellmig")
}
\format{
  A data frame with cells (rows) having with the following features.
  \describe{
    \item{\code{v}}{Numeric. Cell migration speed (micrometers per minute).}
    \item{\code{well}}{Character. Well identifier.}
    \item{\code{compound}}{Character. Compound identifier.}
    \item{\code{dose}}{Character. Treatment dose.}
    \item{\code{plate}}{Character. Plate identifier.}
    \item{\code{offset}}{Integer (0 or 1). Indicates control samples used for 
    batch correction across plates (1 = control, 0 = non-control).}
  }
}
\details{
Both datasets were generated using a simulation script located at
\code{inst/scripts/sim_s.R}.
}
\source{
Simulated data generated using \code{inst/scripts/sim_s.R}.
}
\references{
Simulated data generated using \code{inst/scripts/sim_s.R}.
}
\examples{
data(d, package = "cellmig")
data(d_mini, package = "cellmig")
str(d)
str(d_mini)
}
\keyword{datasets}
