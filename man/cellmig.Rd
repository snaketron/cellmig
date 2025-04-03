\name{cellmig}
\alias{cellmig}
\title{Model-based quantification of cell migration speed}
\description{
The \code{cellmig} function estimates cell migration speed using a 
Bayesian model implemented in \code{rstan}. It takes a data.frame as 
input and allows users to configure the Markov Chain Monte Carlo (MCMC) 
sampling procedure via the \code{control} parameter. The function 
returns a list containing:

\itemize{
\item \code{f} - the fitted model as an \code{rstan} object.
\item \code{x} - the processed input data.
\item \code{s} - a summary of model parameters, including means, medians, and 
95\% credible intervals.
}
}
\usage{
cellmig(x, control = NULL)
}
\arguments{
\item{x}{A data.frame containing the following columns. Each row represents
the features of a cell:
\itemize{
\item well = character. Well ID (w1, w2, w3, etc.).
\item plate = character. plate ID (p1, p2, p3, etc.). A plate contains 
multiple wells.
\item experiment = character. experiment ID (e1, e2, e3, etc.). An experiment 
typically contains multiple plates.
\item compound = character. Compound name (c1, c2, c3, etc.)
\item dose = character. Treatment dose (0, 1, 5, 10, low, mid, high, etc.)
\item v = numeric. Observed cell migration speed.
\item offset = binary (0 or 1). Indicates whether a treatment should be used 
for batch correction across plates. Default: 0 (no correction). Set to 1 for 
specific treatment groups used as offsets, ensuring they appear on each plate.
}
}
\item{control}{A list configuring the MCMC sampling algorithm with the 
following default values:
\itemize{
\item mcmc_warmup = 500 (Number of warmup iterations).
\item mcmc_steps = 1500 (Number of sampling iterations).
\item mcmc_chains = 4 (Number of Markov chains).
\item mcmc_cores = 1 (Number of CPU cores used for sampling).
\item mcmc_algorithm = "NUTS" (MCMC algorithm).
\item adapt_delta = 0.8 (Target acceptance probability for NUTS sampler).
\item max_treedepth = 10 (Maximum depth of tree exploration in NUTS sampler).
}
}
}
\value{
A list containing:
\itemize{
\item \code{f} = The fitted model as an rstan object.
\item \code{x} = The processed input data.
\item \code{s} = A summary of model parameters, including means and 95\% 
credible intervals.
}
}

\examples{
data(d, package = "cellmig")
#cellmig_out <- cellmig(x = d,
#             control = list(mcmc_warmup = 300,
#                            mcmc_steps = 1500,
#                            mcmc_chains = 3,
#                            mcmc_cores = 3,
#                            mcmc_algorithm = "NUTS",
#                            adapt_delta = 0.8,
#                            max_treedepth = 10))
#
data(cellmig_out, package = "cellmig")
str(cellmig_out)
}
