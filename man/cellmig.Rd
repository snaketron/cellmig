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
\item \code{fit} - the fitted model as an \code{rstan} object.
\item \code{data} - the input/processed input data
\item \code{posteriors} - a summary of model parameters: means, medians, and 
95\% credible intervals.
\item \code{control} - the list of input controls
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
\item compound = character. Compound name (c1, c2, c3, etc.)
\item dose = character. Treatment dose (0, 1, 5, 10, low, mid, high, etc.)
\item v = numeric. Observed cell migration speed.
\item offset = binary (0 or 1). Indicates whether a treatment should be used 
for batch correction across plates. Default: 0 (no correction). Set to 1 for 
specific treatment groups used as offsets, ensuring they appear on each plate.
}
}
\item{control}{A list configuring the MCMC sampling algorithm, and parameter
priors with the following default values:
\itemize{
\item mcmc_warmup = 500 (Number of warmup iterations).
\item mcmc_steps = 1500 (Number of sampling iterations).
\item mcmc_chains = 4 (Number of Markov chains).
\item mcmc_cores = 1 (Number of CPU cores used for sampling).
\item mcmc_algorithm = "NUTS" (MCMC algorithm).
\item adapt_delta = 0.8 (Target acceptance probability for NUTS sampler).
\item max_treedepth = 10 (Maximum depth of tree exploration in NUTS sampler).
\item prior_alpha_p_M = 2, prior_alpha_p_SD = 1 (prior on alpha_p: normal
distribution mean (M) and standard deviation (SD))
\item prior_sigma_bio_M = 0, prior_sigma_bio_SD = 0.5 (prior on sigma_bio: 
normal distribution mean (M) and standard deviation (SD))
\item prior_sigma_tech_M = 0, prior_sigma_tech_SD = 0.5 (prior on sigma_tech: 
normal distribution mean (M) and standard deviation (SD))
\item prior_kappa_mu_M = 2, prior_kappa_mu_SD = 1 (prior on kappa_mu: 
normal distribution mean (M) and standard deviation (SD))
\item prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 1 (prior on kappa_sigma: 
normal distribution mean (M) and standard deviation (SD))
\item prior_mu_group_M = 0, prior_mu_group_SD = 1 (prior on mu_group: 
normal distribution mean (M) and standard deviation (SD))
}
}
}
\value{
A list containing:
\itemize{
\item \code{fit} = The fitted model as an rstan object.
\item \code{data} = The raw and processed input data.
\item \code{posteriors} = A summary of model parameters, including means 
and 95\% credible intervals.
\item \code{control} = The list of input controls
}
}

\examples{
data(d, package = "cellmig")
o <- cellmig(x = d,
             control = list(mcmc_warmup = 300,
                    mcmc_steps = 1000,
                    mcmc_chains = 2,
                    mcmc_cores = 2,
                    mcmc_algorithm = "NUTS",
                    adapt_delta = 0.8,
                    max_treedepth = 10))
str(o)
}
