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
\item v = numeric. Observed cell migration speed (in micrometers/min).
\item offset = binary (0 or 1). Indicates whether a treatment should be used 
as control for batch correction across plates. Default: 0 (no correction). Set 
to 1 for specific treatment groups used as offsets, ensuring they appear on 
each plate.
}
}
\item{control}{A list configuring MCMC sampling settings and parameter
priors with the following default values:
\itemize{
  \item \strong{MCMC controls:}
  \itemize{
  \item \code{mcmc_warmup = 500}: number of warmup iterations.
  \item \code{mcmc_steps = 1500}: number of sampling iterations.
  \item \code{mcmc_chains = 4}: number of Markov chains.
  \item \code{mcmc_cores = 1}: number of CPU cores used for sampling.
  \item \code{mcmc_algorithm = "NUTS"}: MCMC algorithm.
  \item \code{adapt_delta = 0.8}: target acceptance probability for the NUTS sampler.
  \item \code{max_treedepth = 10}: maximum tree depth for the NUTS sampler.
  }
  \item \strong{Priors:}
  \itemize{
  \item \code{prior_alpha_p_M = -0.5}, \code{prior_alpha_p_SD = 1}:
  mean and standard deviation of a normal distribution prior for 
  \code{alpha_p}, representing the plate-specific mean (log-scale) cell speed.
  \item \code{prior_sigma_bio_M = 0}, \code{prior_sigma_bio_SD = 1}:
  mean and standard deviation of a half-normal distribution prior for
  \code{sigma_bio}, representing variability between biological replicates.
  \item \code{prior_sigma_tech_M = 0}, \code{prior_sigma_tech_SD = 1}:
  mean and standard deviation of a half-normal distribution prior for 
  \code{sigma_tech}, representing variability between technical replicates.
  \item \code{prior_kappa_mu_M = 1.5}, \code{prior_kappa_mu_SD = 1}:
  mean and standard deviation of a normal distribution prior for 
  \code{kappa_mu}, the mean of the population distribution of Gamma shape 
  parameters.
  \item \code{prior_kappa_sigma_M = 0}, \code{prior_kappa_sigma_SD = 1}:
  mean and standard deviation of a half-normal distribution prior for 
  \code{kappa_sigma}, the standard deviation of the population distribution 
  of Gamma shape parameters.
  \item \code{prior_sigma_delta_M = 0}, 
  \code{sigma_delta_SD = 1}:
  mean and standard deviation of a half-normal distribution prior for
  \code{sigma_delta}, representing the standard deviation of the population 
  of \code{delta_t} parameters.
  }
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
\itemize{
    \item{\code{alpha_p}: batch effect on plate p}
    \item{\code{delta_t}: overall treatment effects relative to the selected 
    control}
    \item{\code{delta_pt}: plate-specific treatment effects relative to the 
    selected 
    control}
    \item{\code{well_mu}: mean of cell speed distribution per well}
    \item{\code{well_kappa}: shape of cell speed distribution per well}
    \item{\code{kappa_mu}, \code{kappa_sigma}: mean/standard deviation 
    parameters of the population of well_kappas}
    \item{\code{sigma_bio}: variability between biological replicates}
    \item{\code{sigma_tech}: variability between technical replicates}
    \item{\code{sigma_delta}: variability between \code{delta_t} parameters}
    \item{\code{yhat}: posterior predictions}
    }
\item \code{control} = The list of input controls
}
}
\examples{
data(d, package = "cellmig")
o <- cellmig(x = d,
             control = list(mcmc_warmup = 200,
                            mcmc_steps = 700,
                            mcmc_chains = 2,
                            mcmc_cores = 2,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))
str(o)
}
