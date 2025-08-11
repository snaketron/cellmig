\name{gen_partial}
\alias{gen_partial}
\title{Simulating data from partially generative Bayesian model}
\description{
  Simulates cell velocities based on a hierarchical model using a Stan model,
  where some of the model parameter values are defined (fixed) by the user,
  while the rest are drawn from the prior distributions.
}
\usage{
  gen_partial(control = list(N_biorep = 3, 
                             N_techrep = 3, 
                             N_cell = 50, 
                             delta = c(0, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4),
                             sigma_bio = 0.2, 
                             sigma_tech = 0.05, 
                             offset = 1,
                             prior_alpha_p_M = 1.7,
                             prior_alpha_p_SD = 0.5,
                             prior_kappa_mu_M = 1.7,
                             prior_kappa_mu_SD = 0.5,
                             prior_kappa_sigma_M = 0,
                             prior_kappa_sigma_SD = 0.3))
}
\arguments{
  \item{control}{A list configuring the MCMC sampling algorithm, and parameter
priors with the following default values:
\itemize{
  \item{N_biorep: Number of plates (biological replicates) in the simulation.}
  \item{N_techrep: Number of technical replicates: wells on a plate treated 
  with the same treatment.}
  \item{N_cell: Number of cells per well.}
  \item{delta: Effects (on log-scale) on cell velocity for each treatment.}
  \item{sigma_bio: Variability between biological replicates.}
  \item{sigma_tech: Variability between technical replicates.}
  \item{offset: Index of the control treatment (which treatment should be used
  for batch correction between plates).}
  \item{prior_alpha_p_M, prior_alpha_p_SD: Mean (M) and standard deviation 
  (SD) of a normal distribution describing plate-specific batch effects 
  (on log-scale).}
  \item{prior_kappa_mu_M, prior_kappa_mu_SD: Mean (M) and standard deviation 
  (SD) of a normal distribution describing the mean of the population of 
  well-specific shape (kappa) parameters.}
  \item{prior_kappa_sigma_M, prior_kappa_sigma_SD: Mean (M) and standard 
  deviation (SD) of a normal distribution describing the standard deviation 
  of the population of well-specific shape (kappa) parameters.}
}
}
}
\details{
  This function constructs a hierarchical dataset by simulating values 
  using a partially generative Stan model. It generates metadata for each 
  well and simulates data using the `sampling` function from the `rstan` 
  package.
}
\value{
  A data frame containing:
  \itemize{
    \item{iteration = Simulation iteration index.}
    \item{well_id = Unique well identifier.}
    \item{y = Simulated response value.}
    \item{group_id = Experimental group identifier.}
    \item{plate_id = Plate identifier.}
  }
}
\examples{
g <- gen_partial(control = list(N_biorep = 3, 
                                N_techrep = 3, 
                                N_cell = 50, 
                                delta=c(0, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4),
                                sigma_bio = 0.2, 
                                sigma_tech = 0.05, 
                                offset = 1,
                                prior_alpha_p_M = 1.7,
                                prior_alpha_p_SD = 0.5,
                                prior_kappa_mu_M = 1.7,
                                prior_kappa_mu_SD = 0.5,
                                prior_kappa_sigma_M = 0,
                                prior_kappa_sigma_SD = 0.3))
str(g)
}
