\name{gen_full}
\alias{gen_full}
\title{Simulating data from fully generative Bayesian model}
\description{
  Simulates cell velocities based on a hierarchical model using a Stan model,
  where some of the model parameters are drawn from their prior distributions.
}
\usage{
  gen_full(control = list(N_biorep = 3, 
                          N_techrep = 3, 
                          N_cell = 50, 
                          prior_alpha_p_M = 1.7,
                          prior_alpha_p_SD = 1.0,
                          prior_kappa_mu_M = 1.7,
                          prior_kappa_mu_SD = 1.0,
                          prior_kappa_sigma_M = 0,
                          prior_kappa_sigma_SD = 1.0,
                          prior_sigma_bio_M = 0.0,
                          prior_sigma_bio_SD = 1.0,
                          prior_sigma_tech_M = 0.0,
                          prior_sigma_tech_SD = 1.0,
                          prior_mu_group_M = 0.0,
                          prior_mu_group_SD = 1.0))
}
\arguments{
  \item{N_biorep}{Number of plates (biological replicates) in the simulation.}
  \item{N_techrep}{Number of technical replicates: wells on a plate treated 
  with the same treatment.}
  \item{N_cell}{Number of cells per well.}
  \item{sigma_bio}{Variability between biological replicates.}
  \item{sigma_tech}{Variability between technical replicates.}
  \item{prior_alpha_p_M, prior_alpha_p_SD}{Mean (M) and standard deviation 
  (SD) of a normal distribution describing plate-specific batch effects 
  (on log-scale).}
  \item{prior_kappa_mu_M, prior_kappa_mu_SD}{Mean (M) and standard deviation 
  (SD) of a normal distribution describing the mean of the population of 
  well-specific shape (kappa) parameters.}
  \item{prior_kappa_sigma_M, prior_kappa_sigma_SD}{Mean (M) and standard 
  deviation (SD) of a normal distribution describing the standard deviation 
  of the population of well-specific shape (kappa) parameters.}
}
\details{
  This function constructs a hierarchical dataset by simulating values 
  using a partially generative Stan model. It generates metadata for each 
  well and simulates data using the `sampling` function from the `rstan` 
  package.
}
\value{
  A data.frame containing cells as rows with:
  \itemize{
    \item{v = Simulated cell velocity.}
    \item{well_id = Unique well identifier.}
    \item{group_id = Experimental group identifier.}
    \item{plate_id = Biological replicate (plate) identifier.}
    \item{trep_id = Technical replicate identifier.}
  }
}
\examples{
f <- gen_full(control = list(N_biorep = 3, 
                             N_techrep = 3, 
                             N_cell = 50, 
                             prior_alpha_p_M = 1.7,
                             prior_alpha_p_SD = 1.0,
                             prior_kappa_mu_M = 1.7,
                             prior_kappa_mu_SD = 1.0,
                             prior_kappa_sigma_M = 0,
                             prior_kappa_sigma_SD = 1.0,
                             prior_sigma_bio_M = 0.0,
                             prior_sigma_bio_SD = 1.0,
                             prior_sigma_tech_M = 0.0,
                             prior_sigma_tech_SD = 1.0,
                             prior_mu_group_M = 0.0,
                             prior_mu_group_SD = 1.0))
str(f)
}
