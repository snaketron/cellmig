\name{gen_full}
\alias{gen_full}
\title{Simulate data from a fully generative hierarchical Bayesian model}
\description{
Simulates cell migration speed data from a hierarchical Bayesian model
implemented in Stan, where model parameters are drawn from their prior
distributions.
}
\usage{
  gen_full(control = list(N_biorep = 3, 
                          N_techrep = 3,
                          N_cell = 50,
                          N_group = 5,
                          prior_alpha_p_M = -0.5,
                          prior_alpha_p_SD = 1,
                          prior_kappa_mu_M = 1.5,
                          prior_kappa_mu_SD = 1,
                          prior_kappa_sigma_M = 0,
                          prior_kappa_sigma_SD = 1,
                          prior_sigma_bio_M = 0,
                          prior_sigma_bio_SD = 1,
                          prior_sigma_tech_M = 0,
                          prior_sigma_tech_SD = 1,
                          prior_sigma_delta_M = 0,
                          prior_sigma_delta_SD = 1))
}
\arguments{
\item{control}{
A named list specifying simulation settings and prior distributions.
Default values are:
\itemize{
  \item \code{N_biorep}: Number of biological replicates (plates).
  \item \code{N_techrep}: Number of technical replicates (wells per plate
  receiving the same treatment).
  \item \code{N_cell}: Number of cells per well.
  \item \code{N_group}: Number of treatment groups.
  \item \code{prior_alpha_p_M}, \code{prior_alpha_p_SD}:
  Mean and standard deviation of the normal prior for plate-specific
  batch effects (log-scale).
  \item \code{prior_kappa_mu_M}, \code{prior_kappa_mu_SD}:
  Mean and standard deviation of the normal prior for the population mean
  of well-specific Gamma shape parameters.
  \item \code{prior_kappa_sigma_M}, \code{prior_kappa_sigma_SD}:
  Mean and standard deviation of the half-normal prior for the population
  standard deviation of well-specific Gamma shape parameters.
  \item \code{prior_sigma_bio_M}, \code{prior_sigma_bio_SD}:
  Mean and standard deviation of the half-normal prior describing variability
  in treatment effects between biological replicates.
  \item \code{prior_sigma_tech_M}, \code{prior_sigma_tech_SD}:
  Mean and standard deviation of the half-normal prior describing variability
  in treatment effects between technical replicates.
  \item \code{prior_sigma_delta_M}, \code{prior_sigma_delta_SD}:
  Mean and standard deviation of the half-normal prior describing variability
  in treatment effects.
}
}
}
\details{
The function constructs a hierarchical synthetic dataset by simulating
cell migration speed values from a fully generative Stan model. Metadata
for plates, wells, treatment groups, and replicates are generated, and
cell-level observations are simulated using the \code{sampling} function
from the \pkg{rstan} package.
}
\value{
A data frame with one row per cell and the following columns:
\describe{
  \item{\code{v}}{Simulated cell migration speed.}
  \item{\code{well_id}}{Unique well identifier.}
  \item{\code{group_id}}{Treatment group identifier.}
  \item{\code{plate_id}}{Biological replicate (plate) identifier.}
  \item{\code{trep_id}}{Technical replicate identifier.}
}
}
\examples{
f <- gen_full(control = list(N_biorep = 3, 
                             N_techrep = 3, 
                             N_cell = 50,
                             N_group = 5,
                             prior_alpha_p_M = -0.5,
                             prior_alpha_p_SD = 1.0,
                             prior_kappa_mu_M = 1.5,
                             prior_kappa_mu_SD = 1.0,
                             prior_kappa_sigma_M = 0,
                             prior_kappa_sigma_SD = 1.0,
                             prior_sigma_bio_M = 0.0,
                             prior_sigma_bio_SD = 1.0,
                             prior_sigma_tech_M = 0.0,
                             prior_sigma_tech_SD = 1.0,
                             prior_sigma_delta_M = 0.0,
                             prior_sigma_delta_SD = 1.0))
str(f)
}
