\name{gen_partial}
\alias{gen_partial}
\title{Simulate data from a partially generative hierarchical Bayesian model}
\description{
Simulates cell migration speed data from a hierarchical Bayesian model
implemented in Stan, where some model parameters are fixed by the user
and the remaining parameters are drawn from their prior distributions.
}
\usage{
  gen_partial(control = list(N_biorep = 3, 
                             N_techrep = 3, 
                             N_cell = 50, 
                             delta,
                             sigma_bio = 0.1, 
                             sigma_tech = 0.05, 
                             offset = 1,
                             prior_alpha_p_M = -0.5,
                             prior_alpha_p_SD = 1.0,
                             prior_kappa_mu_M = 1.5,
                             prior_kappa_mu_SD = 1.0,
                             prior_kappa_sigma_M = 0,
                             prior_kappa_sigma_SD = 0.3))
}
\arguments{
\item{control}{
A named list specifying simulation settings, fixed parameter values, and
prior distributions. Default values are:
\itemize{
  \item \code{N_biorep}: Number of biological replicates (plates).
  \item \code{N_techrep}: Number of technical replicates (wells per plate
  receiving the same treatment).
  \item \code{N_cell}: Number of cells per well.
  \item \code{delta}: Numeric vector of treatment effects on cell migration
  speed (log-scale), one value per treatment group.
  \item \code{sigma_bio}: Variability in treatment effects between biological
  replicates.
  \item \code{sigma_tech}: Variability in treatment effects between technical
  replicates.
  \item \code{offset}: Index of the control treatment used for batch
  correction across plates.
  \item \code{prior_alpha_p_M}, \code{prior_alpha_p_SD}:
  Mean and standard deviation of the normal prior for plate-specific
  batch effects (log-scale).
  \item \code{prior_kappa_mu_M}, \code{prior_kappa_mu_SD}:
  Mean and standard deviation of the normal prior for the population mean 
  (log-scale) of well-specific Gamma shape parameters.
  \item \code{prior_kappa_sigma_M}, \code{prior_kappa_sigma_SD}:
  Mean and standard deviation of the normal prior for the population
  standard deviation of well-specific Gamma shape parameters.
}
}
}
\details{
The function constructs a hierarchical synthetic dataset by simulating
cell migration speed values from a partially generative Stan model.
User-specified parameters (e.g., treatment effects and variability terms)
are held fixed, while remaining parameters are sampled from their prior
distributions. Cell-level observations are simulated using the
\code{sampling} function from the \pkg{rstan} package.
}
\value{
A data frame with one row per simulated observation and the following
columns:
\describe{
  \item{\code{iteration}}{Simulation iteration index.}
  \item{\code{well_id}}{Unique well identifier.}
  \item{\code{y}}{Simulated cell migration speed.}
  \item{\code{group_id}}{Treatment group identifier.}
  \item{\code{plate_id}}{Biological replicate (plate) identifier.}
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
                                prior_alpha_p_M = -0.5,
                                prior_alpha_p_SD = 0.5,
                                prior_kappa_mu_M = 1.7,
                                prior_kappa_mu_SD = 0.5,
                                prior_kappa_sigma_M = 0,
                                prior_kappa_sigma_SD = 0.3))
str(g)
}
