\docType{package}
\name{cellmig-package}
\alias{cellmig-package}
\title{cellmig: quantifying cell migration soeed with hierarchical Bayesian 
models}
\description{
The \code{cellmig} package implements Bayesian hierarchical models for the
analysis of cell migration speed. It quantifies condition-specific effects
(e.g., drug treatments) on cell migration speed while explicitly modeling
nested data structures and technical artifacts, and provides uncertainty-aware
estimates via posterior credible intervals.
}
\details{
This package contains functions for modeling, clustering, and visualization
of cell migration speed data derived from high-throughput migration assays.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/snaketron/cellmig}
  \item Report bugs at \url{https://github.com/snaketron/cellmig/issues}
}
}

\author{
Authors and maintainers:
\itemize{
  \item Simo Kitanovski \email{simokitanovski@gmail.com} 
  (\href{https://orcid.org/0000-0003-2909-5376}{ORCID})
}
}

\examples{
# load package input data
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
