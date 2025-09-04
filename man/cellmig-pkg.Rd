\name{cellmig-pkg}
\alias{cellmig-pkg}
\title{cellmig: quantifying cell migration with hierarchical Bayesian models}
\description{
\code{cellmig} implements Bayesian hierarchical models for migration analysis
and quantifies condition-specific (drug effects) effects on cell velocity 
while modeling nested data structures and technical artifacts, providing 
uncertainty-aware estimates through credible intervals.
}
\details{
This package contains functions for modeling, clustering, and visualization
of cell migration velocity from high-throughput migration assays.
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
  \item Simo Kitanovski \email{simokitanovski@uni-due.de} 
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
