\name{get_ppc_violins}
\alias{get_ppc_violins}
\title{Posterior Predictive Checks (PPC) Visualization}

\usage{
get_ppc_violins(x, wrap = FALSE, ncol = 4)
}

\arguments{
  \item{x}{The function takes as its main input (\code{x}) the output of the 
\code{cellmig} function, containing posterior samples from a Bayesian model, 
specifically the output of a fitted Stan model with posterior predictive 
samples stored in \code{x$f}.}
\item{wrap}{logical, if wrap = FALSE (default) we will use facet_grid, 
and if wrap = TRUE we will use facet_wrap to split data into compound x 
plate pairs.}
\item{ncol}{In case wrap = TRUE, how many columns should the plot have?}
}

\value{
A \code{ggplot2} object displaying:
\itemize{
  \item Violin plots of posterior predictive distributions for each compound 
  and plate.
  \item Overlaid scatter plots (sina plot) of observed values.
}
}

\description{
The \code{get_ppc_violins} function visualizes posterior predictive checks (PPC) by 
comparing observed data with posterior predictive distributions. The function
extracts simulated data from the fitted Stan model and overlays them with 
observed values, grouped by compound and plate.
}

\details{
The function extracts posterior predictive samples using \code{rstan::extract}, 
reshapes them into long format using \code{reshape2::melt}, and merges them 
with metadata from the input object. It then generates a \code{ggplot2} 
visualization, where:
\itemize{
  \item Violin plots (dashed red) represent posterior predictive distributions.
  \item Overlaid sina plots (black) represent observed data points.
  \item Facets are arranged by compound and plate for better comparison.
}
}

\examples{
data(d_mini, package = "cellmig")
o <- cellmig(x = d_mini,
             control = list(mcmc_warmup = 200,
                            mcmc_steps = 500,
                            mcmc_chains = 2,
                            mcmc_cores = 2,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))

p <- get_ppc_violins(x = o, wrap = TRUE, ncol = 4)
print(p)
}
