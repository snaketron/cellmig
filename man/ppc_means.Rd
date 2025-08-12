\name{get_ppc_means}
\alias{get_ppc_means}
\title{Observed vs. predicted cell migration means in each well}

\usage{
get_ppc_means(x)
}

\arguments{
  \item{x}{The function takes as its main input (\code{x}) the output of the 
  \code{cellmig} function, containing posterior samples from a Bayesian model. 
  The posterior predictive means are stored in \code{x$s$yhat}, while the 
  observed migration values are extracted from \code{x$x$d}.}
}

\value{
A \code{ggplot2} object displaying:
\itemize{
  \item A scatter plot comparing observed vs. predicted migration means.
  \item 95\% highest density interval (HDI) error bars for the predicted means.
  \item A dashed reference line (\code{y = x}) for ideal agreement.
}
}

\description{
The \code{get_ppc_means} function visualizes the comparison between observed and 
posterior predicted migration means. The function aggregates observed migration 
values by well, merges them with posterior predictive means, and generates a 
scatter plot with error bars representing the 95\% HDI.
}

\details{
The function follows these steps:
\enumerate{
  \item Compute mean observed (scaled) cell migration in each well
  \item Compute mean predicted (scaled) cell migration in each well + 95\% credible 
  intervals
  \item Generate a \code{ggplot2} visualization with:
    \itemize{
      \item A dashed diagonal indicating perfect agreement.
      \item Scatter points representing observed vs. predicted means.
      \item Error bars (dark gray) showing the 95\% HDI for predictions.
    }
}
}

\examples{
# Load example dataset
data(d, package = "cellmig")

# Run cellmig with MCMC sampling
o <- cellmig(x = d,
             control = list(mcmc_warmup = 300,
                            mcmc_steps = 1000,
                            mcmc_chains = 2,
                            mcmc_cores = 2,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))

# Generate PPC means plot
p <- get_ppc_means(x = o)
print(p)
}
