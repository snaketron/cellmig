# cellmig: Uncertainty-aware quantitative analysis of high-throughput live cell migration data

# Background

High-throughput tracking of cells with time-lapse microscopy followed by
the acquisition of images at ﬁxed time intervals facilitates the
analysis of cell migration across many wells treated under different
biological conditions. These workflows generate considerable technical
noise and biological variability, and therefore technical and biological
replicates are necessary, leading to large, hierarchically structured
datasets, i.e., cells are nested within technical replicates that are
nested within biological replicates.

Current statistical analyses of such data usually ignore the
hierarchical structure of the data and fail to explicitly quantify
uncertainty arising from technical or biological variability. To address
this gap, we present **cellmig**, an R package implementing
Bayesian hierarchical models for migration analysis.
**cellmig** quantifies condition-specific velocity changes
(e.g., drug effects) while modeling nested data structures and technical
artifacts, providing uncertainty-aware estimates through credible
intervals.

Furthermore, **cellmig** includes functionality for simulating
synthetic datasets. These simulations are invaluable for experimental
planning, allowing researchers to assess how different experimental
designs (e.g., varying numbers of biological replicates, technical
replicates, or cells per well) affect the precision of treatment effect
estimates.


# Installation
To install this package, start R (version "4.5") and enter:

```{r, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cellmig")
```

Case studies are provided in the directory /vignettes


# How to cite
>Kitanovski S, Conroy S, Sonneck J, Claas L, Dorsch M, Urban S, et al. (2026) 
>Uncertainty-aware quantitative analysis of high-throughput live cell migration data. 
>PLoS Comput Biol 22(7): e1014472. https://doi.org/10.1371/journal.pcbi.1014472
