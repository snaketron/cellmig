###### testing x #######

test_that("x parameter takes invalid", {
  expect_error(cellmig(), 
               regexp = "missing x")
  expect_error(cellmig(x = NULL), 
               regexp = "x must be a data.frame")
  expect_error(cellmig(x = NA), 
               regexp = "x must be a data.frame")
  expect_error(cellmig(x = 0), 
               regexp = "x must be a data.frame")
  expect_error(cellmig(x = F), 
               regexp = "x must be a data.frame")
  expect_error(cellmig(x = data.frame()), 
               regexp = "zero or one row in x")
})


test_that("x parameter takes column values", {
  data("d", package = "cellmig")
  # use only one well (one compound)
  expect_error(cellmig(x = d[d$well=="1",]), 
               regexp = "only one compound included")
  
  expect_error(cellmig(x = d[d$compound=="C1",]), 
               regexp = "only one compound included")
  
  data("d", package = "cellmig")
  d$offset[d$offset==1] <- 2
  expect_error(cellmig(x = d),
               regexp = "column offset can contain only 0 or 1")
  
  data("d", package = "cellmig")
  d$offset[d$offset==1] <- NA
  expect_error(cellmig(x = d),
               regexp = "column offset can contain only 0 or 1")
  
  data("d", package = "cellmig")
  d$offset[d$offset==1] <- Inf
  expect_error(cellmig(x = d),
               regexp = "column offset can contain only 0 or 1")
  
  data("d_mini", package = "cellmig")
  d_mini$offset <- 0
  expect_no_error(cellmig(x = d_mini,
                       control = list(mcmc_warmup = 200,
                                      mcmc_steps = 900,
                                      mcmc_chains = 1,
                                      mcmc_cores = 1)))
})


test_that("x parameter takes column values", {
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2"),]
  d <- d[d$plate == "1", ]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  expect_length(object = o, n = 4)
  
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2"),]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  expect_length(object = o, n = 4)
})


