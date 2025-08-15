###### testing others #######

test_that("test ppc.R error", {
  expect_error(get_ppc_means(), 
               regexp = "missing x")
  expect_error(get_ppc_means(x = NULL), 
               regexp = "missing x")
  
  expect_error(get_ppc_violins(), 
               regexp = "missing x")
  expect_error(get_ppc_violins(x = NULL), 
               regexp = "missing x")
})


test_that("test ppc.R correct", {
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2"),]
  d <- d[d$plate == "1", ]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  
  expect_no_error(get_ppc_means(x = o))
  expect_no_error(get_ppc_violins(x = o, wrap = FALSE))
  
  expect_error(get_ppc_violins(x = o, wrap = TRUE),
               regexp = "missing ncol")
  expect_error(get_ppc_violins(x = o, wrap = NULL), 
               regexp = "missing wrap")
  expect_error(get_ppc_violins(x = o, ncol = NULL),
               regexp = "missing wrap")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = -1),
               regexp = "ncol must be positive integer")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = TRUE),
               regexp = "ncol must be positive integer")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = NA),
               regexp = "ncol must be positive integer")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = NULL),
               regexp = "missing ncol")
})

