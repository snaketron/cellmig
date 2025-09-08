###### testing others #######

test_that("test ppc.R error", {
  expect_error(get_ppc_means(), 
               regexp = "x is missing or NULL")
  expect_error(get_ppc_means(x = NULL), 
               regexp = "x is missing or NULL")
  expect_error(get_ppc_violins(), 
               regexp = "x is missing or NULL")
  expect_error(get_ppc_violins(x = NULL), 
               regexp = "x is missing or NULL")
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
  
  expect_no_error(get_ppc_violins(x = o, wrap = TRUE))
  expect_error(get_ppc_violins(x = o, wrap = NULL), 
               regexp = "wrap is missing or NULL")
  expect_no_error(get_ppc_violins(x = o, ncol = NULL))
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = -1),
               regexp = "ncol must be >0")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = TRUE),
               regexp = "ncol must be numeric")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = NA),
               regexp = "ncol must be numeric")
  expect_error(get_ppc_violins(x = o, wrap = TRUE, ncol = NULL),
               regexp = "ncol is missing or NULL")
  
  
  expect_no_error(get_ppc_violins(x = o, wrap = TRUE, ncol = 1))
})

