###### testing downstream #######

test_that("downstream functions with invalid inputs", {
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2"),]
  d <- d[d$plate == "1", ]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  
  groups <- get_groups(x = o)
  expect_no_error(get_groups(x = o))
  expect_error(get_groups())
  expect_error(get_groups(x = NULL))
  
  expect_no_error(get_violins(x = o, exponentiate = TRUE, 
                              from_group = "C2|D1",to_group = "C2|D2"))
  expect_no_error(get_violins(x = o, exponentiate = FALSE, 
                              from_group = "C2|D1",to_group = "C2|D2"))
  expect_error(get_violins(x = o, from_groups = "C2|D1",to_group = "C2|D2"))
  expect_error(get_violins(x = o, from_groups = "C2|D1",to_group = "C2|D2"))
  
  
  expect_warning(get_pairs(x = o, exponentiate = TRUE))
  expect_warning(get_pairs(x = o, exponentiate = FALSE))
  
  expect_error(get_pairs())
  expect_error(get_pairs(x = o))
  expect_error(get_pairs(x = o, groups = 1))
  expect_error(get_pairs(x = o, groups = NA))
  expect_error(get_pairs(x = o, groups = NULL))
  expect_no_error(get_pairs(x = o, groups = groups$group, exponentiate = TRUE))
  expect_no_error(get_pairs(x = o, groups = groups$group, exponentiate = FALSE))
})
