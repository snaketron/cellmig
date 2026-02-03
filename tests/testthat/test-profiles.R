###### testing profiles #######

test_that("profiles functions with invalid inputs", {
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2"),]
  d <- d[d$plate == "1", ]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  groups <- get_groups(x = o)$group
  expect_error(get_dose_response_profile(x = o, groups = groups,
                                         exponentiate = TRUE),
               regexp = "must have n >= 2 objects to cluster")
  expect_error(get_dose_response_profile(), regexp = "x is missing or NULL")
  expect_error(get_dose_response_profile(x = o, 
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = NA,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = NULL,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = character(length = 1),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = c("average", "complete"),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = NA,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = NULL,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = character(length = 1),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = c("average", "complete"),
                                         groups = groups,
                                         exponentiate = TRUE))
})



###### testing profiles #######

test_that("profiles functions with invalid inputs", {
  data("d", package = "cellmig")
  d <- d[d$compound %in% c("C1", "C2", "C3"),]
  d <- d[d$plate == "1", ]
  o <- suppressWarnings(cellmig(x = d, control = list(mcmc_chains = 1,
                                                      mcmc_warmup = 200,
                                                      mcmc_steps = 500)))
  groups <- get_groups(x = o)$group
  
  
  expect_warning(get_dose_response_profile(x = o, 
                                            groups = groups,
                                            exponentiate = TRUE))
  expect_error(get_dose_response_profile(), regexp = "x is missing or NULL")
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = NA,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = NULL,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = character(length = 1),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_link = c("average", "complete"),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = NA,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = NULL,
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = character(length = 1),
                                         groups = groups,
                                         exponentiate = TRUE))
  expect_error(get_dose_response_profile(x = o, 
                                         hc_dist = c("average", "complete"),
                                         groups = groups,
                                         exponentiate = TRUE))
})

