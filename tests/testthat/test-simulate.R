###### testing simulations #######

test_that("get_partial with invalid control", {
  expect_error(gen_partial(),
               regexp = "object 'delta' not found")
  expect_error(gen_partial(control = list(N_biorep = 3, 
                                          N_techrep = 3, 
                                          N_cell = 50, 
                                          delta,
                                          sigma_bio = 0.2, 
                                          sigma_tech = 0.05, 
                                          offset = 1,
                                          prior_alpha_p_M = 1.7,
                                          prior_alpha_p_SD = 0.5,
                                          prior_kappa_mu_M = 1.7,
                                          prior_kappa_mu_SD = 0.5,
                                          prior_kappa_sigma_M = 0,
                                          prior_kappa_sigma_SD = 0.3)),
               regexp = "object 'delta' not found")
  expect_error(gen_partial(control = list(N_biorep = 3, 
                                          N_techrep = 3, 
                                          N_cell = 50, 
                                          delta = 0,
                                          sigma_bio = 0.2, 
                                          sigma_tech = 0.05, 
                                          offset = 1,
                                          prior_alpha_p_M = 1.7,
                                          prior_alpha_p_SD = 0.5,
                                          prior_kappa_mu_M = 1.7,
                                          prior_kappa_mu_SD = 0.5,
                                          prior_kappa_sigma_M = 0,
                                          prior_kappa_sigma_SD = 0.3)),
               regexp = "delta must have length > 1")
  expect_no_error(gen_partial(control = list(N_biorep = 3, 
                                             N_techrep = 3, 
                                             N_cell = 50, 
                                             delta = c(0,2),
                                             sigma_bio = 0.2, 
                                             sigma_tech = 0.05, 
                                             offset = 1,
                                             prior_alpha_p_M = 1.7,
                                             prior_alpha_p_SD = 0.5,
                                             prior_kappa_mu_M = 1.7,
                                             prior_kappa_mu_SD = 0.5,
                                             prior_kappa_sigma_M = 0,
                                             prior_kappa_sigma_SD = 0.3)))
  
  expect_error(gen_partial(control = list(N_biorep = 3, 
                                          N_techrep = 3, 
                                          N_cell = 50, 
                                          delta = c(0,2),
                                          sigma_bio = 0.2, 
                                          sigma_tech = 0.05, 
                                          offset = 3,
                                          prior_alpha_p_M = 1.7,
                                          prior_alpha_p_SD = 0.5,
                                          prior_kappa_mu_M = 1.7,
                                          prior_kappa_mu_SD = 0.5,
                                          prior_kappa_sigma_M = 0,
                                          prior_kappa_sigma_SD = 0.3)),
               regexp = "offset must be between 1 and length of delta")
  
  expect_error(gen_partial(control = NA),
               regexp = "control must be a list")
  expect_error(gen_partial(control = NULL),
               regexp = "missing control")
  expect_error(gen_partial(control = NA),
               regexp = "control must be a list")
})



# 
# test_that("get_partial with invalid control", {
#   expect_error(gen_fpartial(),
#                regexp = "object 'delta' not found")
#   expect_error(gen_partial(control = list(N_biorep = 3, 
#                                           N_techrep = 3, 
#                                           N_cell = 50, 
#                                           delta,
#                                           sigma_bio = 0.2, 
#                                           sigma_tech = 0.05, 
#                                           offset = 1,
#                                           prior_alpha_p_M = 1.7,
#                                           prior_alpha_p_SD = 0.5,
#                                           prior_kappa_mu_M = 1.7,
#                                           prior_kappa_mu_SD = 0.5,
#                                           prior_kappa_sigma_M = 0,
#                                           prior_kappa_sigma_SD = 0.3)),
#                regexp = "object 'delta' not found")
#   expect_error(gen_partial(control = list(N_biorep = 3, 
#                                           N_techrep = 3, 
#                                           N_cell = 50, 
#                                           delta = 0,
#                                           sigma_bio = 0.2, 
#                                           sigma_tech = 0.05, 
#                                           offset = 1,
#                                           prior_alpha_p_M = 1.7,
#                                           prior_alpha_p_SD = 0.5,
#                                           prior_kappa_mu_M = 1.7,
#                                           prior_kappa_mu_SD = 0.5,
#                                           prior_kappa_sigma_M = 0,
#                                           prior_kappa_sigma_SD = 0.3)),
#                regexp = "delta must have length > 1")
#   expect_no_error(gen_partial(control = list(N_biorep = 3, 
#                                              N_techrep = 3, 
#                                              N_cell = 50, 
#                                              delta = c(0,2),
#                                              sigma_bio = 0.2, 
#                                              sigma_tech = 0.05, 
#                                              offset = 1,
#                                              prior_alpha_p_M = 1.7,
#                                              prior_alpha_p_SD = 0.5,
#                                              prior_kappa_mu_M = 1.7,
#                                              prior_kappa_mu_SD = 0.5,
#                                              prior_kappa_sigma_M = 0,
#                                              prior_kappa_sigma_SD = 0.3)))
#   expect_error(gen_partial(
#     control = list(N_biorep = 3, 
#                    N_techrep = 3, 
#                    N_cell = 50, 
#                    delta = c(0,2),
#                    sigma_bio = 0.2, 
#                    sigma_tech = 0.05, 
#                    offset = 3,
#                    prior_alpha_p_M = 1.7,
#                    prior_alpha_p_SD = 0.5,
#                    prior_kappa_mu_M = 1.7,
#                    prior_kappa_mu_SD = 0.5,
#                    prior_kappa_sigma_M = 0,
#                    prior_kappa_sigma_SD = 0.3)),
#     regexp = "offset must be an integer in the range 1 to length\\(deltas\\)")
#   
#   expect_error(gen_partial(control = NA),
#                regexp = "control must be a list")
#   expect_error(gen_partial(control = NULL),
#                regexp = "missing control")
#   expect_error(gen_partial(control = NA),
#                regexp = "control must be a list")
# })
