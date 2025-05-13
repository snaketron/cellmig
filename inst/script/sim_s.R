# Simulate data with severe replicate effects
control <- list(N_biorep = 4, 
                N_techrep = 3, 
                N_cell = 20, 
                delta = c(+0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0,
                          -0.8, -0.4, -0.2, +0.0, +0.2, +0.4, +0.8,
                          -0.4, -0.2, -0.1, +0.0, +0.1, +0.2, +0.4,
                          -0.2, -0.1, -0.1, +0.0, +0.2, +0.4, +0.5,
                          +0.4, +0.2, +0.1, +0.0, -0.1, -0.2, -0.4),
                sigma_bio = 0.1, 
                sigma_tech = 0.05, 
                offset = 1,
                prior_alpha_p_M = 1.7,
                prior_alpha_p_SD = 0.5,
                prior_kappa_mu_M = 1.7,
                prior_kappa_mu_SD = 0.5,
                prior_kappa_sigma_M = 0,
                prior_kappa_sigma_SD = 0.3)

y <- gen_partial(control = control)

d <- y$y
d <- d[,c("well", "plate", "compound", "dose", "v")]

d$compound <- paste0("C", rep(x = 1:5, each = 7*20*4*3))
d$dose <- paste0("D", rep(rep(x = 1:7, each = 20*4*3), times = 5))

d$offset <- 0
d$offset[d$compound=="C1"] <- 1

save(d, file = "data/d.RData", compress = TRUE)
