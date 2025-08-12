cellmig <- function(x, control = NULL) {
  
  # check inputs
  x <- process_input(x)
  
  # check control
  control <- process_control(control_in = control)
  
  # fit model
  f <- get_fit(x = x, control = control)
  
  # get summary
  s <- get_summary(x = x, f = f)
  
  return(list(fit = f, posteriors = s, data = x, control = control))
}

get_fit <- function(x, control) {
  message("model fitting... \n")
  
  M <- stanmodels$M
  
  pars <- c("z_1", "z_2", "z_3", "mu_well")
  
  # fit model
  fit <- sampling(object=M,
                  data=list(y=x$proc_x$v, 
                            N=x$proc_x$N[1], 
                            N_well=x$proc_x$N_well[1], 
                            N_plate=x$proc_x$N_plate[1],
                            N_plate_group=x$proc_x$N_plate_group[1],
                            N_group=x$proc_x$N_group[1],
                            well_id=x$proc_x$well_id,
                            plate_id=x$map_w$plate_id,
                            plate_group_id=x$map_w$plate_group_id,
                            offset=x$map_w$offset,
                            group_id=x$map_pg$group_id,
                            prior_alpha_p_M=control$prior_alpha_p_M,
                            prior_alpha_p_SD=control$prior_alpha_p_SD,
                            prior_sigma_bio_M=control$prior_sigma_bio_M,
                            prior_sigma_bio_SD=control$prior_sigma_bio_SD,
                            prior_sigma_tech_M=control$prior_sigma_tech_M,
                            prior_sigma_tech_SD=control$prior_sigma_tech_SD,
                            prior_kappa_mu_M=control$prior_kappa_mu_M,
                            prior_kappa_mu_SD=control$prior_kappa_mu_SD,
                            prior_kappa_sigma_M=control$prior_kappa_sigma_M,
                            prior_kappa_sigma_SD=control$prior_kappa_sigma_SD,
                            prior_mu_group_M=control$prior_mu_group_M,
                            prior_mu_group_SD=control$prior_mu_group_SD),
                  pars=pars,
                  include = FALSE,
                  chains=control$mcmc_chains,
                  cores=control$mcmc_cores,
                  iter=control$mcmc_steps,
                  warmup=control$mcmc_warmup,
                  algorithm=control$mcmc_algorithm,
                  control=list(adapt_delta=control$adapt_delta,
                               max_treedepth=control$max_treedepth),
                  refresh=200)
  
  return(fit)
}

get_summary <- function(x, f) {
  message("computing posterior summaries...\n")
  
  x <- x$proc_x
  
  # get meta data
  l <- x[, c("well_id", "well", "group_id", "group", 
             "compound", "dose", "plate_id", "plate", 
             "plate_group_id", "plate_group", "offset")]
  meta_well <- l[duplicated(l)==FALSE, ]
  meta_plate <- l[duplicated(l[, c("plate", "plate_id")])==FALSE, 
                  c("plate", "plate_id")]
  meta_group <- l[duplicated(l[, c("group", "group_id")])==FALSE, 
                  c("group", "group_id", "compound", "dose", 
                    "plate_id", "plate")]
  meta_plate_group <- l[duplicated(l[, c("group", "group_id", 
                                         "plate", "plate_id",
                                         "plate_group", 
                                         "plate_group_id")])==FALSE, 
                        c("well_id", "group", "group_id", 
                          "plate", "plate_id",
                          "compound", "dose", 
                          "plate_group", "plate_group_id")]
  
  # par: alpha_p
  alpha_p <- data.frame(summary(f, par = "alpha_p")$summary)
  alpha_p$plate_id <- seq_len(nrow(alpha_p))
  alpha_p <- merge(x = alpha_p, y = meta_plate, 
                   by = "plate_id", all.x = TRUE)
  
  # par: mu_group
  mu_group <- data.frame(summary(f, par = "mu_group")$summary)
  mu_group$group_id <- seq_len(nrow(mu_group))
  mu_group <- merge(x = mu_group, y = meta_group, by = "group_id", all.x = TRUE)
  
  # par: mu_plate_group
  mu_plate_group <- data.frame(summary(f, par = "mu_plate_group")$summary)
  mu_plate_group$plate_group_id <- seq_len(nrow(mu_plate_group))
  mu_plate_group <- merge(x = mu_plate_group, y = meta_plate_group, 
                          by = "plate_group_id", all.x = TRUE)
  
  # par: mu
  mu_well <- data.frame(summary(f, par = "mu")$summary)
  mu_well$well_id <- seq_len(nrow(mu_well))
  mu_well <- merge(x = mu_well, y = meta_well, by = "well_id", all.x = TRUE)
  
  # par: kappa
  kappa_well <- data.frame(summary(f, par = "kappa")$summary)
  kappa_well$well_id <- seq_len(nrow(kappa_well))
  kappa_well <- merge(x = kappa_well, y = meta_well, 
                      by = "well_id", all.x = TRUE)
  
  # par: sigma_bio
  sigma_bio <- data.frame(summary(f, par = "sigma_bio")$summary)
  sigma_tech <- data.frame(summary(f, par = "sigma_tech")$summary)
  
  # kappa mu and sigma
  kappa_mu <- data.frame(summary(f, par = "kappa_mu")$summary)
  kappa_sigma <- data.frame(summary(f, par = "kappa_sigma")$summary)
  
  # par: y_hat_sample
  yhat <- data.frame(summary(f, par = "y_hat_sample")$summary)
  yhat$well_id <- seq_len(nrow(yhat))
  yhat <- merge(x = yhat, y = meta_well, by = "well_id", all.x = TRUE)
  
  return(list(alpha_p = alpha_p,
              mu_group = mu_group, 
              mu_plate_group = mu_plate_group,
              mu_well = mu_well,
              kappa_well = kappa_well,
              kappa_mu = kappa_mu,
              kappa_sigma = kappa_sigma,
              sigma_bio = sigma_bio,
              sigma_tech = sigma_tech,
              yhat = yhat))
}
