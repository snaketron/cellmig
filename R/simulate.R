
gen_partial <- function(control = list(N_biorep = 3, 
                                       N_techrep = 3, 
                                       N_cell = 50, 
                                       delta,
                                       sigma_bio = 0.1, 
                                       sigma_tech = 0.05, 
                                       offset = 1,
                                       prior_alpha_p_M = -0.5,
                                       prior_alpha_p_SD = 1.0,
                                       prior_kappa_mu_M = 1.5,
                                       prior_kappa_mu_SD = 1.0,
                                       prior_kappa_sigma_M = 0,
                                       prior_kappa_sigma_SD = 0.3)) {
    
    check_sim_control(control = control, partial = TRUE)
    control <- get_sim_control(control_in = control, partial = TRUE)
    o <- get_partial_sim_data(control = control)
    return(list(y = o$yhat, par = o$pars, meta = o$meta, control = control))
}

gen_full <- function(control = list(N_biorep = 3, 
                                    N_techrep = 3, 
                                    N_cell = 50, 
                                    N_group = 5,
                                    prior_alpha_p_M = -0.5,
                                    prior_alpha_p_SD = 1.0,
                                    prior_kappa_mu_M = 1.5,
                                    prior_kappa_mu_SD = 1.0,
                                    prior_kappa_sigma_M = 0,
                                    prior_kappa_sigma_SD = 1.0,
                                    prior_sigma_bio_M = 0.0,
                                    prior_sigma_bio_SD = 1.0,
                                    prior_sigma_tech_M = 0.0,
                                    prior_sigma_tech_SD = 1.0,
                                    prior_mu_group_M = 0.0,
                                    prior_mu_group_SD = 1.0)) {
    
    check_sim_control(control = control, partial = FALSE)
    control <- get_sim_control(control_in = control, partial = FALSE)
    # get meta
    meta <- get_meta_control(control = control)
    
    # sample
    s <- sampling(object = stanmodels$gen_F, algorithm = "Fixed_param", 
                  chains = 1, iter = control$N_cell+10, warmup = 10, 
                  data = control, refresh = -1)
    # extract
    y <- extract(object = s, par = "y")$y
    
    # get data.frame
    y <- melt(data = y)
    colnames(y) <- c("i", "well_id", "v")
    y <- merge(x = y, y = meta, by = "well_id", all.x = TRUE)
    y <- y[y$i <= control$N_cell, ]
    
    y$well <- as.character(y$well_id)
    y$plate <- as.character(y$plate_id)
    y$group <- as.character(y$group_id)
    y$compound <- y$group
    y$dose <- "X"
    y <- y[,c("i", "v", "well", "plate", "group", 
              "compound", "dose", "trep_id")]
    
    return(list(y = y, control = control))
}

check_sim_control <- function(control, partial) {
    check_control_list(control_in = control, partial = partial)
    
    check_positive_integer(y = control$N_biorep, par = "N_biorep")
    check_positive_integer(y = control$N_techrep, par = "N_techrep")
    check_positive_integer(y = control$N_cell, par = "N_cell")
    
    check_sd(y = control$prior_alpha_p_SD, par = "prior_alpha_p_SD")
    check_sd(y = control$prior_kappa_mu_SD, par = "prior_kappa_mu_SD")
    check_sd(y = control$prior_kappa_sigma_SD, par="prior_kappa_sigma_SD")
    
    check_loc(y = control$prior_alpha_p_M, par = "prior_alpha_p_M")
    check_loc(y = control$prior_kappa_mu_M, par = "prior_kappa_mu_M")
    check_sd(y = control$prior_kappa_sigma_M, par = "prior_kappa_sigma_M")
    
    if(partial) {
        check_sd(y = control$sigma_bio, par = "sigma_bio")
        check_sd(y = control$sigma_tech, par = "sigma_tech")
        check_delta(x = control$delta, o = control$offset)
    } 
    else {
        check_sd(y = control$prior_sigma_bio_M, par = "prior_sigma_bio_M")
        check_sd(y = control$prior_sigma_bio_SD, par = "prior_sigma_bio_SD")
        check_sd(y = control$prior_sigma_tech_M, par = "prior_sigma_tech_M")
        check_sd(y = control$prior_sigma_tech_SD, par = "prior_sigma_tech_SD")
        check_sd(y = control$prior_mu_group_SD, par = "prior_mu_group_SD")
        check_sd(y = control$prior_mu_group_M, par = "prior_mu_group_M")
        check_positive_integer(y = control$N_group, par = "N_group")
    }
}

get_sim_control <- function(control_in, partial) {
    if(partial) {
        control <- get_default_control_partial_sim()
        if(all(names(control_in) %in% names(control)) == FALSE) {
            stop("unrecognized elements found in control")
        }
        control[names(control_in)] <- control_in
        control$N_group <- length(control$delta)
        control$N_plate <- control$N_biorep
        control$N_well_reps <- control$N_techrep
        control$mu_group <- control$delta
    } 
    else {
        control <- get_default_control_full_sim()
        if(all(names(control_in) %in% names(control)) == FALSE) {
            stop("unrecognized elements found in control")
        }
        control[names(control_in)] <- control_in
        control$N_plate <- control$N_biorep
        control$N_well_reps <- control$N_techrep
    }
    return(control)
}

get_partial_sim_data <- function(control) {
    control$N_group <- length(control$delta)
    control$N_well <- control$N_biorep*control$N_group*control$N_techrep
    
    well_grid <- expand.grid(g = seq_len(control$N_group),
                             p = seq_len(control$N_biorep),
                             w = seq_len(control$N_techrep))
    mu_plate_group <- t(vapply(X = seq_len(control$N_group), 
                               FUN = function(g) {
                                   return(rnorm(n = control$N_biorep, 
                                                mean = control$delta[g], 
                                                sd = control$sigma_bio))
                               }, FUN.VALUE = numeric(control$N_biorep)))
    
    kappa_mu <- rnorm(n = 1, mean = control$prior_kappa_mu_M, 
                      sd = control$prior_kappa_mu_SD)
    kappa_sigma <- abs(rnorm(n = 1, mean = control$prior_kappa_sigma_M, 
                             sd = control$prior_kappa_sigma_SD))
    alpha_p <- rnorm(n = control$N_biorep, mean = control$prior_alpha_p_M, 
                     sd = control$prior_alpha_p_SD)
    
    
    # Use mapply to compute for each well
    results <- mapply(function(g, p, w, well_id) {
        if(g == control$offset) {
            mu_well <- rnorm(n = 1, mean = alpha_p[p], sd = control$sigma_tech)
        } else {
            mu_well <- rnorm(n = 1, mean = alpha_p[p] + mu_plate_group[g, p], 
                             sd = control$sigma_tech)
        }
        kappa <- exp(rnorm(n = 1, mean = kappa_mu, sd = kappa_sigma))
        mu <- exp(mu_well)
        y <- rgamma(n = control$N_cell, shape = kappa, rate = kappa/mu)
        list(yhat = data.frame(v = y, well = well_id, plate = p, group = g, 
                               compound = g, dose = "X"),
             meta = data.frame(w = well_id, g = g, p = p),
             mu_well = mu_well,
             kappa = kappa,
             mu = mu)
    }, well_grid$g, well_grid$p, well_grid$w, seq_len(nrow(well_grid)), 
    SIMPLIFY = FALSE)
    
    yhat <- do.call(rbind, lapply(results, `[[`, "yhat"))
    meta <- do.call(rbind, lapply(results, `[[`, "meta"))
    mu_well <- vapply(X = results, `[[`, "mu_well", FUN.VALUE = numeric(1))
    kappa <- vapply(X = results, `[[`, "kappa", FUN.VALUE = numeric(1))
    mu <- vapply(X = results, `[[`, "mu", FUN.VALUE = numeric(1))
    
    pars <- list(alpha_p = alpha_p, mu_plate_group = mu_plate_group, 
                 mu = mu, mu_well = mu_well, kappa = kappa)
    return(list(yhat = yhat, pars = pars, meta = meta))
}

get_meta_control <- function(control) {
    meta <- expand.grid(group_id = seq_len(control$N_group),
                        plate_id = seq_len(control$N_plate),
                        trep_id  = seq_len(control$N_well_reps))
    meta$well_id <- seq_len(nrow(meta))
    meta <- meta[, c("well_id", "group_id", "plate_id", "trep_id")]
    return(meta)
}

check_control_list <- function(control_in, partial) {
    control <- get_default_control_full_sim()
    if(partial) {
        control <- get_default_control_partial_sim()
    }
    
    # if missing control_in -> use default values
    check_missing(y = control_in, par = "control")
    check_list(y = control_in, par = "control")
    if(all(names(control_in) %in% names(control))==FALSE) {
        stop("control has missing elements")
    }
}

get_default_control_partial_sim <- function() {
    control <- list(N_biorep = 3, N_techrep = 3, N_cell = 50, 
                    delta = c(0, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4),
                    sigma_bio = 0.2, sigma_tech = 0.05, offset = 1,
                    prior_alpha_p_M = -0.5, prior_alpha_p_SD = 1.0,
                    prior_kappa_mu_M = 1.5, prior_kappa_mu_SD = 1.0,
                    prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 0.3,
                    prior_sigma_bio_M = 0.0,
                    prior_sigma_bio_SD = 0.0,
                    prior_sigma_tech_M = 0.0,
                    prior_sigma_tech_SD = 0.0,
                    prior_mu_group_M = 0.0,
                    prior_mu_group_SD = 0.0)
    return(control)
}

get_default_control_full_sim <- function() {
    control <- list(N_biorep = 3, N_techrep = 3, N_cell = 50, 
                    N_group = 5,
                    prior_alpha_p_M = -0.5, prior_alpha_p_SD = 1.0,
                    prior_kappa_mu_M = 1.5, prior_kappa_mu_SD = 1.0,
                    prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 1.0,
                    prior_sigma_bio_M = 0.0, prior_sigma_bio_SD = 1.0,
                    prior_sigma_tech_M = 0.0, prior_sigma_tech_SD = 1.0,
                    prior_mu_group_M = 0.0, prior_mu_group_SD = 1.0,
                    offset = 1)
    return(control)
}

