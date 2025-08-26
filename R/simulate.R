
gen_partial <- function(control = list(N_biorep = 3, 
                                       N_techrep = 3, 
                                       N_cell = 50, 
                                       delta,
                                       sigma_bio = 0.1, 
                                       sigma_tech = 0.05, 
                                       offset = 1,
                                       prior_alpha_p_M = 1.7,
                                       prior_alpha_p_SD = 0.5,
                                       prior_kappa_mu_M = 1.7,
                                       prior_kappa_mu_SD = 0.5,
                                       prior_kappa_sigma_M = 0,
                                       prior_kappa_sigma_SD = 0.3)) {
    
    check_sim_control(control = control, partial = TRUE)
    control <- get_sim_control(control_in = control, partial = TRUE)
    yhat <- get_partial_sim_data(control = control)
    
    return(list(y = yhat$yhat, 
                par = yhat$pars, 
                meta = yhat$meta, 
                control = control))
}

gen_full <- function(control = list(N_biorep = 3, 
                                    N_techrep = 3, 
                                    N_cell = 50, 
                                    N_group = 5,
                                    prior_alpha_p_M = 1.7,
                                    prior_alpha_p_SD = 1.0,
                                    prior_kappa_mu_M = 1.7,
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
    check_control_list(control_in = control, partial)
    
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
        control <- list(N_biorep = 3, N_techrep = 3, N_cell = 50, 
                        delta = c(0, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4),
                        sigma_bio = 0.2, sigma_tech = 0.05, offset = 1,
                        prior_alpha_p_M = 1.7, prior_alpha_p_SD = 0.5,
                        prior_kappa_mu_M = 1.7, prior_kappa_mu_SD = 0.5,
                        prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 0.3)
        if(all(names(control_in) %in% names(control)) == FALSE) {
            stop("unrecognized elements found in control")
        }
        ns <- names(control_in)
        for (i in seq_len(length(control_in))) {
            control[[ns[i]]] <- control_in[[ns[i]]]
        }
        control$N_group <- length(control$delta)
        control$N_plate <- control$N_biorep
        control$N_well_reps <- control$N_techrep
        control$mu_group <- control$delta
        control$prior_sigma_bio_M <- 0.0
        control$prior_sigma_bio_SD <- 0.0
        control$prior_sigma_tech_M <- 0.0
        control$prior_sigma_tech_SD <- 0.0
        control$prior_mu_group_M <- 0.0
        control$prior_mu_group_SD <- 0.0
    } else {
        control <- list(N_biorep = 3, N_techrep = 3, N_cell = 50, N_group = 5,
                        prior_alpha_p_M = 1.7, prior_alpha_p_SD = 1.0,
                        prior_kappa_mu_M = 1.7, prior_kappa_mu_SD = 1.0,
                        prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 1.0,
                        prior_sigma_bio_M = 0.0, prior_sigma_bio_SD = 1.0,
                        prior_sigma_tech_M = 0.0, prior_sigma_tech_SD = 1.0,
                        prior_mu_group_M = 0.0, prior_mu_group_SD = 1.0)
        if(all(names(control_in) %in% names(control)) == FALSE) {
            stop("unrecognized elements found in control")
        }
        ns <- names(control_in)
        for (i in seq_len(length(control_in))) {
            control[[ns[i]]] <- control_in[[ns[i]]]
        }
        control$N_plate <- control$N_biorep
        control$N_well_reps <- control$N_techrep
        control$offset <- 1
    }
    return(control)
}

get_partial_sim_data <- function(control) {
    control$N_group <- length(control$delta)
    control$N_well <- control$N_biorep*control$N_group*control$N_techrep
    yhat <- vector(mode = "list", length = control$N_well)
    kappa <- numeric(length = control$N_well)
    mu <- numeric(length = control$N_well)
    mu_well <- numeric(length = control$N_well)
    kappa_mu <- rnorm(n = 1, mean = control$prior_kappa_mu_M, 
                      sd = control$prior_kappa_mu_SD)
    kappa_sigma <- abs(rnorm(n = 1, mean = control$prior_kappa_sigma_M, 
                             sd = control$prior_kappa_sigma_SD))
    alpha_p <- rnorm(n = control$N_biorep, mean = control$prior_alpha_p_M,
                     sd = control$prior_alpha_p_SD)
    mu_plate_group <- matrix(data = NA, nrow = control$N_group, 
                             ncol = control$N_biorep)
    meta <- c()
    well_id <- 1
    for(g in seq_len(control$N_group)) {
        mu_plate_group[g,] <- rnorm(n = control$N_biorep, mean=control$delta[g],
                                    sd = control$sigma_bio)
        for(p in seq_len(control$N_biorep)) {
            for(w in seq_len(control$N_techrep)) {
                if(g==control$offset) {
                    mu_well[well_id] <- rnorm(n = 1, mean = alpha_p[p], 
                                              sd = control$sigma_tech)
                }
                if(g!=control$offset) {
                    mu_well[well_id] <- rnorm(n = 1,
                                              mean = alpha_p[p] + 
                                                  mu_plate_group[g,p], 
                                              sd = control$sigma_tech)
                }
                kappa[well_id] <- exp(rnorm(n = 1, mean = kappa_mu, 
                                            sd = kappa_sigma))
                
                mu[well_id] <- exp(mu_well[well_id])
                y <- rgamma(n = control$N_cell, shape = kappa[well_id], 
                            rate = kappa[well_id]/mu[well_id])
                yhat[[well_id]] <- data.frame(v = y, well = well_id, plate = p, 
                                              group = g, compound = g, dose="X")
                meta <- rbind(meta, data.frame(w = well_id, g = g, p = p))
                well_id <- well_id + 1
            }
        }
    }
    pars <- list(alpha_p = alpha_p, mu_plate_group = mu_plate_group, mu = mu,
                 mu_well = mu_well, kappa = kappa)
    return(list(yhat = do.call(rbind, yhat), pars = pars, meta = meta))
}

get_meta_control <- function(control) {
    wi <- 1
    m <- c()
    for(g in seq_len(control$N_group)) {
        for(p in seq_len(control$N_plate)) {
            for(w in seq_len(control$N_well_reps)) {
                m <- rbind(m, data.frame(well_id = wi, 
                                         group_id = g, 
                                         plate_id = p,
                                         trep_id = w))
                wi <- wi + 1
            }
        }
    }
    return(m)
}

check_control_list <- function(control_in, partial) {
    if(partial) {
        control <- list(N_biorep = NA, N_techrep = NA, N_cell = NA, 
                        delta = NA, sigma_bio = NA, sigma_tech = NA, 
                        offset = NA, 
                        prior_alpha_p_M = NA, prior_alpha_p_SD = NA, 
                        prior_kappa_mu_M = NA, prior_kappa_mu_SD = NA,
                        prior_kappa_sigma_M = NA, prior_kappa_sigma_SD = NA)
    } else {
        control <- list(N_biorep = 3, N_techrep = 3, N_cell = 50, 
                        N_group = 5,
                        prior_alpha_p_M = 1.7, prior_alpha_p_SD = 1.0,
                        prior_kappa_mu_M = 1.7, prior_kappa_mu_SD = 1.0,
                        prior_kappa_sigma_M = 0, prior_kappa_sigma_SD = 1.0,
                        prior_sigma_bio_M = 0.0, prior_sigma_bio_SD = 1.0,
                        prior_sigma_tech_M = 0.0, prior_sigma_tech_SD = 1.0,
                        prior_mu_group_M = 0.0, prior_mu_group_SD = 1.0)
    }
    
    # if missing control_in -> use default values
    if(missing(control_in) || is.null(control_in)) {
        stop("missing control")
    }
    if(is.list(control_in) == FALSE) {
        stop("control must be a list")
    }
    if(all(names(control) %in% names(control_in))==FALSE) {
        stop("control has missing elements")
    }
}
