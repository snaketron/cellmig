
gen_partial <- function(control = list(N_biorep = 3, 
                                       N_techrep = 3, 
                                       N_cell = 50, 
                                       delta = c(0, -0.4, -0.2, -0.1, 
                                                 0, 0.1, 0.2, 0.4),
                                       sigma_bio = 0.2, 
                                       sigma_tech = 0.05, 
                                       offset = 1,
                                       prior_alpha_p_M = 1.7,
                                       prior_alpha_p_SD = 0.5,
                                       prior_kappa_mu_M = 1.7,
                                       prior_kappa_mu_SD = 0.5,
                                       prior_kappa_sigma_M = 0,
                                       prior_kappa_sigma_SD = 0.3)) {
    
    get_meta <- function(control) {
        wi = 1
        m <- c()
        for(g in 1:control$N_group) {
            for(p in 1:control$N_plate) {
                for(w in 1:control$N_well_reps) {
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
    
    get_control <- function(control_in) {
        control <- list(N_biorep = 3, 
                        N_techrep = 3, 
                        N_cell = 50, 
                        delta = c(0, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4),
                        sigma_bio = 0.2, 
                        sigma_tech = 0.05, 
                        offset = 1,
                        prior_alpha_p_M = 1.7,
                        prior_alpha_p_SD = 0.5,
                        prior_kappa_mu_M = 1.7,
                        prior_kappa_mu_SD = 0.5,
                        prior_kappa_sigma_M = 0,
                        prior_kappa_sigma_SD = 0.3)
        
        # if missing control_in -> use default values
        if(missing(control_in) || is.null(control_in)) {
            return(control)
        }
        if(is.list(control_in) == FALSE) {
            stop("control must be a list")
        }
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
        
        return(control)
    }
    
    get_sim <- function(control) {
        control$N_group <- length(control$delta)
        control$N_well <- control$N_biorep*control$N_group*control$N_techrep
        yhat <- vector(mode = "list", length = control$N_well)
        kappa <- numeric(length = control$N_well)
        mu <- numeric(length = control$N_well)
        mu_well <- numeric(length = control$N_well)
        
        kappa_mu <- rnorm(n = 1, 
                          mean = control$prior_kappa_mu_M, 
                          sd = control$prior_kappa_mu_SD)
        
        kappa_sigma <- abs(rnorm(n = 1, 
                                 mean = control$prior_kappa_sigma_M, 
                                 sd = control$prior_kappa_sigma_SD))
        
        alpha_p <- rnorm(n = control$N_biorep,
                         mean = control$prior_alpha_p_M,
                         sd = control$prior_alpha_p_SD)
        
        mu_plate_group <- matrix(data = NA, 
                                 nrow = control$N_group, 
                                 ncol = control$N_biorep)
        
        meta <- c()
        well_id <- 1
        for(g in 1:control$N_group) {
            mu_plate_group[g,] <- rnorm(n = control$N_biorep, 
                                        mean = control$delta[g], 
                                        sd = control$sigma_bio)
            for(p in 1:control$N_biorep) {
                for(w in 1:control$N_techrep) {
                    if(g==control$offset) {
                        mu_well[well_id] <- rnorm(n = 1, 
                                                  mean = alpha_p[p], 
                                                  sd = control$sigma_tech)
                    }
                    if(g!=control$offset) {
                        mu_well[well_id] <- rnorm(n = 1,
                                                  mean = alpha_p[p] + mu_plate_group[g,p], 
                                                  sd = control$sigma_tech)
                    }
                    kappa[well_id] = exp(rnorm(n = 1, 
                                               mean = kappa_mu, 
                                               sd = kappa_sigma))
                    
                    mu[well_id] <- exp(mu_well[well_id])
                    y <- rgamma(n = control$N_cell,
                                shape = kappa[well_id], 
                                rate = kappa[well_id]/mu[well_id])
                    
                    yhat[[well_id]] <- data.frame(v = y, 
                                                  well = well_id, 
                                                  plate = p, 
                                                  group = g,
                                                  compound = g, 
                                                  dose = "X")
                    
                    meta <- rbind(meta, data.frame(w = well_id, g = g, p = p))
                    well_id = well_id + 1
                }
            }
        }
        yhat <- do.call(rbind, yhat)
        
        pars <- list(alpha_p = alpha_p,
                     mu_plate_group = mu_plate_group,
                     mu = mu,
                     mu_well = mu_well,
                     kappa = kappa)
        
        return(list(yhat = yhat, pars = pars, meta = meta))
    }
    
    yhat <- get_sim(control = control)
    
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
    
    get_meta <- function(control) {
        wi = 1
        m <- c()
        for(g in 1:control$N_group) {
            for(p in 1:control$N_plate) {
                for(w in 1:control$N_well_reps) {
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
    
    get_control <- function(control_in) {
        control <- list(N_biorep = 3, 
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
                        prior_mu_group_SD = 1.0)
        
        # if missing control_in -> use default values
        if(missing(control_in) || is.null(control_in)) {
            return(control)
        }
        if(is.list(control_in) == FALSE) {
            stop("control must be a list")
        }
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
        
        return(control)
    }
    
    control <- get_control(control_in = control)
    
    # get meta
    meta <- get_meta(control = control)
    
    # sample
    s <- sampling(object = stanmodels$gen_F,
                  algorithm = "Fixed_param", 
                  chains = 1, 
                  iter = control$N_cell+10, 
                  warmup = 10, 
                  data = control,
                  refresh = -1)
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

