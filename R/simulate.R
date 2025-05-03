
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
                    wi <- wi + 1;
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
    
    control <- get_control(control_in = control)
    
    # get meta
    meta <- get_meta(control = control)
    
    # sample
    s <- sampling(object = stanmodels$gen_P,
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
                    wi <- wi + 1;
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

