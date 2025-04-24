
sim <- function(N_well_cells,
                N_plate, 
                N_group, 
                N_well_reps, 
                kappa_mu,
                kappa_sigma,
                sigma_bio, 
                sigma_tech, 
                alpha_p,
                mu_group,
                offset) {
    
    data_in <- list(N_plate = N_plate,
                    N_group = N_group,
                    N_well_reps = N_well_reps,
                    N_well_cells = N_well_cells,
                    kappa_mu = kappa_mu,
                    kappa_sigma = kappa_sigma,
                    sigma_bio = sigma_bio,
                    sigma_tech = sigma_tech,
                    alpha_p = alpha_p,
                    mu_group = mu_group,
                    offset = offset)
    
    message("simulation... \n")
    
    y <- sim_d(data_in = data_in)
    
    return(y)
}




sim_d <- function(data_in) {
    
    get_gamma <- function(x, w, n_cells) {
        y <- rgamma(n = n_cells, shape = w$kappa[x], 
                    rate = w$kappa[x]/w$mu[x])
        y <- data.frame(v = y,
                        well = paste0("w", w$well_id[x]),
                        plate = paste0("p", w$p[x]), 
                        group = paste0("g", w$g[x]), 
                        offset = w$offset[x])
        return(y)
    }
    
    m <- matrix(data = 0, nrow = data_in$N_group, ncol = data_in$N_plate)
    for(g in 1:data_in$N_group) {
        m[g,] <- rnorm(n = data_in$N_plate, mean = data_in$mu_group[g],
                       sd = data_in$sigma_bio)
        # mu_raw <- rt(n = data_in$N_plate, df = data_in$mu_df)
        # sd_raw <- sqrt(data_in$mu_df / (data_in$mu_df - 2))
        # m[g,] <- mu_raw + (data_in$sigma_bio / sd_raw) * mu_raw
    }
    m <- reshape2::melt(m)
    colnames(m) <- c("g", "p", "m")
    m$offset <- ifelse(m$g == data_in$offset, yes = 1, no = 0)

    
    w <- c()
    for(i in 1:nrow(m)) {
        if(m$offset[i]==1) {
            y <- rnorm(n = data_in$N_well_reps, 
                       mean = data_in$alpha_p[m$p[i]], 
                       sd = data_in$sigma_tech)
            y <- data.frame(mu = y, p = m$p[i], g = m$g[i], offset = m$offset[i])
            w <- rbind(w, y)
        } else {
            y <- rnorm(n = data_in$N_well_reps, 
                       mean = data_in$alpha_p[m$p[i]]+m$m[i], 
                       sd = data_in$sigma_tech)
            y <- data.frame(mu = y, p = m$p[i], g = m$g[i], offset = m$offset[i])
            w <- rbind(w, y)
        }
    }
    w$mu <- exp(w$mu)
    w$kappa <- exp(rnorm(n = nrow(w), 
                         mean = data_in$kappa_mu, 
                         sd = data_in$kappa_sigma))
    w$well_id <- 1:nrow(w)
    u <- do.call(rbind, lapply(X = 1:nrow(w), FUN = get_gamma, 
                               w = w, n_cells = data_in$N_well_cells))
    u$compound <- NA
    u$dose <- NA
    return(u)
}
