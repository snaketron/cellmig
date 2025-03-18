
sim <- function(N_well_cells,
                N_plate, 
                N_group, 
                N_well_reps, 
                shape, 
                sigma_bplate, 
                sigma_wplate, 
                alpha_plate,
                mu_group) {
    
    get_meta <- function(data_in) {
        meta <- c()
        w_i <- 1
        for(p_i in 1:data_in$N_plate) {
            for(g_i in 1:data_in$N_group) {
                for(w in 1:data_in$N_well_reps) {
                    r <- data.frame(well = paste0("p",p_i,"|g",g_i,"|w",w_i),
                                    group = paste0("g", g_i),
                                    plate = paste0("p", p_i))
                    meta <- rbind(meta, r)
                    w_i = w_i + 1;
                }
            }
        }
        meta$group_id <- as.numeric(as.factor(meta$group))
        meta$well_id <- as.numeric(as.factor(meta$well))
        meta$plate_id <- as.numeric(as.factor(meta$plate))
        return(meta)
    }
    
    data_in <- list(N_plate = N_plate,
                    N_group = N_group,
                    N_well_reps = N_well_reps,
                    N_well_cells = N_well_cells,
                    shape = shape,
                    sigma_bplate = sigma_bplate,
                    sigma_wplate = sigma_wplate,
                    alpha_plate = alpha_plate,
                    mu_group = mu_group)
    
    message("simulation... \n")
    
    # sim data from model
    f <- sampling(object = stanmodels$S,
                  data = data_in,
                  chains = 1, 
                  cores = 1,
                  warmup = 1,
                  iter = N_well_cells+1,
                  algorithm = "Fixed_param")
    
    # prep data
    y <- extract(object = f, par = "y_hat_sample")$y_hat_sample
    y <- melt(y)
    colnames(y) <- c("iteration", "well_id", "v")
    y$iteration <- NULL
    
    # meta
    m <- get_meta(data_in = data_in)
    y <- merge(x = y, y = m, by.x = "well_id", by.y = "well_id", all.x = TRUE)
    
    return(y)
}
