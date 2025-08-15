
get_ppc_violins <- function(x, wrap = FALSE, ncol = 4) {
    
    if(missing(x) || is.null(x)) {
        stop("missing x")
    }
    if(missing(wrap) || is.null(wrap)) {
        stop("missing wrap")
    }
    
    if(length(wrap)!=1) {
        stop("wrap must be TRUE or FALSE")
    }
    if(is.logical(wrap)==FALSE | is.na(wrap)) {
        stop("wrap must be logical")
    }

    if(wrap) {
        if(missing(ncol) || is.null(ncol)) {
            stop("missing ncol")
        }
        if(length(ncol)!=1) {
            stop("ncol must be positive integer")
        }
        if(is.numeric(ncol)==FALSE) {
            stop("ncol must be positive integer")
        }
        if(ncol<=0 | is.na(ncol) | is.infinite(ncol)) {
            stop("ncol must be positive integer")
        }
    }
    
    
    
    e <- extract(object = x$fit, par = "y_hat_sample")$y_hat_sample
    e <- melt(data = e)
    colnames(e) <- c("iter", "well_id", "yhat")
    q <- x$data$proc_x[, c("well_id", "compound", "dose", 
                           "group", "plate", "well")]
    q <- q[duplicated(q)==FALSE,]
    e <- merge(x = e, y = q, all.x = TRUE)
    
    if(wrap==FALSE) {
        g <- ggplot()+
            facet_grid(compound~plate, scales = "free")+
            geom_sina(data = x$data$proc_x, aes(x = as.factor(dose), 
                                        y = v, group = well), 
                      col = "black", size = 0.3)+
            geom_violin(data = e, aes(x = as.factor(dose), 
                                      y = yhat, group = well), 
                        fill = NA, col = "#f75ea3", alpha = 0.35)+
            theme_bw(base_size = 10)+
            theme(legend.position = "none")+
            theme(strip.text.x = element_text(
                margin = margin(0.02,0,0.02,0, "cm")))+
            xlab(label = "dose")+
            scale_y_continuous(name = "Cell migration speed", 
                               breaks = pretty_breaks(3))
    }
    
    if(wrap) {
        g <- ggplot()+
            facet_wrap(~compound+plate, scales = "free", ncol = ncol)+
            geom_sina(data = x$data$proc_x, aes(x = as.factor(dose), 
                                        y = v, group = well), 
                      col = "black", size = 0.3)+
            geom_violin(data = e, aes(x = as.factor(dose), 
                                      y = yhat, group = well), 
                        fill = NA, col = "#f75ea3", alpha = 0.35)+
            theme_bw(base_size = 10)+
            theme(legend.position = "none")+
            theme(strip.text.x = element_text(
                margin = margin(0.02,0,0.02,0, "cm")))+
            xlab(label = "dose")+
            scale_y_continuous(name = "Cell migration speed", 
                               breaks = pretty_breaks(3))
    }
    
    return(g)
}


get_ppc_means <- function(x) {
    
    if(missing(x) || is.null(x)) {
        stop("missing x")
    }
    
    y <- aggregate(v~well_id, data = x$data$proc_x, FUN = mean)
    yhat <- x$posteriors$yhat
    y <- merge(x = y, y = yhat, by = "well_id")
    
    g <- ggplot(data = y)+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
        geom_errorbar(aes(x = v, y = mean, ymin = X2.5., ymax = X97.5.), 
                      col = "darkgray")+
        geom_point(aes(x = v, y = mean), shape = 21, 
                   fill = "white", alpha = 0.75)+
        theme_bw(base_size = 10)+
        xlab(label = "Observed migration + 95% HDI")+
        ylab(label = "Predicted migration")
    
    return(g)
}
