
get_ppc_violins <- function(x, wrap = FALSE, ncol = 4) {
    e <- extract(object = x$f, par = "y_hat_sample")$y_hat_sample
    e <- melt(data = e)
    colnames(e) <- c("iter", "well_id", "yhat")
    q <- x$x$d[, c("well_id", "compound", "dose", "group", "plate", "well")]
    q <- q[duplicated(q)==F,]
    e <- merge(x = e, y = q, all.x = T)
    
    if(wrap==FALSE) {
        g <- ggplot()+
            facet_grid(compound~plate, scales = "free")+
            geom_sina(data = x$x$d, aes(x = as.factor(dose), 
                                        y = sv, group = well), 
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
            geom_sina(data = x$x$d, aes(x = as.factor(dose), 
                                        y = sv, group = well), 
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
    y <- aggregate(sv~well_id, data = x$x$d, FUN = mean)
    yhat <- x$s$yhat
    y <- merge(x = y, y = yhat, by = "well_id")
    
    g <- ggplot(data = y)+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
        geom_errorbar(aes(x = sv, y = mean, ymin = X2.5., ymax = X97.5.), 
                      col = "darkgray")+
        geom_point(aes(x = sv, y = mean), shape = 21, 
                   fill = "white", alpha = 0.75)+
        theme_bw(base_size = 10)+
        xlab(label = "Observed migration + 95% HDI [scaled]")+
        ylab(label = "Predicted migration [scaled]")+
        xlim(c(0,1))+
        ylim(c(0,1))
    
    return(g)
}
