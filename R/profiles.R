
get_dose_response_profile <- function(x,
                                      hc_link = "average",
                                      hc_dist = "euclidean",
                                      groups,
                                      B = 1000, 
                                      exponentiate) {
    
    check_profile_inputs(x = x, hc_link = hc_link, hc_dist = hc_dist, B = B)
    check_logical_val(y = exponentiate, par = "exponentiate")
    
    gmap <- get_groups(x = x)
    if(missing(groups)) {
        warning("groups not specified, we will use all groups")
        groups <- gmap$group
    }
    if(any(is.na(groups))|any(is.nan(groups))) {
        warning("groups not specified, we will use all groups")
        groups <- gmap$group
    }
    if(any(is.character(groups)==FALSE)) {
        warning("groups must be characters")
    }
    if(length(groups)==1) {
        stop("only one treatment groups provided, length(groups)>1")
    }
    eg <- x$posteriors$delta_t
    es <- x$posteriors$delta_tp
    if(any(!groups %in% unique(eg$group))) {
        stop("selected groups not found in data")
    }
    eg <- eg[eg$group %in% groups, ]
    es <- es[es$group %in% groups, ]
    
    bt <- get_boot_drp(x = x, hc_dist = hc_dist, hc_link = hc_link, 
                       groups = groups, B = B)
    
    tree <- ggtree(bt$main, linetype='solid')+
        geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col="black")+
        geom_tippoint(size = 2, fill = "white", shape = 21)+
        geom_tiplab(color='black', as_ylab = TRUE, align = TRUE)+
        layout_rectangular()+
        theme_bw(base_size = 10)+
        scale_x_continuous(labels = abs)+
        geom_nodelab(geom='text', color = "#4c4c4c" ,size = 2.75, hjust=-0.2,
                     mapping = aes(label=label,subset=isTip==FALSE))
    tree <- revts(tree)
    
    t <- tree$data
    t <- t[order(t$y, decreasing = FALSE), ]
    tips <- t$label[t$isTip==TRUE]
    
    q <- eg
    q$compound <- factor(q$compound, levels = rev(tips))
    
    ylab <- expression(delta[t])
    yhline <- 0
    if(exponentiate) {
        q$mean <- exp(q$mean)
        q$X2.5. <- exp(q$X2.5.)
        q$X97.5. <- exp(q$X97.5.)
        ylab <- expression(delta[t]^"'")
        yhline <- 1
    }
    
    g <- ggplot(data = q)+
        facet_grid(compound~., switch = "y")+
        geom_hline(yintercept = yhline, linetype = "dashed", col = "gray")+
        geom_point(aes(x = dose, y = mean), size = 1)+
        geom_errorbar(aes(x=dose, y=mean, ymin=X2.5., ymax = X97.5.),width=0)+
        scale_y_continuous(name = ylab,
                           position = "right", breaks = pretty_breaks(n = 4))+
        theme_bw(base_size = 10)+
        theme(strip.text.y=element_text(
            margin=margin(0.01,0.01,0.01,0.01,"cm")))+
        xlab(label = "Dose")
    
    q <- es[es$compound %in% q$compound, ]
    q$compound <- factor(q$compound, levels = rev(tips))
    ylab <- expression(delta[t])
    yhline <- 0
    if(exponentiate) {
        q$mean <- exp(q$mean)
        q$X2.5. <- exp(q$X2.5.)
        q$X97.5. <- exp(q$X97.5.)
        ylab <- expression(delta[t]^"'")
        yhline <- 1
    }
    g2 <- ggplot(data = q)+
        geom_hline(yintercept = yhline, linetype = "dashed", col = "gray")+
        facet_grid(compound~., switch = "y")+
        geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.,
                          group = plate), position = position_dodge(width=0.8), 
                      width = 0, alpha = 0.5)+
        geom_point(aes(x = dose, y = mean, group = plate),
                   position = position_dodge(width = 0.8), size = 1)+
        scale_y_continuous(name = ylab,
                           position = "right", breaks = pretty_breaks(n = 4))+
        theme_bw(base_size = 10)+
        theme(legend.position = "none", strip.text.y = element_text(
            margin = margin(0.01,0.01,0.01,0.01, "cm")))+
        xlab(label = "Dose")
    return((tree|g|g2)+
               plot_layout(widths = c(1,1,2))+
               plot_annotation(tag_levels = 'A'))
}


check_profile_inputs <- function(x, 
                                 hc_link, 
                                 hc_dist, 
                                 B) {
    check_generic(y = x)
    
    check_missing(y = hc_link, par = "hc_link")
    check_length_one(y = hc_link, par = "hc_link")
    if(is.character(hc_link)==FALSE) {
        stop("hclust_method must be one of: ward.D, single, complete, 
    average, mcquitty, median, centroid or ward.D2")
    }
    if(!hc_link %in% c("ward.D", "ward.D2", "single", 
                       "complete", "average", "mcquitty",
                       "median", "centroid")) {
        stop("hc_link must be one of: ward.D, single, complete, 
    average, mcquitty, median, centroid or ward.D2")
    }
    
    check_missing(y = hc_dist, par = "hc_dist")
    check_length_one(y = hc_dist, par = "hc_dist")
    if(is.character(hc_dist)==FALSE) {
        stop("hc_dist must be one of: euclidean or manhattan")
    }
    if(!hc_dist %in% c("euclidean", "manhattan")) {
        stop("hc_dist must be one of: euclidean or manhattan")
    }
    
    check_positive_integer(y = B, par = "B")
}

get_boot_drp <- function(x, 
                         hc_dist, 
                         hc_link, 
                         groups,
                         B) {
    
    eg <- x$posteriors$delta_t
    
    if(missing(groups)) { 
        groups <- eg$group
    }
    
    if(any(!groups %in% unique(eg$group))) {
        stop("selected treatment groups not found in data")
    }
    eg <- eg[eg$group %in% groups, ]
    
    # meta
    meta <- eg[, c("group_id", "compound", "dose")]
    meta <- meta[order(meta$group_id, decreasing = FALSE),]
    
    # extract posterior
    e <- extract(x$f, par = "delta_t")$delta_t[, eg$group_id]
    e <- e[sample(x = seq_len(nrow(e)), size = min(nrow(e), B),replace = TRUE),]
    
    # summaries
    q <- acast(data = eg, formula = compound~dose, value.var = "mean")

    hc <- hclust(dist(q, method = hc_dist), method = hc_link)
    main_ph <- as.phylo(x = hc)
    
    boot_phs <- lapply(seq_len(nrow(e)), function(i) {
        u <- data.frame(g = seq_len(ncol(e)), mean = e[i, ])
        u <- cbind(u, meta)
        q <- acast(data = u, formula = compound~dose, value.var = "mean")
        hc <- hclust(dist(q, method = hc_dist), method = hc_link)
        return(as.phylo(x = hc))
    })
    
    clades <- prop.clades(phy = main_ph, x = c(boot_phs), part = NULL,
                          rooted = is.rooted(main_ph))
    
    # add bootstrap
    main_ph$node.label <- clades
    
    # b = 0 for these nodes
    na_nodes <- which(is.na(main_ph$node.label))
    if(length(na_nodes)!=0) {
        main_ph$node.label[na_nodes] <- 0
    }
    
    return(list(main_ph = main_ph, boot_ph = c(boot_phs)))
}

