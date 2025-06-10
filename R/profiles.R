
get_dose_response_profile <- function(x,
                                      hc_link = "average",
                                      hc_dist = "euclidean",
                                      select_gs,
                                      B = 1000) {
  
  
  get_boot <- function(x, hc_dist, hc_link, B) {
    
    eg <- x$posteriors$mu_group
    gs <- eg$group_id
    
    # hclust
    q <- acast(data = eg, formula = compound~dose, value.var = "mean")
    hc <- hclust(dist(q, method = hc_dist), method = hc_link)
    main_ph <- as.phylo(x = hc)
    
    # meta
    meta <- x$posteriors$mu_group[, c("group_id", "compound", "dose")]
    meta <- meta[order(meta$group_id, decreasing = F),]
    meta <- meta[meta$group_id %in% gs, ]
    meta$group_id <- NULL
    
    # extract posterior
    e <- extract(x$f, par = "mu_group")$mu_group[, gs]
    e <- e[sample(x = 1:nrow(e), size = min(nrow(e), B), replace = TRUE),]
    
    boot_ph <- c()
    for(i in 1:nrow(e)) {
      
      u <- data.frame(g = 1:ncol(e), mean = e[i, ])
      u <- cbind(u, meta)
      
      # hclust
      q <- acast(data = eg, formula = compound~dose, value.var = "mean")
      hc <- hclust(dist(q, method = hc_dist), method = hc_link)
      ph <- as.phylo(x = hc)
      
      if(i == 1) {
        boot_ph <- ph
      }
      else {
        boot_ph <- c(boot_ph, ph)
      }
    }
    
    clades <- prop.clades(phy = main_ph, x = boot_ph, part = NULL,
                          rooted = is.rooted(main_ph))
    
    # add bootstrap
    main_ph$node.label <- clades
    
    # b = 0 for these nodes
    na_nodes <- which(is.na(main_ph$node.label))
    if(length(na_nodes)!=0) {
      main_ph$node.label[na_nodes] <- 0
    }
    
    return(list(main_ph = main_ph, boot_ph = boot_ph))
  }
  
  eg <- x$posteriors$mu_group
  es <- x$posteriors$mu_plate_group
  if(missing(select_gs)==FALSE) {
    if(any(!select_gs %in% unique(eg$group))) {
      stop("selected treatment groups not found in data")
    }
    eg <- eg[eg$group %in% select_gs, ]
    es <- es[es$group %in% select_gs, ]
  }
  
  bt <- get_boot(x = x, hc_dist = hc_dist, hc_link = hc_link, B = B)
  
  tree <- ggtree(bt$main, linetype='solid')+
    geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
    geom_tippoint(size = 2, fill = "white", shape = 21)+
    geom_tiplab(color='black', as_ylab = T, align = TRUE)+
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
  
  g <- ggplot(data = q)+
    facet_grid(compound~., switch = "y")+
    geom_hline(yintercept = 0, linetype = "dashed", col = "gray")+
    geom_point(aes(x = dose, y = mean))+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), 
                  width = 0)+
    scale_y_continuous(position = "right", 
                       breaks = pretty_breaks(n = 5))+
    theme_bw(base_size = 10)+
    theme(strip.text.y = element_text(
      margin = margin(0.01,0.01,0.01,0.01, "cm")))
  
  
  q <- es[es$compound %in% q$compound, ]
  q$compound <- factor(q$compound, levels = rev(tips))
  g2 <- ggplot(data = q)+
    geom_hline(yintercept = 0, linetype = "dashed", col = "gray")+
    facet_wrap(facets = compound~plate, 
               nrow = length(unique(q$compound)), 
               switch = "y", scales = "free_y")+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), 
                  width = 0, alpha = 0.5)+
    geom_line(aes(x = dose, y = mean))+
    geom_point(aes(x = dose, y = mean))+
    scale_y_continuous(position = "right", 
                       breaks = pretty_breaks(n = 3))+
    theme_bw(base_size = 10)+
    theme(legend.position = "none", strip.text.y = element_text(
      margin = margin(0.01,0.01,0.01,0.01, "cm")))
  
  
  gs <- (tree|g|g2)+
    plot_annotation(tag_levels = 'A')
  
  gs
  return(gs)
}

get_treatment_profile <- function(x, 
                                  hc_link = "average", 
                                  hc_dist = "euclidean",
                                  select_gs,
                                  B = 1000) {
  
  get_boot <- function(x, hc_dist, hc_link, select_gs, B) {
    
    eg <- x$posteriors$mu_group
    
    if(missing(select_gs)==FALSE) {
      if(any(!select_gs %in% unique(eg$group))) {
        stop("selected treatment groups not found in data")
      }
      eg <- eg[eg$group %in% select_gs, ]
    }
    
    # hclust -> main tree
    gs <- eg$group_id
    gns <- eg$group
    v <- eg$mean
    names(v) <- eg$group_id
    hc <- hclust(dist(v, method = hc_dist), method = hc_link)
    main_ph <- as.phylo(x = hc)
    
    
    # extract posterior
    e <- extract(x$f, par = "mu_group")$mu_group[, gs]
    e <- e[sample(x = 1:nrow(e), size = min(nrow(e), B), replace = TRUE),]
    
    boot_ph <- c()
    for(i in 1:nrow(e)) {
      
      # hclust
      hc <- hclust(dist(e[i,], method = hc_dist), method = hc_link)
      ph <- as.phylo(x = hc)
      
      if(i == 1) {
        boot_ph <- ph
      }
      else {
        boot_ph <- c(boot_ph, ph)
      }
    }
    
    clades <- prop.clades(phy = main_ph, x = boot_ph, part = NULL,
                          rooted = is.rooted(main_ph))
    
    # add bootstrap
    main_ph$node.label <- clades
    
    if(all(main_ph$tip.label == gs)) {
      main_ph$tip.label <- gns
    }
    
    # b = 0 for these nodes
    na_nodes <- which(is.na(main_ph$node.label))
    if(length(na_nodes)!=0) {
      main_ph$node.label[na_nodes] <- 0
    }
    
    return(list(main_ph = main_ph, boot_ph = boot_ph))
  }
  
  bt <- get_boot(x = x, hc_dist = hc_dist, hc_link = hc_link,
                 select_gs = select_gs, B = B)
  
  tree <- ggtree(bt$main, linetype='solid')+
    geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
    geom_tippoint(size = 2, fill = "white", shape = 21)+
    geom_tiplab(color='black', as_ylab = T, align = TRUE)+
    layout_rectangular()+
    theme_bw(base_size = 10)+
    scale_x_continuous(labels = abs)+
    geom_nodelab(geom='text', color = "#4c4c4c" ,size = 2.75, hjust=-0.2,
                 mapping = aes(label=label,subset=isTip==FALSE))
  
  tree <- revts(tree)
  return(tree)
}
