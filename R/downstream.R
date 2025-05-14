


get_groups <- function(x) {
  m <- x$posteriors$mu_group[, c("group_id", "group",
                        "compound", "dose")]
  
  return(m)
}



get_pairs <- function(x, groups = NA, exponentiate = FALSE) {
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
  
  p <- extract(x$fit, par = "mu_group")$mu_group
  if(ncol(p)==1) {
    stop("only one treatment group: nothing to compare")
  }
  if(all(groups %in% gmap$group)==FALSE) {
    stop("unknown group in groups")
  }
  
  gmap <- gmap[gmap$group %in% groups,]
  gmap <- gmap[match(groups, gmap$group),]
  
  ds <- vector(mode = "list", length = nrow(gmap)*nrow(gmap))
  ct <- 1
  for(i in 1:nrow(gmap)) {
    for(j in 1:nrow(gmap)) {
      d <- p[,gmap$group_id[i]]-p[,gmap$group_id[j]]
      pmax <- get_pmax(d)
      d_M <- mean(d)
      d_HDI <- get_hdi(vec = d, hdi_level = 0.95)
      ds[[ct]] <- data.frame(group_id_x = gmap$group_id[i], 
                             group_id_y = gmap$group_id[j], 
                             group_x = gmap$group[i], 
                             group_y = gmap$group[j], 
                             compound = gmap$compound[i], 
                             dose = gmap$dose[i],
                             rho_M = d_M,
                             rho_L95 = d_HDI[1],
                             rho_H95 = d_HDI[2],
                             rho_M_exp = exp(d_M),
                             rho_L95_exp = exp(d_HDI[1]),
                             rho_H95_exp = exp(d_HDI[2]),
                             pmax = pmax)
      ct <- ct + 1
    }
  }
  ds <- do.call(rbind, ds)
  
  ds$group_x <- factor(x = ds$group_x, levels = groups)
  ds$group_y <- factor(x = ds$group_y, levels = groups)
  
  g_pi <- ggplot(data = ds)+
    geom_tile(aes(x = group_x, y = group_y, fill = pmax), col = "white")+
    geom_text(aes(x = group_x, y = group_y, 
                  label = round(x = pmax, digits = 2)), size = 2)+
    scale_fill_gradient(name = expression(pi), low = "white",high = "darkgray")+
    theme_bw(base_size = 10)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(label = '')+
    ylab(label = '')
  
  if(exponentiate==FALSE) {
    g <- ggplot(data = ds)+
      geom_tile(aes(x = group_x, y = group_y, fill = rho_M), col = "white")+
      geom_text(aes(x = group_x, y = group_y, 
                    label = round(x = rho_M, digits = 1)), size = 2)+
      scale_fill_distiller(name = expression(rho), palette = "Spectral")+
      scale_radius(name = expression(rho))+
      theme_bw(base_size = 10)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab(label = '')+
      ylab(label = '')
  } 
  else {
    g <- ggplot(data = ds)+
      geom_tile(aes(x = group_x, y = group_y, fill = rho_M_exp), col = "white")+
      geom_text(aes(x = group_x, y = group_y, 
                    label = round(x = rho_M_exp, digits = 1)), size = 2)+
      scale_fill_distiller(name = expression(rho*"'"), palette = "Spectral")+
      scale_radius(name = expression(rho*"'"))+
      theme_bw(base_size = 10)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab(label = '')+
      ylab(label = '')
  }
  
  plot <- g|g_pi
  
  return(list(ds = ds, plot = plot))
}



get_violins <- function(x, from_groups, to_group, exponentiate = FALSE) {
    if(length(to_group)!=1) {
        stop("only one to_group allowed")
    }
    
    p <- extract(x$fit, par = "mu_group")$mu_group
    if(ncol(p)==1) {
        stop("only one treatment group: nothing to compare")
    }
    
    gmap <- get_groups(x = x)
    
    if(all(from_groups %in% gmap$group)==FALSE) {
        stop("unknown group in from_groups")
    }
    if(all(to_group %in% gmap$group)==FALSE) {
        stop("unknown group in to_group")
    }
    
    gmap_from <- gmap[gmap$group %in% from_groups,]
    gmap_to <- gmap[gmap$group %in% to_group,]
    
    ds <- vector(mode = "list", length = nrow(gmap_from))
    ct <- 1
    for(i in 1:nrow(gmap_from)) {
        for(j in 1:nrow(gmap_to)) {
            
            d <- p[,gmap_from$group_id[i]]-p[,gmap_to$group_id[j]]
            pmax <- get_pmax(d)
            ds[[ct]] <- data.frame(rho = d,
                                   group_id_x = gmap_from$group_id[i], 
                                   group_id_y = gmap_to$group_id[j], 
                                   group_x = gmap_from$group[i], 
                                   group_y = gmap_to$group[j], 
                                   compound = gmap_from$compound[i], 
                                   dose = gmap_from$dose[i],
                                   pmax = pmax)
            ct <- ct + 1
        }
    }
    ds <- do.call(rbind, ds)
    ds$contrast <- paste0(ds$group_x, "-vs-", ds$group_y) 
    ds$contrast <- factor(x = ds$contrast, 
                          levels = paste0(from_groups, "-vs-", to_group))
    
    ds_pmax <- ds[duplicated(ds[, c("group_x", "group_y")])==FALSE,]
    
    if(exponentiate==FALSE) {
      g <- ggplot(data = ds)+
        facet_wrap(facets = ~compound, scales = "free_x", nrow = 1)+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_violin(aes(x = dose, y = rho), 
                    col = "steelblue", fill = "steelblue", alpha = 0.7)+
        geom_text(data = ds_pmax,
                  aes(x = dose, y = max(ds$rho)+0.15, 
                      label = round(x = pmax, digits = 2)), size = 2.25)+
        theme_bw(base_size = 10)+
        theme(strip.text.x = element_text(margin = margin(0.03,0,0.03,0,"cm")))+
        xlab(label = 'Dose')+
        ylab(label = expression("LFC ("*rho*")"))+
        ggtitle(label = paste0("treatments / ", to_group))
    } else {
      g <- ggplot(data = ds)+
        facet_wrap(facets = ~compound, scales = "free_x", nrow = 1)+
        geom_hline(yintercept = 1, linetype = "dashed")+
        geom_violin(aes(x = dose, y = exp(rho)), 
                    col = "steelblue", fill = "steelblue", alpha = 0.7)+
        geom_text(data = ds_pmax,
                  aes(x = dose, y = max(exp(ds$rho))+0.15, 
                      label = round(x = pmax, digits = 2)), size = 2.25)+
        theme_bw(base_size = 10)+
        theme(strip.text.x = element_text(margin = margin(0.03,0,0.03,0,"cm")))+
        xlab(label = 'Dose')+
        ylab(label = expression("FC ("*rho*"')"))+
        ggtitle(label = paste0("treatments / ", to_group))
    }
    
    g
    
    return(list(ds = ds, plot = g))
}



get_pmax <- function(x) {
  if(all(x==0)) {
    return(0)
  }
  l <- length(x)
  return(2*max(sum(x<0)/l, sum(x>0)/l)-1)
}

# Description:
# Computes HDI for vector vec and hdi_level (e.g. 0.95)
# Taken (and renamed) from "Doing Bayesian Analysis", section 25.2.3 R code
# for computing HDI of a MCMC sample
get_hdi <- function(vec, hdi_level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi_level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}

