
get_groups <- function(x) {
  check_generic(y = x)
  m <- x$posteriors$delta_t[, c("group_id", "group", "compound", "dose")]
  return(m)
}

get_pairs <- function(x, groups, exponentiate) {
  check_generic(y = x)
  check_logical(y = exponentiate, par = "exponentiate")
  
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
  
  idx <- expand.grid(i = seq_len(nrow(gmap)), j = seq_len(nrow(gmap)))
  
  ds <- do.call(rbind, mapply(FUN = function(i, j) {
    d <- p[, gmap$group_id[i]] - p[, gmap$group_id[j]]
    pmax <- get_pmax(d)
    d_M <- mean(d)
    d_HDI <- get_hdi(vec = d, hdi_level = 0.95)
    return(data.frame(group_id_x = gmap$group_id[i], 
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
                      pmax = pmax))}, 
    idx$i, idx$j, SIMPLIFY = FALSE))
  ds$group_x <- factor(x = ds$group_x, levels = groups)
  ds$group_y <- factor(x = ds$group_y, levels = groups)
  
  g_pi <- ggplot(data = ds)+
    geom_tile(aes(y = group_x, x = group_y, fill = pmax), col = "white")+
    geom_text(aes(y = group_x, x = group_y, 
                  label = round(x = pmax, digits = 2)), size = 2)+
    scale_fill_gradient(name = expression(pi), low = "white",high = "darkgray")+
    theme_bw(base_size = 10)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(label = '')+
    ylab(label = '')
  
  if(exponentiate==FALSE) {
    g <- ggplot(data = ds)+
      geom_tile(aes(y = group_x, x = group_y, fill = rho_M), col = "white")+
      geom_text(aes(y = group_x, x = group_y, 
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
      geom_tile(aes(y = group_x, x = group_y, fill = rho_M_exp), col = "white")+
      geom_text(aes(y = group_x, x = group_y, 
                    label = round(x = rho_M_exp, digits = 1)), size = 2)+
      scale_fill_distiller(name = expression(rho*"'"), palette = "Spectral")+
      scale_radius(name = expression(rho*"'"))+
      theme_bw(base_size = 10)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab(label = '')+
      ylab(label = '')
  }
  return(list(ds = ds, plot_rho = g, plot_pi = g_pi))
}

get_violins <- function(x, from_groups, to_group, exponentiate) {
  check_generic(y = x)
  check_logical(y = exponentiate, par = "exponentiate")
  
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
  
  idx <- expand.grid(i = seq_len(nrow(gmap_from)), j = seq_len(nrow(gmap_to)))
  ds <- do.call(rbind, mapply(function(i, j) {
    d <- p[, gmap_from$group_id[i]] - p[, gmap_to$group_id[j]]
    pmax <- get_pmax(d)
    return(data.frame(rho = d,
                      group_id_x = gmap_from$group_id[i], 
                      group_id_y = gmap_to$group_id[j], 
                      group_x = gmap_from$group[i], 
                      group_y = gmap_to$group[j], 
                      compound = gmap_from$compound[i], 
                      dose = gmap_from$dose[i],
                      pmax = pmax))
  }, idx$i, idx$j, SIMPLIFY = FALSE))
  
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
      geom_text(data = ds_pmax, aes(x = dose, y = max(ds$rho)+0.15, 
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
      geom_text(data = ds_pmax, aes(x = dose, y = max(exp(ds$rho))+0.15, 
                    label = round(x = pmax, digits = 2)), size = 2.25)+
      theme_bw(base_size = 10)+
      theme(strip.text.x = element_text(margin = margin(0.03,0,0.03,0,"cm")))+
      xlab(label = 'Dose')+
      ylab(label = expression("FC ("*rho*"')"))+
      ggtitle(label = paste0("treatments / ", to_group))
  }
  return(list(ds = ds, plot = g))
}
