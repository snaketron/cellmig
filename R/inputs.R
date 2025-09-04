
process_input <- function(x) {

  check_x <- function(x) {
    if(missing(x)) {
      stop("missing x")
    }
    if(is.data.frame(x)==FALSE) {
      stop("x must be a data.frame")
    }
    if(nrow(x)<=1) {
      stop("zero or one row in x")
    }
  }

  check_x(x=x)
  check_input_columns(x=x)
  
  if(any(x$offset == 1)) {
    return(get_input_with_offset(x = x))
  } else {
    return(get_input_without_offset(x = x))
  }
}

process_control <- function(control_in) {
  control <- get_default_control()
  
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
  control[names(control_in)] <- control_in

  check_positive_integer(y = control$mcmc_chains, par = "mcmc_chains")
  check_positive_integer(y = control$mcmc_cores, par = "mcmc_cores")
  check_positive_integer(y = control$mcmc_steps, par = "mcmc_steps")
  check_positive_integer(y = control$mcmc_warmup, par = "mcmc_warmup")
  
  check_mcmc_steps(mcmc_steps = control$mcmc_steps, 
                   mcmc_warmup = control$mcmc_warmup)

  return(control)
}

# MCMC steps vs. warmup
check_mcmc_steps <- function(mcmc_steps, mcmc_warmup) {
  if(as.integer(x = mcmc_steps) < 500 |
     as.integer(x = mcmc_warmup) < 100) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }

  if(as.integer(x = mcmc_steps) <= as.integer(x = mcmc_warmup)) {
    stop("mcmc_steps > mcmc_warmup")
  }
}

get_default_control <- function() {
  return(list(mcmc_warmup = 500,
              mcmc_steps = 1500,
              mcmc_chains = 4,
              mcmc_cores = 1,
              mcmc_algorithm = "NUTS",
              adapt_delta = 0.8,
              max_treedepth = 10,
              prior_alpha_p_M = 1.7,
              prior_alpha_p_SD = 1.0,
              prior_sigma_bio_M = 0.0,
              prior_sigma_bio_SD = 1.0,
              prior_sigma_tech_M = 0.0,
              prior_sigma_tech_SD = 1.0,
              prior_kappa_mu_M = 1.7,
              prior_kappa_mu_SD = 1.0,
              prior_kappa_sigma_M = 0.0,
              prior_kappa_sigma_SD = 1.0,
              prior_mu_group_M = 0.0,
              prior_mu_group_SD = 1.0))
}

check_input_columns <- function(x) {
    if(!"compound" %in% colnames(x)) {
      stop("x does not have a column compound")
    }
    if(!"dose" %in% colnames(x)) {
      stop("x does not have a column dose")
    }
    if(!"well" %in% colnames(x)) {
      stop("x does not have a column well")
    }
    if(!"plate" %in% colnames(x)) {
      stop("x does not have a column plate")
    }
    if(!"v" %in% colnames(x)) {
      stop("x does not have a column v")
    }
    if(!"offset" %in% colnames(x)) {
      stop("x does not have a column offset")
    }
    
    if(is.character(x[,"compound"])==FALSE) {
      stop("column compound must be character")
    }
    if(is.character(x[,"dose"])==FALSE) {
      stop("column dose must be character")
    }
    if(is.character(x[,"well"])==FALSE) {
      stop("column well must be character")
    }
    if(is.character(x[,"plate"])==FALSE) {
      stop("column plate must be character")
    }
    if(is.numeric(x[,"v"])==FALSE) {
      stop("column v must be numeric")
    }
    if(is.numeric(x[,"offset"])==FALSE) {
      stop("column offset must be numeric")
    }
    
    if(any(is.na(x[,"compound"]))) {
      stop("column compound contains NAs")
    }
    if(length(unique(x$compound))==1) {
      stop("only one compound included")
    }
    
    if(any(is.na(x[,"dose"]))) {
      stop("column dose contains NAs")
    }
    if(any(is.na(x[,"plate"]))) {
      stop("column plate contains NAs")
    }
    if(any(is.na(x[,"well"]))) {
      stop("column well contains NAs")
    }
    if(any(is.na(x[,"v"]))) {
      stop("column v contains NAs")
    }
    if(any(is.na(x[,"offset"]))) {
      stop("column offset contains NAs")
    }
    if(all(x[,"offset"] %in% c(0,1))==FALSE) {
      stop("column offset can contain only 0 or 1")
    }
  }

get_input_with_offset <- function(x) {
  org_x <- x
  
  x$plate_id <- as.numeric(as.factor(x$plate))
  x$group <- paste0(x$compound, '|', x$dose)
  x$well <- paste0(x$plate, '|', x$well, '|', x$group)
  x$plate_group <- paste0(x$plate, '|', x$group)
  x$well_id <- as.numeric(as.factor(x$well))
  
  xs <- x[x$offset == 0,]
  xr <- x[x$offset == 1,]
  
  xs$group_id <- as.numeric(as.factor(xs$group))
  xr$group_id <- 0
  
  xs$plate_group_id <- as.numeric(as.factor(xs$plate_group))
  xr$plate_group_id <- 0
  
  x <- rbind(xr, xs)
  
  x$N <- nrow(x)
  x$N_well <- max(x$well_id)
  x$N_plate <- max(x$plate_id)
  x$N_plate_group <- max(x$plate_group_id, na.rm = TRUE)
  x$N_group <- max(x$group_id)
  
  q <- x[, c("well_id", "plate_id", "plate_group_id", "offset")]
  q <- q[duplicated(q) == FALSE, ]
  q <- q[order(q$well_id, decreasing = FALSE),]
  
  # group_id -> plate_group_id 
  z <- x[, c("plate_group_id", "group_id")]
  z <- z[duplicated(z) == FALSE, ]
  z <- z[order(z$plate_group_id, decreasing = FALSE),]
  z <- z[z$group_id !=0,]
  
  return(list(proc_x = x, org_x = org_x, map_w = q, map_pg = z))
}

get_input_without_offset <- function(x) {
  org_x <- x
  
  x$plate_id <- as.numeric(as.factor(x$plate))
  
  x$group <- paste0(x$compound, '|', x$dose)
  x$group_id <- as.numeric(as.factor(x$group))
  
  x$well <- paste0(x$plate, '|', x$well, '|', x$group)
  x$well_id <- as.numeric(as.factor(x$well))
  
  x$plate_group <- paste0(x$plate, '|', x$group)
  x$plate_group_id <- as.numeric(as.factor(x$plate_group))
  
  x$N <- nrow(x)
  x$N_well <- max(x$well_id)
  x$N_plate <- max(x$plate_id)
  x$N_plate_group <- max(x$plate_group_id, na.rm = TRUE)
  x$N_group <- max(x$group_id)
  
  q <- x[, c("well_id", "plate_id", "plate_group_id", "offset")]
  q <- q[duplicated(q) == FALSE, ]
  q <- q[order(q$well_id, decreasing = FALSE),]
  
  # group_id -> plate_group_id 
  z <- x[, c("plate_group_id", "group_id")]
  z <- z[duplicated(z) == FALSE, ]
  z <- z[order(z$plate_group_id, decreasing = FALSE),]
  
  return(list(proc_x = x, org_x = org_x, map_w = q, map_pg = z))
}