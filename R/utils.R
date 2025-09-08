check_positive_integer <- function(y, par) {
    check_missing(y = y, par = par)
    check_length_one(y = y, par = par)
    check_na_val(y = y, par = par)
    check_numeric_val(y = y, par = par)
    check_inf_val(y = y, par = par)
    if(y <= 0) {
        str <- paste0(par, " must be >0")
        stop(str)
    }
    if(any(!(abs(y - round(y)) < .Machine$double.eps^0.5))) {
        str <- paste0(par, " must be an integer")
        stop(str)
    }
}
check_sd <- function(y, par) {
    check_missing(y = y, par = par)
    check_length_one(y = y, par = par)
    check_numeric_val(y = y, par = par)
    check_na_val(y = y, par = par)
    if(y < 0) {
        str <- paste0(par, " must be >0")
        stop(str)
    }
    check_inf_val(y = y, par = par)
}
check_loc <- function(y, par) {
    check_missing(y = y, par = par)
    check_length_one(y = y, par = par)
    check_na_val(y = y, par = par)
    check_numeric_val(y = y, par = par)
    check_inf_val(y = y, par = par)
}
check_delta <- function(x, o) {
    check_missing(y = x, par = "delta")
    if(length(x)<=1) {
        stop("delta must have length > 1")
    }
    if(any(is.na(x) | is.infinite(x))) {
        stop("some deltas are NA")
    }
    check_numeric_vec(y = x, par = "delta")
    
    # offset
    check_missing(y = o, par = "offset")
    check_length_one(y = o, par = "offset")
    check_na_val(y = o, par = "offset")
    check_numeric_val(y = o, par = "offset")
    if(o <= 0) {
        stop("offset must be greater than 0")
    }
    if(any(!(abs(o - round(o)) < .Machine$double.eps^0.5))) {
        stop("offset must be an integer")
    }
    if(o > length(x)) {
        stop("offset must be between 1 and length of delta")
    }
}
check_generic <- function(y) {
    check_missing(y = y, par = "x")
    check_list(y = y, par = "x")
    if(any(names(y)=="posteriors")==FALSE) {
        stop("wrong y")
    }
}

check_missing <- function(y, par) {
    if(missing(y) || is.null(y)) {
        str <- paste0(par, " is missing or NULL")
        stop(str)
    }
}
check_length_one <- function(y, par) {
    if(length(y)!=1) {
        str <- paste0(par, " must have length of 1")
        stop(str)
    }
}
check_list <- function(y, par) {
    if(is.list(y) == FALSE) {
        str <- paste0(par, " must be list")
        stop(str)
    }
}
check_dataframe <- function(y, par) {
    if(is.data.frame(y) == FALSE) {
        str <- paste0(par, " must be data.frame")
        stop(str)
    }
}
check_numeric_val <- function(y, par) {
    if(is.numeric(y)==FALSE) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
}
check_numeric_vec <- function(y, par) {
    if(any(is.numeric(y)==FALSE)) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
}
check_inf_val <- function(y, par) {
    if(is.infinite(y)) {
        str <- paste0(par, " is infinite")
        stop(str)
    }
}
check_na_val <- function(y, par) {
    if(is.na(y)) {
        str <- paste0(par, " is NA")
        stop(str)
    }
}
check_na_vec <- function(y, par) {
    if(any(is.na(y))) {
        str <- paste0(par, " contains NA")
        stop(str)
    }
}
check_logical_val <- function(y, par) {
    if(is.logical(y) == FALSE) {
        str <- paste0(par, " must be logical")
        stop(str)
    }
}
check_character_val <- function(y, par) {
    if(is.character(y) == FALSE) {
        str <- paste0(par, " must be character")
        stop(str)
    }
}
check_character_vec <- function(y, par) {
    if(any(is.character(y)) == FALSE) {
        str <- paste0(par, " must be character")
        stop(str)
    }
}
check_x_in_y <- function(x, y, e) {
    if(!x %in% y) {
        stop(e)
    }
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
    nCIs <- length(sortedPts) - ciIdxInc
    # vapply is preferred for type safety
    ciWidth <- vapply(seq_len(nCIs), function(i) {
        sortedPts[i + ciIdxInc] - sortedPts[i]
    }, numeric(1))
    HDImin <- sortedPts[which.min(ciWidth)]
    HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc]
    HDIlim <- c(HDImin, HDImax)
    return(HDIlim)
}


