check_positive_integer <- function(y, par) {
    if(missing(y) || is.null(y)) {
        str <- paste0(par, " is missing or NULL")
        stop(str)
    }
    if(length(y)!=1) {
        str <- paste0(par, " must have length of 1")
        stop(str)
    }
    if(is.na(y)) {
        str <- paste0(par, " is NA")
        stop(str)
    }
    if(is.numeric(y)==FALSE) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
    if(is.infinite(y)) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
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
    if(missing(y) || is.null(y)) {
        str <- paste0(par, " is missing or NULL")
        stop(str)
    }
    if(is.na(y)) {
        str <- paste0(par, " is NA")
        stop(str)
    }
    if(length(y)!=1) {
        str <- paste0(par, " must have length of 1")
        stop(str)
    }
    if(is.numeric(y)==FALSE) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
    if(y < 0) {
        str <- paste0(par, " must be >0")
        stop(str)
    }
}
check_loc <- function(y, par) {
    if(missing(y) || is.null(y)) {
        str <- paste0(par, " is missing or NULL")
        stop(str)
    }
    if(is.na(y)) {
        str <- paste0(par, " is NA")
        stop(str)
    }
    if(length(y)!=1) {
        str <- paste0(par, " must have length of 1")
        stop(str)
    }
    if(is.numeric(y)==FALSE) {
        str <- paste0(par, " must be numeric")
        stop(str)
    }
}
check_delta <- function(x, o) {
    # delta
    if(missing(x) || is.null(x)) {
        stop("delta is missing or NULL")
    }
    if(any(is.na(x) | is.infinite(x))) {
        stop("some deltas are NA")
    }
    if(length(x)<=1) {
        stop("delta must have length > 1")
    }
    if(any(is.numeric(x)==FALSE)) {
        stop("delta must be numeric")
    }
    
    
    # offset
    if(missing(o) || is.null(o)) {
        stop("offset is missing or NULL")
    }
    if(is.na(o)) {
        stop("offset is NA")
    }
    if(length(o)!=1) {
        stop("offset must have length of 1")
    }
    if(is.numeric(o)==FALSE) {
        stop("offset must be numberic")
    }
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
check_logical <- function(y, par) {
    if(missing(y) || is.null(y)) {
        str <- paste0(par, " is missing or NULL")
        stop(str)
    }
    if(length(y) != 1) {
        str <- paste0(par, " must be have length one")
        stop(str)
    }
    if(is.logical(y)==FALSE) {
        str <- paste0(par, " must be logcal")
        stop(str)
    }
}
check_generic <- function(y) {
    if(missing(y) || is.null(y)) {
        stop("missing x")
    }
    if(is.list(y) == FALSE) {
        stop("wrong x")
    }
    if(any(names(y)=="posteriors")==FALSE) {
        stop("wrong y")
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
    ciWidth <- rep(0 , nCIs)
    for (i in seq_len(nCIs)) {
        ciWidth[i] <- sortedPts[i + ciIdxInc] - sortedPts[i]
    }
    HDImin <- sortedPts[which.min(ciWidth)]
    HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc]
    HDIlim <- c(HDImin, HDImax)
    return(HDIlim)
}
