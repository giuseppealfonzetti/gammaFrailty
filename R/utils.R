
# theta reparameterisation
#'@export
partorepar <- function(par){
    repar <- par
    repar[1] <- -log( repar[1])
    #repar[1] <- log( repar[1])
    #repar[1] <- zofr( repar[1]-1.1)
    repar[2] <- zofr_cpp(repar[2])
    return(repar)
}

# inverse theta reparameterisation
#'@export
repartopar <- function(repar){
    par <- repar
    # par[1] <- exp(par[1])
    par[1] <- exp(-par[1])
    #par[1] <- rofz(par[1])+1.1
    par[2] <- rofz_cpp(par[2])
    return(par)
}

#'@export
check_SCSD_args <- function(ARGS, N){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N^1.5,0)
    if(is.null(ARGS$BURN)) out$BURN <- 0
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1e-3
    if(is.null(ARGS$NU)) out$NU <- 1
    if(is.null(ARGS$SEED)) out$SEED <- 123

    return(out)
}

#'@export
check_GD_args <- function(ARGS, N){

    out <- list('MAXT' = ARGS$MAXT, 'STEPSIZE' = ARGS$STEPSIZE)

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N^.75,0)
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1e-2

    return(out)
}
