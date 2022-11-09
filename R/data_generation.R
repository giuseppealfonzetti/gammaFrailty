#### functions #####
#'@export
generate_C <- function(rho, p){
    out <- diag(rep(1,p))
    for (j in 2:p) {
        for (k in 1:(j-1)) {
            out[j,k] <- rho^(abs(j-k)/2)
        }
    }
    out <- out+t(out)-diag(rep(1,p))

    return(out)
}

#'@export
generate_mgamma <- function(q, C, seed){
    set.seed(seed)
    out <- colMeans((mvtnorm::rmvnorm(n = q, sigma = C))^2)
    return(out)
}

#'@export
generate_data <- function(INTERCEPT, BETA, X, Q, RHO, SEED){

    if(!is.matrix(X)){cat('Error: X must be a matrix!\n'); return(NULL)}
    if(length(BETA)!=ncol(X)){cat('Error: Dimensions of beta and X not compatible!\n'); return(NULL)}
    if(RHO < 0){cat('Error: rho must be non-negative!\n'); return(NULL)}

    p <- length(INTERCEPT)
    C <- generate_C(rho = RHO, p = p)
    n <- nrow(X)
    Z <- t(sapply(1:n, function(x) generate_mgamma(Q, C, seed = SEED + x)))
    u <- t(sapply(1:n, function(i) exp(INTERCEPT+ as.numeric(crossprod(BETA, X[i,])))))



    set.seed(SEED)
    #Z2 <- matrix(rep(Z[1,], n), nrow = n, ncol =  p, byrow = T)
    out <- purrr::modify2(Z, u, ~rpois(1, lambda = .x * .y))

    return(out)
}

#'@export
uni_dnbinom <- function(u, n, eps){
    nf <- factorial(n)

    pr <- 1
    for (i in 0:(n-1)) {
        pr <- pr * (1+i*eps)
    }

    frac <- (u^n)/(nf *(1+eps*u)^(n+1/eps))
    out <- frac * pr

    return(as.numeric(out))
}
