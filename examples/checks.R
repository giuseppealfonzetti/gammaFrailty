#### libraries #####
library(tidyverse)

#### choose true model ####

p <- 50
q <- 2

xi <- 2/q
rho <- .8

m <- 20
n <- 1000
int <- runif(p, -.5, .5)#rep(.5, p)#
b <-  rnorm(m, 0, .1) #runif(m,0,.5)#rep(0, m)
set.seed(1)
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)
d <- length(theta)

##### generate the data ###
seed <- 3
dt <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = seed
)

#### Derivatives #####
dt

i <- 1
j <- 2
jp <- 0
n_j <- dt[i + 1, j+1]
n_jp <- dt[i + 1, jp+1]
p
alpha_j <- int[j+1]
alpha_jp <- int[jp+1]
x <- X[i+1,]
beta <- b
xi <- 2/q
rho

Rwrapper_obj_ind <- function(par, ind){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2], ind = ind, verboseS = F, verboseSind = T
    )

    obj
}
Rwrapper_obj <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2], verboseS = T
    )

    obj
}
Rwrapper_obj(repar_theta)


Rwrapper_obj_fun <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2]
    )

    obj$ll
}
Rwrapper_obj_der <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2]
    )

    obj$gradient
}

#### Derivatives cl ####
Rwrapper_obj_fun(repar_theta)
numDeriv::grad(Rwrapper_obj_fun, repar_theta)
Rwrapper_obj_der(repar_theta)

ncl(theta, dt, X, T)
Rwrapper_ncl <- function(par){
    ncl(par, dt, X)$nll
}

Rwrapper_ngr <- function(par){
    ncl(par, dt, X)$ngradient
}

par_init <- rep(0, length(theta))
par_init <- repar_theta + runif(length(repar_theta), -1, 1);

Rwrapper_ncl(par = par_init)
numDeriv::grad(Rwrapper_ncl, par_init)
Rwrapper_ngr(par = par_init)

#### Opt ####
par_init <- repar_theta + runif(length(repar_theta), -1, 1);
system.time(
    opt <- ucminf::ucminf(par_init, Rwrapper_ncl, Rwrapper_ngr)
)
repartopar(opt$par)
theta

mean((repartopar(opt$par)-theta)^2)
mean((repartopar(par_init)-theta)^2)

Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
repartopar(Opt_u$theta)
theta


ctrl_sgd <- list(
    MAXT = 5000,
    BURN = 1000,
    STEPSIZE = .0001,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = seed
)
Opt_sgd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SCSD',
    CPP_CONTROL = ctrl_sgd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL#trajSub
)
repartopar(Opt_sgd$theta)
theta

Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)

repartopar(Opt_u$theta)
theta

Opt_u$theta
repar_theta

mean((repartopar(opt$par)-theta)^2)

maxT <- n; set.seed(1)
trajSub <- c(50, 100, 1000, 2000)
trajSub0 <- c(1, trajSub[trajSub < maxT]+1, maxT+1)
Opt <- gammaFrailty(
    THETA_INIT = par_init,
    DATA = dt,
    X = X,
    MAXT = maxT,
    BURN = maxT*(2/3),
    STEPSIZE = .001,
    STEPSIZE0 = .001,
    NU = 1,
    METHODFLAG = 1,
    VERBOSEFLAG = F,
    STEPSIZEFLAG = F
)
Opt$path_av_theta
