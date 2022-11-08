#### libraries #####
library(tidyverse)
######

p <- 12
q <- 8

eps <- 2/q
rho <- .6

m <- 2
n <- 250
int <- rep(log(4), p)#runif(p, 0, 1)
b <- rep(0, m) #rnorm(m, 0, .5) #
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
#X <- matrix(runif(m*n, 0, 1), n, m)

#X <- matrix(0, n, m)
Xb <- X%*%b
cbind(id=1:n,Xb)[which.max(Xb),]

dt <- generate_data(
    intercept = int,
    beta = b,
    X = X,
    q = q,
    rho = rho,
    seed = 1
)

theta <- c(eps, rho, b, int)
repar_theta <- partorepar(theta)
# check marginal of latent variables
C <- generate_C(rho = rho, p = p)
s <- 123
Z <- t(sapply(1:n, function(x) generate_mgamma(q, C, seed = s + x)))
Z2 <- matrix(rep(Z[1,], n), nrow = n, ncol =  p, byrow = T)

gg_density <- as_tibble(Z) %>%
    gather(key = 'response', value = 'frailty', starts_with('V')) %>%
    mutate(response = as.factor(response)) %>%
    ggplot()+
    geom_density(aes(x = frailty, group = response, col = response), alpha = .5)+
    geom_line(
        data = tibble(frailty = seq(0.01,3, by = 0.1)) %>%
            mutate(density = dgamma(frailty, 1/eps, 1/eps)),
        aes(x = frailty, y = density), linetype = 'dashed', col = 'black'
    )+
    theme_minimal()+
    scale_color_viridis_d()
plotly::ggplotly(gg_density, dynamicTicks = T)
(cor(Z) - C^2)/cor(Z)

# check marginals count responses
u <- t(sapply(1:n, function(i) exp(int + as.numeric(crossprod(b, X[i,])))))
tib <- tibble(response = as.factor(paste0('V', 1:p)), mu = u[1,], Z = Z[1,]) %>%
    expand_grid(count = 0:max(dt)) %>%
    mutate(
        density = dpois(x = count, lambda = Z*mu)
    )
tib$density3 <- map2_dbl(tib$mu,tib$count, ~uni_dnbinom(.x, .y, eps = eps))
gg_density <- as_tibble(dt) %>%
    gather(key = 'response', value = 'count', starts_with('V')) %>%
    mutate(response = as.factor(response)) %>%
    ggplot()+
    geom_density(aes(x = count, group = response, col = response), alpha = .5)+
    geom_line(
        data = tib,
        aes(x = count, y = density3), linetype = 'dashed', col = 'black'
    )+
    facet_wrap(vars(response), scales = 'free')+
    theme_minimal()+
    scale_color_viridis_d()

plotly::ggplotly(gg_density, dynamicTicks = T)


#### Derivatives #####
dt

i <- 1
j <- 0
jp <- 1
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


system.time(
    opt <- ucminf::ucminf(par_init, Rwrapper_ncl, Rwrapper_ngr)
)
repartopar(opt$par)
theta

mean((repartopar(opt$par)-theta)^2)
mean((repartopar(par_init)-theta)^2)

#root <- nleqslv::nleqslv(par_init, fn = Rwrapper_ngr)

maxT <- n; set.seed(1)
trajSub <- c(50, 100, 1000, 2000)
trajSub0 <- c(1, trajSub[trajSub < maxT]+1, maxT+1)
system.time(
    Opt <- gammaFrailty(
        THETA_INIT = par_init,
        DATA = dt,
        X = X,
        MAXT = maxT,
        BURN = maxT*(2/3),
        STEPSIZE = .01,
        STEPSIZE0 = .001,
        NU = 1,
        METHODFLAG = 1,
        VERBOSEFLAG = F,
        STEPSIZEFLAG = F
    )
)
#Opt
vec <- Opt$path_theta[c(51,101),]
Opt$path_theta[1:4,]
mean((repartopar(Opt$path_av_theta[maxT+1,])-theta)^2)
mean((repartopar(opt$par)-theta)^2)
mean((repartopar(par_init)-theta)^2)

probl <- Opt$path_theta[9,]
repartopar(probl)
ncl(probl, dt, X, T)
i <- 249; j <- 10; jp <- 6
n_j  <- dt[i + 1, j+1 ]; n_jp <- dt[i + 1, jp+1]
x <- X[i+1,];
Rwrapper_obj(probl)
Rwrapper_obj_ind(probl, ind = 0)

J <- sampleJ(repar_theta, dt, X)
H <- sampleH(par_init, dt, X, invertFLAG = T)
diag(H)
Hnum <- numDeriv::jacobian(Rwrapper_ngr, par_init)
Hnum-H
diag(solve(Hnum))
diag(solve(H))

par_init <- repar_theta + runif(length(repar_theta), -1, 1);
maxT <- n
ctrl <- list(
    MAXT = maxT,
    BURN = 0,
    STEPSIZE = .01,
    STEPSIZE0 = .005,
    NU = 1,
    SEED = 123
)
ctrl <- list(
    MAXT = n^1.5,
    BURN = 100,
    STEPSIZE = .001,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = 123
)

trajSub <- c(500, 1000, 2000)
ctrl <- list()
trajSub <- c(50, 100, 1000, 2000)
Opt0 <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SGD',
    CPP_CONTROL = ctrl,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
Opt1$fit$path_theta[trajSub0,]
Opt0$fit$path_theta

Opt0$time
mean((repartopar(Opt0$theta)-theta)^2)
mean((repartopar(Opt0$theta_init)-theta)^2)
Opt0$control
Opt0$method
Opt0$fit
repartopar(Opt0$theta)
repartopar(par_init)

theta

start_time <- Sys.time()

end_time <- Sys.time()

d <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
d
x <- matrix(rnorm(n*m), n, m)
x <- X
varList <- sampleVar(
    THETA = Opt0$theta,
    DATA = dt,
    X = x,
    NU = 1,
    METHOD = 2,
    RANGE = Opt0$control$MAXT - Opt0$control$BURN,
    TOTFLAG = TRUE,
    PRINTFLAG = F
    )

diag(varList$var$var_stoc*Opt0$control$MAXT/n)
diag(varList$var$var_stat)

diag(varList$var$var_stoc*Opt0$control$MAXT*((1/Opt0$control$MAXT) + (1/n)))
diag(varList$var$var_tot)

##################
#### libraries #####
library(tidyverse)
###### choose true model

p <- 50
q <- 8

xi <- 2/q
rho <- .6

m <- 2
n <- 2000
int <- rep(log(4), p)#runif(p, 0, 1)
b <- rep(0, m) #rnorm(m, 0, .5) #
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)

##### generate the data ###
dt <- generate_data(
    intercept = int,
    beta = b,
    X = X,
    q = q,
    rho = rho,
    seed = 1
)

##### estimation ####
par_init <- repar_theta + runif(length(repar_theta), -1, 1)


ctrl_scsd <- list(
    MAXT = n,
    BURN = 100,
    STEPSIZE = .001,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = 123
)


ctrl_gd <- list(MAXT = n^.75,     STEPSIZE = .001)
trajSub <- c(50, 100, 500, 1000)

mean((repartopar(Opt0$theta_init)-theta)^2)

Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = ctrl,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
mean((repartopar(Opt_u$theta)-theta)^2)

Opt_sgd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SGD',
    CPP_CONTROL = ctrl_scsd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
mean((repartopar(Opt_sgd$theta)-theta)^2)

Opt_scsd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SCSD',
    CPP_CONTROL = ctrl_scsd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
mean((repartopar(Opt_scsd$theta)-theta)^2)

Opt_gd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'GD',
    CPP_CONTROL = ctrl_gd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
mean((repartopar(Opt_gd$theta)-theta)^2)
