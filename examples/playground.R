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
Z <- t(sapply(1:n, function(x) generate_mgamma(q, C, SEED = s + x)))
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
#### choose true model ####

p <- 12
q <- 8

xi <- 2/q
rho <- .6

m <- 10
n <- 1000
int <- rep(log(4), p)#runif(p, 0, 1)
b <- rep(0, m) #rnorm(m, 0, .5) #
set.seed(1)
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)

##### generate the data ###
dt <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = 1
)

##### estimation ####
set.seed(1); par_init <- repar_theta + runif(length(repar_theta), -1, 1)


ctrl_scsd <- list(
    MAXT = 5000,
    BURN = 1000,
    STEPSIZE = .002,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = 1
)


ctrl_gd <- list(MAXT = n^.75,     STEPSIZE = .001)
trajSub <- c(50, 100, 500, 1000)
trajSub <- c(750, 1000, 2000, 3000, 5000)

mean((repartopar(par_init)-theta)^2)

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
    ITERATIONS_SUBSET = NULL#trajSub
)
mean((repartopar(Opt_sgd$theta)-theta)^2)
mean((Opt_sgd$theta-repar_theta)^2)

Opt_scsd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SCSD',
    CPP_CONTROL = ctrl_scsd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
mean((repartopar(Opt_scsd$theta)-theta)^2)
mean((Opt_scsd$theta-repar_theta)^2)

Opt_gd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'GD',
    CPP_CONTROL = ctrl_gd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = trajSub
)
mean((repartopar(Opt_gd$theta)-theta)^2)

sampleVar(
    THETA = Opt_sgd$theta,
    DATA = dt,
    X = X,
    NU = 1,
    METHOD = 1,
    RANGE = ctrl_scsd$MAXT-ctrl_scsd$BURN,
    TOTFLAG = T,
    PRINTFLAG = F
)

#### check hessian ####
Rwrapper_ncl <- function(par){
    ncl(par, dt, X)$nll
}
Rwrapper_ngr <- function(par){
    ncl(par, dt, X)$ngradient
}

chosen_par <- Opt_sgd$fit$path_av_theta[1001,]#Opt_u$theta#par_init #
chosen_par <- par_init
chosen_par <- Opt_u$theta

H <- sampleH(chosen_par, dt, X, INVERTFLAG = F)
Hnum <- numDeriv::jacobian(Rwrapper_ngr, chosen_par)
H2 <- H
diag(H2) <- diag(Hnum)
diag(H)
diag(Hnum)
Hnum %>% solve() %>% diag() %>% round(4)
H2 %>% solve() %>% diag() %>% round(4)
H %>% solve() %>% diag() %>% round(4)

it <- 5000
sgd_se <- sampleVar(
    THETA = Opt_sgd$fit$path_av_theta[it,],
    DATA = dt,
    X = X,
    NU = 1,
    METHOD = 1,
    RANGE = it-ctrl_scsd$BURN,
    TOTFLAG = T,
    PRINTFLAG = F
)

scsd_se <- sampleVar(
    THETA = Opt_scsd$fit$path_av_theta[it,],
    DATA = dt,
    X = X,
    NU = 1,
    METHOD = 2,
    RANGE = it-ctrl_scsd$BURN,
    TOTFLAG = T,
    PRINTFLAG = F
)
sgd_se$se$se_stoc
scsd_se$se$se_stoc

Hinv <- sampleH(Opt_scsd$fit$path_av_theta[it,], dt, X, INVERTFLAG = T)
J <- sampleJ(Opt_scsd$fit$path_av_theta[it,], dt, X)
(Hinv%*%J%*%Hinv/(it-ctrl_scsd$BURN)) %>% diag() %>% sqrt()
(Hinv/(it-ctrl_scsd$BURN)) %>% diag() %>% sqrt()
scsd_se$var$var_stoc %>% diag() %>% sqrt()
H %>%  diag(); (scsd_se$var$var_stoc * (it-ctrl_scsd$BURN)) %>% diag()
diag(J);diag(Hinv); diag(solve(Hinv))
gg <- get_tidy_path(Opt_sgd, 'path_av_theta', F) %>%
    mutate( mod = 'SGD') %>%
    bind_rows(
        get_tidy_path(Opt_scsd, 'path_av_theta', F) %>%
            mutate( mod = 'SCSD')
    ) %>%
    mutate(
        mse = map_dbl(path_av_theta, ~mean((.x-repar_theta)^2))
    ) %>%
    ggplot( aes(x = iter, y = (mse), col = mod))+
    geom_line()+
    geom_hline(yintercept = (mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    theme_minimal()
plotly::ggplotly(gg)

vec <- 1:10
{
    start_time2 <- Sys.time()                                                 # Save starting time
    ls2 <- pbapply::pblapply(vec, function(id){

        # INTERCEPT <- int
        # BETA <- b
        # Q <- q
        # RHO <- rho
        # SEED <-id
        # C <- generate_C(rho = RHO, p = p)
        # set.seed(SEED)
        # rmvn(SAMPLE_SIZE = Q, VAR = C)
        #MASS::mvrnorm(n = Q, mu = rep(0, ncol(C)), Sigma = C)
        #mvtnorm::rmvnorm(n = Q, sigma = C)
        #out <- colMeans((rmvnorm(n = Q, sigma = C))^2)
        #Z <- purrr::reduce(purrr::map(1:n, ~generate_mgamma(Q, C, seed = SEED + .x)), rbind)


        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X[1:250,],
            Q = q,
            RHO = rho,
            SEED = id
        )

        suppressMessages(
            mod_obj <- fit_gammaFrailty(
                DATA_LIST = list('DATA' = data, 'X' = X[1:250,]),
                METHOD = 'SGD',
                CPP_CONTROL = list(
                    MAXT = 1000,
                    BURN = 100,
                    STEPSIZE = .001,
                    NU = 1,
                    SEED = id
                ),
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = NULL,
            )
        )

        return(C)
    }, cl = 10)
    end_time2 <- Sys.time()
    time_diff2 <- end_time2 - start_time2
    time_diff2

}

a <- t(sapply(1:n, function(x) generate_mgamma(q, C = generate_C(rho = rho, p = p), seed = 1 + x)))


 sapply(1:n, function(x) generate_mgamma(q, C = generate_C(rho = rho, p = p), seed = 1 + x))



generate_mgamma(q, C = generate_C(rho = rho, p = p), seed = 1 + 2)
b <- reduce(map(1:n, ~generate_mgamma(q, C = generate_C(rho = rho, p = p), seed = 1 + .x)), rbind)


a[1,]
b[1,]

c <- generate_C(rho = rho, p = p)
a <- matrix(rnorm(q*p), p, q)
t(t(chol(c))%*%a)

dim <- 2
sample_size <- 50000
Var <- matrix(.3, dim, dim); diag(Var) <- 1
sample <- t(t(chol(Var))%*%matrix(rnorm(dim*sample_size), dim, sample_size))
cor(sample)
sample0 <- mvtnorm::rmvnorm(sample_size, sigma = Var)
cor(sample0)



sample <- rmvn(sample_size, Var)
cor(sample)
