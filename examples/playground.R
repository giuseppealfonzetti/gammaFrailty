#### libraries #####
library(tidyverse)
######

p <- 12
q <- 8

eps <- 2/q
rho <- .6

m <- 2
n <- 100
int <- rep(log(4), p)#runif(p, 0, 1)
b <- rep(0, m) #rnorm(m, 0, .5) #
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
#X <- matrix(runif(m*n, 0, 1), n, m)

#X <- matrix(0, n, m)
Xb <- X%*%b
cbind(id=1:n,Xb)[which.max(Xb),]

dt <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = seed
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
H <- sampleH(repar_theta, dt, X, INVERTFLAG = T)
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
int <- runif(p, -.5, .5)#rep(.5, p)#
b <- rnorm(m, 0, .5) #runif(m,0,.5)#rep(0, m) #
set.seed(1)
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)

seed <- 3
##### generate the data ###
dt <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = seed
)

##### estimation ####
set.seed(1); par_init <- repar_theta + runif(length(repar_theta), -1, 1)




ctrl_gd <- list(MAXT = n^.75,     STEPSIZE = .001)
trajSub <- c(50, 100, 500, 1000)
trajSub <- c(750, 1000, 2000, 3000, 5000, 7500, 10000)

mean((par_init-repar_theta)^2)
H0 <- sampleH(THETA = par_init, DATA = dt, X = X, F, F)
Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
mean((Opt_u$theta-repar_theta)^2)

sgd_var <- sampleVar(
    THETA = repar_theta,
    DATA = dt,
    X = X,
    NU = 1,
    METHOD = 2,
    RANGE = ctrl_scsd$MAXT-ctrl_scsd$BURN,
    TOTFLAG = T,
    PRINTFLAG = F
)
eig <- eigen(solve(sgd_var$var$var_stoc))$values
 2/max(eig) ;1/(2*min(eig));
ctrl_sgd <- list(
    MAXT = 1000,
    BURN = 200,
    STEPSIZE = 0.0001265499,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = seed
)
Opt_sgd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SGD',
    CPP_CONTROL = ctrl_sgd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL#trajSub
)
#mean((repartopar(Opt_sgd$theta)-theta)^2)
mean((Opt_sgd$theta-repar_theta)^2)
#Opt_sgd$fit$path_av_theta
dim(Opt_sgd$fit$path_theta)
length(Opt_sgd$fit$path_nll)
Rwrapper_ncl(Opt_sgd$theta_init)
Rwrapper_ncl(Opt_sgd$theta)

ctrl_scsd <- list(
    MAXT = 1000,
    BURN = 200,
    STEPSIZE = 5e-4,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = seed
)
Opt_scsd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SCSD',
    CPP_CONTROL = ctrl_scsd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
#mean((repartopar(Opt_scsd$theta)-theta)^2)
mean((Opt_scsd$theta-repar_theta)^2)
Opt_scsd$fit


# test <- Opt
# #ncl(par_init, dt, X, T)
# Opt_u$theta
# test$theta
# sum(Opt_scsd$fit$weights[[2]]!=0)
# sum(Opt_scsd$fit$weights[[2]] == Opt_scsd$fit$weights[[3]])
# which(Opt_scsd$fit$weights[[2]]!=0, arr.ind = T)
# which(Opt_scsd$fit$weights[[4]]!=0, arr.ind = T)

# Opt_gd <- fit_gammaFrailty(
#     DATA_LIST = list('DATA' = dt, 'X' = X),
#     METHOD = 'GD',
#     CPP_CONTROL = ctrl_gd,
#     VERBOSEFLAG= 0,
#     INIT = par_init,
#     ITERATIONS_SUBSET = trajSub
# )
# mean((repartopar(Opt_gd$theta)-theta)^2)

# sampleVar(
#     THETA = Opt_sgd$theta,
#     DATA = dt,
#     X = X,
#     NU = 1,
#     METHOD = 1,
#     RANGE = ctrl_scsd$MAXT-ctrl_scsd$BURN,
#     TOTFLAG = T,
#     PRINTFLAG = F
# )
#
# #### check hessian ####
# Rwrapper_ncl <- function(par){
#     ncl(par, dt, X)$nll
# }
# Rwrapper_ngr <- function(par){
#     ncl(par, dt, X)$ngradient
# }
#
# chosen_par <- Opt_sgd$fit$path_av_theta[1001,]#Opt_u$theta#par_init #
# chosen_par <- par_init
# chosen_par <- Opt_u$theta
#
# H <- sampleH(chosen_par, dt, X, INVERTFLAG = F)
# Hnum <- numDeriv::jacobian(Rwrapper_ngr, chosen_par)
# H2 <- H
# diag(H2) <- diag(Hnum)
# diag(H)
# diag(Hnum)
# Hnum %>% solve() %>% diag() %>% round(4)
# H2 %>% solve() %>% diag() %>% round(4)
# H %>% solve() %>% diag() %>% round(4)
#
# it <- 5000
# sgd_se <- sampleVar(
#     THETA = Opt_sgd$fit$path_av_theta[it,],
#     DATA = dt,
#     X = X,
#     NU = 1,
#     METHOD = 1,
#     RANGE = it-ctrl_scsd$BURN,
#     TOTFLAG = T,
#     PRINTFLAG = F
# )
#
# scsd_se <- sampleVar(
#     THETA = Opt_scsd$fit$path_av_theta[it,],
#     DATA = dt,
#     X = X,
#     NU = 1,
#     METHOD = 2,
#     RANGE = it-ctrl_scsd$BURN,
#     TOTFLAG = T,
#     PRINTFLAG = F
# )
# sgd_se$se$se_stoc
# scsd_se$se$se_stoc
#
# Hinv <- sampleH(Opt_scsd$fit$path_av_theta[it,], dt, X, INVERTFLAG = T)
# J <- sampleJ(Opt_scsd$fit$path_av_theta[it,], dt, X)
# (Hinv%*%J%*%Hinv/(it-ctrl_scsd$BURN)) %>% diag() %>% sqrt()
# (Hinv/(it-ctrl_scsd$BURN)) %>% diag() %>% sqrt()
# scsd_se$var$var_stoc %>% diag() %>% sqrt()
# H %>%  diag(); (scsd_se$var$var_stoc * (it-ctrl_scsd$BURN)) %>% diag()
# diag(J);diag(Hinv); diag(solve(Hinv))
lab <- 'path_av_theta'
gg <- get_tidy_path(Opt_sgd, lab, F) %>%
    mutate( mod = 'SGD') %>%
    bind_rows(
        get_tidy_path(Opt_scsd, lab, F) %>%
            mutate( mod = 'SCSD')
    ) %>%
    mutate(
        mse = map_dbl(path_av_theta, ~mean((.x-repar_theta)^2))
    ) %>%
    ggplot( aes(x = iter, y = log(mse), col = mod))+
    geom_line()+
    geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    theme_bw()+
    scale_color_viridis_d()
plotly::ggplotly(gg)

lab <- 'path_nll'
get_tidy_path(Opt_sgd, lab, F) %>%
    mutate( mod = 'SGD') %>%
    bind_rows(
        get_tidy_path(Opt_scsd, lab, F) %>%
            mutate( mod = 'SCSD')
    )%>%
    ggplot( aes(x = iter, y = path_nll, col = mod))+
    geom_line()+
    #geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    theme_bw()+
    scale_color_viridis_d()

lab <- 'path_grad'
get_tidy_path(Opt_sgd, lab, F) %>%
    mutate( mod = 'SGD') %>%
    bind_rows(
        get_tidy_path(Opt_scsd, lab, F) %>%
            mutate( mod = 'SCSD')
    )%>%
    mutate(
        grad_norm = map_dbl(path_grad, ~norm(as.matrix(.x)))
    ) %>%
    ggplot( aes(x = iter, y = grad_norm, col = mod))+
    geom_line()+
    #geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    theme_bw()+
    scale_color_viridis_d()

get_tidy_path(Opt_sgd, lab, F) %>%
    mutate( mod = 'SGD') %>%
    bind_rows(
        get_tidy_path(Opt_scsd, lab, F) %>%
            mutate( mod = 'SCSD')
    )%>%
    mutate(
        grad_norm = map_dbl(path_grad, ~norm(as.matrix(.x)))
    ) %>%
    filter(iter == 1)




#############
ctrl_sgd <- list(
    MAXT = 1000,
    BURN = 200,
    STEPSIZE = 0.0001265499,
    #STEPSIZE0 = .0005,
    NU = 1,
    SEED = seed
)
Opt_sgd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SGD',
    CPP_CONTROL = ctrl_sgd,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL#trajSub
)
#mean((repartopar(Opt_sgd$theta)-theta)^2)
mean((Opt_sgd$theta-repar_theta)^2)
#Opt_sgd$fit$path_av_theta
dim(Opt_sgd$fit$path_theta)
length(Opt_sgd$fit$path_nll)




######## simulation test #######
set.seed(1); par_init <- repar_theta + runif(length(repar_theta), -1, 1)
Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
library(tidyverse)
sim_settings <- expand_grid(
    mod = c('SGD', 'SCSD'),
    stepsize = c(5e-5,1.25e-4),
    stoc_seed = 1:5,
    maxiter = 5000,
    burn = 500
)

custom_est_fun <- function(MOD, STEPSIZE, SEED, MAXT, BURN){
    ctrl <- list(
        MAXT = MAXT,
        BURN = BURN,
        STEPSIZE = STEPSIZE,
        NU = 1,
        SEED = SEED
    )
    mod_obj <- fit_gammaFrailty(
        DATA_LIST = list('DATA' = dt, 'X' = X),
        METHOD = MOD,
        CPP_CONTROL = ctrl,
        VERBOSEFLAG= 0,
        INIT = par_init,
        ITERATIONS_SUBSET = NULL#trajSub
    )

    return(mod_obj)
}

#custom_est_fun(MOD = 'SGD', STEPSIZE=1.25e-4, SEED = 123, MAXT = 1000, BURN = 200)
est_tab <- sim_settings %>%
    mutate(
        mod_obj = pmap(
            list(mod, stepsize, stoc_seed, maxiter, burn),
            function(mod_, stepsize_, stoc_seed_, maxiter_, burn_){
                custom_est_fun(
                    MOD = mod_,
                    STEPSIZE = stepsize_,
                    SEED = stoc_seed_,
                    MAXT = maxiter_,
                    BURN = burn_)
            }
            )
    )


metrics_tab <- est_tab %>%
    mutate(
        path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta', F)),
        path_nll = map(mod_obj, ~get_tidy_path(.x, 'path_nll', F)),
        path_grad = map(mod_obj, ~get_tidy_path(.x, 'path_grad', F))
    ) %>%
    select(-mod_obj) %>%
    mutate(
        metrics = pmap(
            list(path_av_theta, path_nll, path_grad),
            function(path_av_theta_, path_nll_, path_grad_){
                path_av_theta_ %>%
                    left_join(path_nll_, by = 'iter') %>%
                    left_join(path_grad_, by = 'iter')
            }
        )
    ) %>%
    select(-c(path_av_theta, path_nll, path_grad)) %>%
    unnest(c(metrics)) %>%
    mutate(
        mse = map_dbl(path_av_theta, ~mean((.x-repar_theta)^2)),
        grad_norm = map_dbl(path_grad, ~norm(as.matrix(.x)))
    ) %>%
    gather(key = 'performance', value = 'val', path_nll, grad_norm, mse)

gg <- metrics_tab %>%
    ggplot( aes(x = iter, y = val, col = factor(stepsize), group = interaction(mod, stepsize, stoc_seed)))+
    geom_line(aes(linetype = mod))+
    #geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    facet_wrap(vars(performance), scales = 'free')+
    theme_bw()+
    scale_color_viridis_d()
plotly::ggplotly(gg, dynamicTicks = T)

########
name_par <- function(par){
    if(par == 1)
        'xi'
    else if(par == 2)
        'correlation'
    else if(par %in% 3:(2+m))
        'regression coefficients'
    else
        'intercepts'
}
true_tib <- tibble(
    par = 1:length(repar_theta), true_val = repar_theta
) %>%
    mutate(
    par_type = map_chr(par, ~name_par(.x)),
    par = as.factor(par))
    # ) %>%
    # mutate(
    #     true_val = map2_dbl(par, true_val, ~if(.x==2){rofz_cpp(.y)}else{.y}),
    #     true_val = map2_dbl(par, true_val, ~if(.x==1){exp(-.y)}else{.y})
    #     )

num_tib <- tibble(
    par = 1:length(repar_theta), num_val = Opt_u$theta
)%>%
    mutate(
        par_type = map_chr(par, ~name_par(.x)),
        par = as.factor(par))
    # ) %>%
    # mutate(
    #     num_val = map2_dbl(par, num_val, ~if(.x==2){rofz_cpp(.y)}else{.y}),
    #     num_val = map2_dbl(par, num_val, ~if(.x==1){exp(-.y)}else{.y})
    # )

par_tab <- est_tab %>%
    mutate(
        path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta', F))
    ) %>%
    unnest(c(path_av_theta)) %>%
    mutate(
        path_av_theta = lapply(path_av_theta, function(x){
            tib <- tibble(
                par = 1:length(x),
                val = x
            )
            tib
        })
    ) %>%
    unnest(c(path_av_theta)) %>%
    mutate(
        par_type = map_chr(par, ~name_par(.x)),
        par = as.factor(par)
    )


gg <- par_tab  %>%
    # mutate(
    #     val = map2_dbl(par, val, ~if(.x==2){rofz_cpp(.y)}else{.y}),
    #     val = map2_dbl(par, val, ~if(.x==1){exp(-.y)}else{.y})
    # ) %>%
    ggplot(aes(x = iter, y = val))+
    geom_line(aes(linetype = mod,  col = factor(stepsize), group = interaction(mod, stepsize, stoc_seed, par))) +
    geom_point(data = num_tib, aes(x = 5000, y = num_val), col = 'red', shape = 4, size = 2)+
    geom_point(data = true_tib, aes(x = 5020, y = true_val), col = 'blue', shape = 4, size = 2)+
    facet_wrap(vars(par_type), scales = 'free') +
    theme_bw()+
    scale_color_viridis_d()
gg
plotly::ggplotly(gg, dynamicTicks = T)


av_par_tab <- est_tab %>%
    mutate(
        path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta', F))
    ) %>%
    unnest(c(path_av_theta)) %>%
    mutate(
        path_av_theta = lapply(path_av_theta, function(x){
            tib <- tibble(
                par = 1:length(x),
                val = x
            )
            tib
        })
    ) %>%
    unnest(c(path_av_theta)) %>%
    select(mod, stepsize, iter, par, val) %>%
    group_by(mod, stepsize, iter, par) %>%
        summarise(av_val = mean(val)) %>%
    mutate(
        par_type = map_chr(par, ~name_par(.x)),
        par = as.factor(par)
    )

gg1 <- av_par_tab  %>%
    filter(iter %in% seq(0, 5000, 100)) %>%
    mutate(av_val = map2_dbl(av_val, par_type, ~if_else(.y == 'correlation', rofz_cpp(.x), .x))) %>%
    ggplot(aes(x = iter, y = av_val))+
    geom_line(aes(linetype = mod,  col = factor(stepsize), group = interaction(mod, stepsize, par))) +
    geom_point(data = num_tib %>% mutate(num_val = map2_dbl(num_val, par_type, ~if_else(.y == 'correlation', rofz_cpp(.x), .x)))
, aes(x = 5000, y = num_val), col = 'red', shape = 4, size = 2)+
    geom_point(data = true_tib, aes(x = 5020, y = true_val), col = 'blue', shape = 4, size = 2)+
    facet_wrap(vars(par_type), scales = 'free') +
    theme_bw()+
    scale_color_viridis_d()

gg1
plotly::ggplotly(gg1, dynamicTicks = T)
