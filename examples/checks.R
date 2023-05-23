#### libraries #####
library(tidyverse)

#### choose true model ####

p <- 20
q <- 4

xi <- 2/q
rho <- .4

m <- 5
n <- 500
int <- rnorm(p, 0, .05)#rep(0, p)#
b <- rnorm(m, 0, .05) #rep(0, m)#
set.seed(1)
X <- matrix(rbinom(m*n, 1, .5), n, m)#matrix(runif(m*n, 0, 1), n, m)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)
d <- length(theta)

##### generate the data ###
seed <- 3
structlab <- 'AR'
dt <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = seed,
    STRUCT = structlab
)

#### Derivatives #####
dt
structlabb <- 0
i <- 13
j <- 7
jp <- 2
n_j <- dt[i + 1, j+1]
n_j <- 19
n_jp <- 30
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
# sm <- sapply(0:min(n_j, n_jp), function(x)Rwrapper_obj_ind(repar_theta, x)$Sind)
# sm1 <- sm
# sm[1]+sm[length(sm)]
# sm1[1] <- 0
# cumsum(sm1)
# cumsum(sm)
# Rwrapper_obj_ind(repar_theta, 0)$Sind
# Rwrapper_obj_ind(repar_theta, 1)$Sind
# Rwrapper_obj_ind(repar_theta, 2)$Sind
# Rwrapper_obj_ind(repar_theta, 3)$Sind
# Rwrapper_obj_ind(repar_theta, 4)$Sind



Rwrapper_obj <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2], verboseS = T, STRUCT = structlabb
    )

    obj
}
structlabb <- 1
Rwrapper_obj(repar_theta)


Rwrapper_obj_fun <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2], STRUCT = structlabb
    )

    obj$ll
}
Rwrapper_obj_fun_stable <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2]
    )

    obj$ll_stable
}
Rwrapper_obj_der <- function(par){
    obj <- pair_wrapper(
        j, jp, n_j, n_jp, p, alpha_j = par[m+2+j+1], alpha_jp= par[m+2+jp+1], x, beta = par[3:(m+2)], lxi = par[1], artanhrho = par[2], STRUCT = structlabb
    )

    obj$gradient
}

#### Derivatives cl ####
Rwrapper_obj_fun(repar_theta)
Rwrapper_obj_fun_stable(repar_theta)
Rwrapper_obj_fun(repar_theta)
numDeriv::grad(Rwrapper_obj_fun, repar_theta)
Rwrapper_obj_der(repar_theta)
p_range = 100
ncl(theta, dt, X, T)
Rwrapper_ncl <- function(par){
    ncl(par, dt, X, PAIRS_RANGE = p_range, STRUCT = structlabb)$nll
}

Rwrapper_ngr <- function(par){
    ncl(par, dt, X, PAIRS_RANGE = p_range, STRUCT = structlabb)$ngradient
}

par_init <- rep(0, length(theta))
par_init <- repar_theta + runif(length(repar_theta), -1, 1);
structlabb <- 1

Rwrapper_ncl(par = par_init)
numDeriv::grad(Rwrapper_ncl, par_init)
Rwrapper_ngr(par = par_init)

numH <- numDeriv::jacobian(Rwrapper_ngr, repar_theta)
H <- sampleH(repar_theta, dt, X, INVERTFLAG = F, F, PAIRS_RANGE = p_range)
diag(numH)-diag(H)
diag(solve(numH))
#### Opt ####
# p <- 9
# p_range <- 10
# counter <- 0
# for (i in 1:(p-1)) {
#     for (j in max(0, i-p_range):(i-1)) {
#          counter <- counter+1
#     }
# }
# counter
#
# (p-p_range)*p_range+(p_range*(p_range-1))/2
par_init <- repar_theta + runif(length(repar_theta), -1, 1);
p_range <- 3
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

scls <- diag(sampleH(repar_theta, dt, X, INVERTFLAG = T, PAIRS_RANGE = range_lag))
ctrl_sgd <- list(
    MAXT = 5000,
    BURN = 1000,
    STEPSIZE = .01,
    #SCALEVEC = scls,
    NU = 1,
    SEED = seed
)
Opt_sgd <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'SCSD',
    CPP_CONTROL = ctrl_sgd,
    PAIRS_RANGE = 3,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL#trajSub
)
repartopar(Opt_sgd$theta)
theta
mean(sum(Opt_sgd$theta-repar_theta)^2)

     Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    PAIRS_RANGE = 5,
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

######## simulation test #######
set.seed(1); par_init <- repar_theta + runif(length(repar_theta), -1, 1)
par_init <- rep(0, d)
range_lag = 100
Hinv <- sampleH(repar_theta, dt, X, INVERTFLAG = T, PAIRS_RANGE = range_lag)
scls <- diag(Hinv)
scls <- rep(.1, d)
structlab <- 'COMPOUND'
Opt_u <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = dt, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list('ctrl' = list(invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]), 'hessian' = 0),
    PAIRS_RANGE = range_lag,
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL,
    STRUCT = structlab
)
library(tidyverse)
sim_settings <- expand_grid(
    mod = c('SGD', 'SCSD'),
    stepsize = c(1e-4, 5e-4, 1e-3, 5e-3,1e-2, 5e-2, 1e-1, 5e-1, 1),
    stoc_seed = 1:5,
    maxiter = 3000,
    burn = 500
)

custom_est_fun <- function(MOD, STEPSIZE, SEED, MAXT, BURN){
    ctrl <- list(
        MAXT = MAXT,
        BURN = BURN,
        STEPSIZE = STEPSIZE,
        SCALEVEC = scls,
        par1 = 1,
        par2 = STEPSIZE,
        par3 = .501,
        NU = 1,
        SEED = SEED,
        STEPSIZEFLAG = 1
    )
    mod_obj <- fit_gammaFrailty(
        DATA_LIST = list('DATA' = dt, 'X' = X),
        METHOD = MOD,
        CPP_CONTROL = ctrl,
        PAIRS_RANGE = range_lag,
        VERBOSEFLAG= 0,
        INIT = par_init,
        ITERATIONS_SUBSET = seq(0, 5000, 100),
        STRUCT = structlab
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


# metrics_tab <- est_tab %>%
#     mutate(
#         path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta', F)),
#         path_nll = map(mod_obj, ~get_tidy_path(.x, 'path_nll', F)),
#         path_grad = map(mod_obj, ~get_tidy_path(.x, 'path_grad', F))
#     ) %>%
#     select(-mod_obj) %>%
#     mutate(
#         metrics = pmap(
#             list(path_av_theta, path_nll, path_grad),
#             function(path_av_theta_, path_nll_, path_grad_){
#                 path_av_theta_ %>%
#                     left_join(path_nll_, by = 'iter') %>%
#                     left_join(path_grad_, by = 'iter')
#             }
#         )
#     ) %>%
#     select(-c(path_av_theta, path_nll, path_grad)) %>%
#     unnest(c(metrics)) %>%
#     mutate(
#         mse = map_dbl(path_av_theta, ~mean((.x-repar_theta)^2)),
#         grad_norm = map_dbl(path_grad, ~norm(as.matrix(.x)))
#     ) %>%
#     gather(key = 'performance', value = 'val', path_nll, grad_norm, mse)
#
# gg <- metrics_tab %>%
#     ggplot( aes(x = iter, y = val, col = factor(stepsize), group = interaction(mod, stepsize, stoc_seed)))+
#     geom_line(aes(linetype = mod))+
#     #geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
#     facet_wrap(vars(performance), scales = 'free')+
#     theme_bw()+
#     scale_color_viridis_d()
# plotly::ggplotly(gg, dynamicTicks = T)

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
    par = 1:length(repar_theta), true_val = repar_theta) %>%
    mutate(par_type = map_chr(par, ~name_par(.x)),
           par = as.factor(par))


num_tib <- tibble(par = 1:length(repar_theta), num_val = Opt_u$theta)%>%
    mutate(par_type = map_chr(par, ~name_par(.x)),
           par = as.factor(par))

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
    #mutate(av_val = map2_dbl(av_val, par_type, ~if_else(.y == 'correlation', rofz_cpp(.x), .x))) %>%
    ggplot(aes(x = iter, y = av_val))+
    geom_line(aes(linetype = mod,  col = factor(stepsize), group = interaction(mod, stepsize, par))) +
    geom_point(data = num_tib
               , aes(x = 3020, y = num_val), col = 'red', shape = 4, size = 2)+
    geom_point(data = true_tib, aes(x = 3040, y = true_val), col = 'blue', shape = 4, size = 2)+
    facet_grid(par_type~stepsize, scales = 'free') +
    theme_bw() +
    scale_color_viridis_d()
plotly::ggplotly(gg1, dynamicTicks = T)

gg1

est_tab %>% pluck('mod_obj', 1) %>% pluck('actrl_args')
