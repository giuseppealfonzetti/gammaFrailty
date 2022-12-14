#'@export
fit_gammaFrailty <- function(
        DATA_LIST = list('DATA', 'X'),
        METHOD,
        CPP_CONTROL = list(
            MAXT = 1000,
            BURN = 500,
            STEPSIZE = .01,
            STEPSIZE0 = NULL,
            NU = 1,
            SEED = 123
        ),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
        INIT = NULL,
        ITERATIONS_SUBSET = NULL,
        VERBOSEFLAG = 0
){

    out <- list()
    start_time <- Sys.time()

    # Identify model dimensions
    p <- ncol(DATA_LIST$DATA)
    n <- nrow(DATA_LIST$DATA)
    m <- ncol(DATA_LIST$X)
    d <- p + m + 2

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
            out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at zero vector.')
            out$theta_init <-  rep(0, d)
    }




    # Check if method entry is correct
    if(!(METHOD %in% c('ucminf','GD', 'SGD', 'SCSD'))) stop('Method not available.')
    out$method <- METHOD


    # Numerical optimisation
    if(METHOD == 'ucminf'){

        message('2. Optimising with ucminf...')
        # R wrapper of cpp function for negative composite log-likelihood
        Rwr_ncl <- function(par){ ncl(par, DATA_LIST$DATA, DATA_LIST$X)$nll }

        # R wrapper of cpp function for negative composite score
        Rwr_ngr <- function(par){ ncl(par, DATA_LIST$DATA, DATA_LIST$X)$ngradient }

        # list of ucminf args
        args <- list(
            'par' = out$theta_init,
            'fn' = Rwr_ncl,
            'gr' = Rwr_ngr,
            'control' = UCMINF_CONTROL$ctrl,
            'hessian' = UCMINF_CONTROL$hessian)

        # optimisation
        opt <- do.call(ucminf::ucminf, args)

        out$control <- UCMINF_CONTROL
        out$theta   <- opt$par

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('3. Done! (', round(out$time,2),' secs)')
        return(out)
    }

    # Stochastic approximation of numerical optimiser
    if(METHOD == 'SGD' | METHOD == 'SCSD' ){

        message(paste0('2. Optimising with ', METHOD, '...'))

        # Check stochastic controlparameters
        cpp_ctrl <- check_SCSD_args(CPP_CONTROL, N = n)

         # Check iterations selected
        if(!is.null(ITERATIONS_SUBSET)){
            out$iterations_subset <- c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT)
        }else{
            out$iterations_subset <- 0:cpp_ctrl$MAXT
        }

        # Guarantee reproducibility stochastic optimisation
        set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

        args$SEED <- NULL

        args$METHODFLAG <- dplyr::if_else(METHOD == 'SGD', 1, 2)
        if(is.null(CPP_CONTROL$STEPSIZE0)){
            args$STEPSIZE0 <-  cpp_ctrl$STEPSIZE
            args$STEPSIZEFLAG <- 0
        } else {
            args$STEPSIZEFLAG <- 1
        }


        fit <- do.call(gammaFrailty, args)
        fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,]
        fit$path_av_theta <- fit$path_av_theta[out$iterations_subset + 1,]
        fit$path_grad     <- fit$path_grad[out$iterations_subset,]

        fit$methodflag <- NULL



        out$control <- cpp_ctrl
        out$fit <- fit
        out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]
        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('\n3. Done! (', round(out$time,2),' secs)')
        return(out)

    }

    # Gradient descent
    if(METHOD == 'GD'){

        message(paste0('2. Optimising with ', METHOD, '...'))

        cpp_ctrl <- check_GD_args(CPP_CONTROL, N = n)

        # Check iterations selected
        if(!is.null(ITERATIONS_SUBSET)){
            out$iterations_subset <- c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT)
        }else{
            out$iterations_subset <- 0:cpp_ctrl$MAXT
        }

        # Collect and rearrange arguments to pass to cpp function
        args <- list(
            'MAXT' = cpp_ctrl$MAXT,
            'METHODFLAG' = 0,
            'BURN' = 0,
            'STEPSIZE'= cpp_ctrl$STEPSIZE,
            'STEPSIZE0' = cpp_ctrl$STEPSIZE,
            'STEPSIZEFLAG' = F,
            'NU' = 1,
            'VERBOSEFLAG' = VERBOSEFLAG
        )
        args <- append(args, append(list( 'THETA_INIT' = out$theta_init), DATA_LIST))

        fit <- do.call(gammaFrailty, args)
        fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,   ]
        fit$path_grad     <- fit$path_grad[out$iterations_subset,  ]
        fit$path_av_theta <- NULL
        fit$methodflag <- NULL


        out$control <- cpp_ctrl
        out$fit <- fit
        out$theta <- fit$path_theta[nrow(fit$path_theta),]
        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('\n3. Done! (', round(out$time,2),' secs)')
        return(out)

    }
}
