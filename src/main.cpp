#include <Rcpp.h>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE

#include <RcppEigen.h>
#include <math.h>
#include <boost/math/special_functions.hpp>
#include "variance.h"
#include "pairClass.h"

//' @export
// [[Rcpp::export]]
Rcpp::List pair_wrapper(
        unsigned int j,
        unsigned int jp,          // periods indices
        unsigned int n_j,
        unsigned int n_jp,        // count values
        unsigned int p,           // number of periods
        double alpha_j,
        double alpha_jp,          // intercepts
        Eigen::VectorXd x,        // external regressors
        Eigen::VectorXd beta,     // regression coefficients
        double lxi,               // log xi = log inverse of gamma parameter
        double artanhrho,         // artanh of correlation of adjacent periods
        unsigned int ind = 0,
        bool verboseS = false,
        bool verboseSind = false

){


    pair_class pair;
    pair.setup(j, jp, n_j, n_jp, x, beta, alpha_j, alpha_jp, lxi, artanhrho, p);
    pair.compute_intermediate_( );
    pair.compute_dintermediate_(verboseS);

    Rcpp::List interm  = pair.get_intermediate_( );
    Rcpp::List dinterm = pair.get_dintermediate_();

    double ll = pair.compute_ll_();
    Eigen::VectorXd gradient = pair.compute_gradient_();


    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("Quantities") = interm,
        Rcpp::Named("Derivatives") = dinterm,
        Rcpp::Named("ll") = ll,
        Rcpp::Named("gradient") = gradient,
        Rcpp::Named("Sind") = pair.compute_Sind_(ind, verboseSind),
        Rcpp::Named("Qind") = pair.compute_Q_(ind, verboseSind)
        );

    // Rcpp::List output =
    //     Rcpp::List::create(
    //         Rcpp::Named("a") = 0
    //     );

    return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List ncl(
    Eigen::VectorXd theta,
    Eigen::MatrixXd data,
    Eigen::MatrixXd X,
    bool printFLAG = false
){
    unsigned int d = theta.size();
    unsigned int n = data.rows();
    unsigned int p = data.cols();
    unsigned int r = X.cols();

    double lxi  = theta(0);
    double artanhrho = theta(1);
    Eigen::VectorXd beta = theta.segment(2, r);



    pair_class pair;
    double nll = 0;
    Eigen::VectorXd ngradient = Eigen::VectorXd::Zero(d);
    for(unsigned int i = 0; i < n; i++){
        Eigen::VectorXd x_i = X.row(i);

        for(unsigned int j = 1; j < p; j++){
            unsigned int n_j = data(i, j);
            double alpha_j = theta(r+2+j);

            for( unsigned int jp = 0; jp < j; jp++){
                unsigned int n_jp = data(i, jp);
                double alpha_jp = theta(r+2+jp);

                pair.setup(j, jp, n_j, n_jp, x_i, beta, alpha_j, alpha_jp, lxi, artanhrho, p);
                pair.compute_intermediate_( );
                pair.compute_dintermediate_();
                double ll = pair.compute_ll_()/n;
                nll -= ll;
                ngradient -= pair.compute_gradient_()/n;

                if(printFLAG &!(ll==ll)) Rcpp::Rcout << "(i, j, j') = ("<<i<<", "<<j<<", "<< jp << ") -> "<< ll << "\n";

            }
        }


    }

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("nll") = nll,
        Rcpp::Named("ngradient") = ngradient
    );

    return output;
}




//' @export
// [[Rcpp::export]]
Rcpp::List gammaFrailty(
        Eigen::Map<Eigen::VectorXd> THETA_INIT,
        Eigen::Map<Eigen::MatrixXd> DATA,
        const Eigen::MatrixXd X,
        const unsigned int MAXT,
        const unsigned int BURN,
        const double STEPSIZE,
        const double STEPSIZE0,
        const double NU,
        const int METHODFLAG = 0,
        const bool VERBOSEFLAG = false,
        const bool STEPSIZEFLAG = false
){

    // Set up clock monitor to export to R session trough RcppClock
    // Rcpp::Clock clock;
    // clock.tick("main");

    // Identify model dimensions
    const unsigned int d = THETA_INIT.size();
    const unsigned int n = DATA.rows();
    const unsigned int p = DATA.cols();
    const unsigned int r = X.cols();
    const unsigned int kk = p*(p-1)/2;


    Eigen::VectorXd stepsizeV(d); stepsizeV.fill(STEPSIZE); stepsizeV(0) = STEPSIZE0;
    // Initialize storage for iterations quantities
    Eigen::MatrixXd path_theta    = Eigen::MatrixXd::Zero(MAXT + 1, d); path_theta.row(0)    = THETA_INIT;
    Eigen::MatrixXd path_av_theta = Eigen::MatrixXd::Zero(MAXT + 1, d); path_av_theta.row(0) = THETA_INIT;
    Eigen::MatrixXd path_grad     = Eigen::MatrixXd::Zero(MAXT,     d);

    // Initialise generic pair-object: it will compute pairwise quantities along the optimisation
    pair_class pair;

    // Compute scaling constant
    double scale;
    switch(METHODFLAG){
    case 0:
        scale = 1/static_cast<double>(n) ;
        break;
    case 1:
        scale = 1/static_cast<double>(NU);
        break;
    case 2:
        scale = 1/static_cast<double>(NU);
        break;
    }

    //Rcpp::Rcout << "Method:" << METHODFLAG << "\n ";
    unsigned int k_counter = 0;
    Eigen::VectorXd theta_t = THETA_INIT;
    for(unsigned int t = 1; t <= MAXT; t++){

        // check user interruption
        Rcpp::checkUserInterrupt();
        Rcpp::Rcout << "\rIteration:" << t << " ";

        /////////////////////
        // SAMPLING STEP   //
        /////////////////////
        Rcpp::NumericMatrix sampling_weights(n,kk);
        std::fill(sampling_weights.begin(), sampling_weights.end(), 0) ;

        double prob;
        switch(METHODFLAG){
        case 0:
            std::fill( sampling_weights.begin(), sampling_weights.end(), 1);
            break;
        case 1:
            prob = 1/static_cast<double>(n);
            sampling_weights = rmultinom_wrapper(prob, n, NU, kk);
            break;
        case 2:
            prob = static_cast<double>(NU)/static_cast<double>(n);
            for(unsigned int i = 0; i < n; i++){
                for(unsigned int k = 0; k < kk; k++){
                    if(R::runif(0,1) < prob ) sampling_weights(i, k) = 1;
                }
            }
            break;
        }


        //////////////////
        /*   GRADIENT   */
        //////////////////
        Eigen::VectorXd ngradient_t = Eigen::VectorXd::Zero(d);

        //if(theta_t(0)>2) theta_t(0) = 2;
        // Read parameter invariant to next loops
        double lxi  = theta_t(0);
        double artanhrho = theta_t(1);
        Eigen::VectorXd beta = theta_t.segment(2, r);

        for(unsigned int i = 0; i < n; i++){

            // Pair counter (used to index weights)
            unsigned int k_counter = 0;

            // Read unit i covariates
            Eigen::VectorXd x_i = X.row(i);

            for(unsigned int j = 1; j < p; j++){

                // Read quantities dependent on j
                unsigned int n_j = DATA(i, j);
                double alpha_j   = theta_t(r + 2 + j);

                for( unsigned int jp = 0; jp < j; jp++){

                    unsigned int weight = sampling_weights(i, k_counter);

                    if(weight != 0){
                        // Read quantities dependent on j'
                        unsigned int n_jp = DATA(i, jp);
                        double alpha_jp   = theta_t(r + 2 + jp);

                        pair.setup(j, jp, n_j, n_jp, x_i, beta, alpha_j, alpha_jp, lxi, artanhrho, p);
                        pair.compute_intermediate_( );
                        pair.compute_dintermediate_();

                        ngradient_t -= weight * pair.compute_gradient_();
                    }


                    k_counter ++;
                }
            }
        }

        ngradient_t *= scale;
        ///////////////////////////
        /*    PARAMETERS UPDATE  */
        ///////////////////////////

        if(STEPSIZEFLAG){
            Eigen::VectorXd stepsize_t = stepsizeV * pow(t, -.501);
            theta_t -= Eigen::VectorXd(stepsize_t.array() * ngradient_t.array());
        }else{
            double stepsize_t = STEPSIZE * pow(t, -.501);
            theta_t -= stepsize_t * ngradient_t;
        }

        ///////////////////////////////
        /* STORE ITERATION QUANTITIES  */
        /////////////////////////////////
        path_theta.row(t ) = theta_t;
        path_grad.row(t-1) = ngradient_t;

        // averaging after burnsize
        if(t <= BURN){
            path_av_theta.row(t) = path_theta.row(t);
        }else{
            path_av_theta.row(t) = ( (t - BURN - 1) * path_av_theta.row(t - 1) + path_theta.row(t) ) / (t - BURN);
        }


    }

    // clock.tock("main");
    // clock.stop("clock");

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("path_theta") = path_theta,
        Rcpp::Named("path_av_theta") = path_av_theta,
        Rcpp::Named("path_grad") = path_grad,
        Rcpp::Named("scale") = scale,
        Rcpp::Named("n") = n,
        Rcpp::Named("methodflag") = METHODFLAG
    );

    return output;
}





