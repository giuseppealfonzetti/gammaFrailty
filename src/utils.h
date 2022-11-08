#ifndef utils_H
#define utils_H

// Fisher transformation for correlation
//' @export
// [[Rcpp::export]]
double zofr_cpp(const double r){
    double z = .5*log((1+r)/(1-r));
    return z;
}

// Inverse of Fisher transformation
//' @export
// [[Rcpp::export]]
double rofz_cpp(const double z){
    double r = (exp(2*z)-1)/(exp(2*z)+1);
    return r;
}

// Derivative of inverse of Fisher transformation
//' @export
// [[Rcpp::export]]
double drofz_cpp(const double z){
    double r = rofz_cpp(z);
    double dr = 2*(1-r)*exp(2*z)/(exp(2*z)+1);
    return dr;
}

//'@export
// [[Rcpp::export]]
Rcpp::NumericMatrix rmultinom_wrapper(const double prob, const unsigned int classes, const unsigned int batch, const unsigned int K) {

    Rcpp::NumericVector probs(classes, prob);
    Rcpp::IntegerVector outcome(classes);
    R::rmultinom(batch, probs.begin(), classes, outcome.begin());


    Rcpp::NumericMatrix out(classes, K);
    for(unsigned int j = 0; j < K; j++){
        out(Rcpp::_,j) = outcome;
    }

    return out;
}

// //' @export
// // [[Rcpp::export]]
// double chooseC(double n, double k) {
//     return Rf_choose(n, k);
// }
#endif
