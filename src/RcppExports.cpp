// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pair_wrapper
Rcpp::List pair_wrapper(unsigned int j, unsigned int jp, unsigned int n_j, unsigned int n_jp, unsigned int p, double alpha_j, double alpha_jp, Eigen::VectorXd x, Eigen::VectorXd beta, double lxi, double artanhrho, unsigned int ind, bool verboseS, bool verboseSind);
RcppExport SEXP _gammaFrailty_pair_wrapper(SEXP jSEXP, SEXP jpSEXP, SEXP n_jSEXP, SEXP n_jpSEXP, SEXP pSEXP, SEXP alpha_jSEXP, SEXP alpha_jpSEXP, SEXP xSEXP, SEXP betaSEXP, SEXP lxiSEXP, SEXP artanhrhoSEXP, SEXP indSEXP, SEXP verboseSSEXP, SEXP verboseSindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type jp(jpSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_j(n_jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_jp(n_jpSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_j(alpha_jSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_jp(alpha_jpSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lxi(lxiSEXP);
    Rcpp::traits::input_parameter< double >::type artanhrho(artanhrhoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ind(indSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseS(verboseSSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseSind(verboseSindSEXP);
    rcpp_result_gen = Rcpp::wrap(pair_wrapper(j, jp, n_j, n_jp, p, alpha_j, alpha_jp, x, beta, lxi, artanhrho, ind, verboseS, verboseSind));
    return rcpp_result_gen;
END_RCPP
}
// ncl
Rcpp::List ncl(Eigen::VectorXd theta, Eigen::MatrixXd data, Eigen::MatrixXd X, bool printFLAG, const int PAIRS_RANGE);
RcppExport SEXP _gammaFrailty_ncl(SEXP thetaSEXP, SEXP dataSEXP, SEXP XSEXP, SEXP printFLAGSEXP, SEXP PAIRS_RANGESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type printFLAG(printFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type PAIRS_RANGE(PAIRS_RANGESEXP);
    rcpp_result_gen = Rcpp::wrap(ncl(theta, data, X, printFLAG, PAIRS_RANGE));
    return rcpp_result_gen;
END_RCPP
}
// gammaFrailty
Rcpp::List gammaFrailty(Eigen::VectorXd THETA_INIT, Eigen::MatrixXd DATA, const Eigen::MatrixXd X, const unsigned int MAXT, const unsigned int BURN, const double STEPSIZE, const double STEPSIZE0, const double NU, const int METHODFLAG, const bool VERBOSEFLAG, const bool STEPSIZEFLAG, const double par1, const double par2, const double par3, int PAIRS_RANGE);
RcppExport SEXP _gammaFrailty_gammaFrailty(SEXP THETA_INITSEXP, SEXP DATASEXP, SEXP XSEXP, SEXP MAXTSEXP, SEXP BURNSEXP, SEXP STEPSIZESEXP, SEXP STEPSIZE0SEXP, SEXP NUSEXP, SEXP METHODFLAGSEXP, SEXP VERBOSEFLAGSEXP, SEXP STEPSIZEFLAGSEXP, SEXP par1SEXP, SEXP par2SEXP, SEXP par3SEXP, SEXP PAIRS_RANGESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA_INIT(THETA_INITSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type DATA(DATASEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type BURN(BURNSEXP);
    Rcpp::traits::input_parameter< const double >::type STEPSIZE(STEPSIZESEXP);
    Rcpp::traits::input_parameter< const double >::type STEPSIZE0(STEPSIZE0SEXP);
    Rcpp::traits::input_parameter< const double >::type NU(NUSEXP);
    Rcpp::traits::input_parameter< const int >::type METHODFLAG(METHODFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSEFLAG(VERBOSEFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type STEPSIZEFLAG(STEPSIZEFLAGSEXP);
    Rcpp::traits::input_parameter< const double >::type par1(par1SEXP);
    Rcpp::traits::input_parameter< const double >::type par2(par2SEXP);
    Rcpp::traits::input_parameter< const double >::type par3(par3SEXP);
    Rcpp::traits::input_parameter< int >::type PAIRS_RANGE(PAIRS_RANGESEXP);
    rcpp_result_gen = Rcpp::wrap(gammaFrailty(THETA_INIT, DATA, X, MAXT, BURN, STEPSIZE, STEPSIZE0, NU, METHODFLAG, VERBOSEFLAG, STEPSIZEFLAG, par1, par2, par3, PAIRS_RANGE));
    return rcpp_result_gen;
END_RCPP
}
// zofr_cpp
double zofr_cpp(const double r);
RcppExport SEXP _gammaFrailty_zofr_cpp(SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(zofr_cpp(r));
    return rcpp_result_gen;
END_RCPP
}
// rofz_cpp
double rofz_cpp(const double z);
RcppExport SEXP _gammaFrailty_rofz_cpp(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(rofz_cpp(z));
    return rcpp_result_gen;
END_RCPP
}
// drofz_cpp
double drofz_cpp(const double z);
RcppExport SEXP _gammaFrailty_drofz_cpp(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(drofz_cpp(z));
    return rcpp_result_gen;
END_RCPP
}
// rmultinom_wrapper
Rcpp::NumericMatrix rmultinom_wrapper(const double prob, const unsigned int classes, const unsigned int batch, const unsigned int K);
RcppExport SEXP _gammaFrailty_rmultinom_wrapper(SEXP probSEXP, SEXP classesSEXP, SEXP batchSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinom_wrapper(prob, classes, batch, K));
    return rcpp_result_gen;
END_RCPP
}
// sampleJ
Eigen::MatrixXd sampleJ(Eigen::VectorXd THETA, Eigen::MatrixXd DATA, Eigen::MatrixXd X, const bool PRINTFLAG);
RcppExport SEXP _gammaFrailty_sampleJ(SEXP THETASEXP, SEXP DATASEXP, SEXP XSEXP, SEXP PRINTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type DATA(DATASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const bool >::type PRINTFLAG(PRINTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleJ(THETA, DATA, X, PRINTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// sampleH
Eigen::MatrixXd sampleH(Eigen::VectorXd THETA, Eigen::MatrixXd DATA, Eigen::MatrixXd X, const bool PRINTFLAG, const bool INVERTFLAG);
RcppExport SEXP _gammaFrailty_sampleH(SEXP THETASEXP, SEXP DATASEXP, SEXP XSEXP, SEXP PRINTFLAGSEXP, SEXP INVERTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type DATA(DATASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const bool >::type PRINTFLAG(PRINTFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type INVERTFLAG(INVERTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleH(THETA, DATA, X, PRINTFLAG, INVERTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// sampleVar
Rcpp::List sampleVar(Eigen::VectorXd THETA, Eigen::MatrixXd DATA, Eigen::MatrixXd X, const unsigned int NU, const unsigned int METHOD, const unsigned int RANGE, const bool TOTFLAG, const bool PRINTFLAG);
RcppExport SEXP _gammaFrailty_sampleVar(SEXP THETASEXP, SEXP DATASEXP, SEXP XSEXP, SEXP NUSEXP, SEXP METHODSEXP, SEXP RANGESEXP, SEXP TOTFLAGSEXP, SEXP PRINTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type DATA(DATASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NU(NUSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type METHOD(METHODSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type RANGE(RANGESEXP);
    Rcpp::traits::input_parameter< const bool >::type TOTFLAG(TOTFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type PRINTFLAG(PRINTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleVar(THETA, DATA, X, NU, METHOD, RANGE, TOTFLAG, PRINTFLAG));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gammaFrailty_pair_wrapper", (DL_FUNC) &_gammaFrailty_pair_wrapper, 14},
    {"_gammaFrailty_ncl", (DL_FUNC) &_gammaFrailty_ncl, 5},
    {"_gammaFrailty_gammaFrailty", (DL_FUNC) &_gammaFrailty_gammaFrailty, 15},
    {"_gammaFrailty_zofr_cpp", (DL_FUNC) &_gammaFrailty_zofr_cpp, 1},
    {"_gammaFrailty_rofz_cpp", (DL_FUNC) &_gammaFrailty_rofz_cpp, 1},
    {"_gammaFrailty_drofz_cpp", (DL_FUNC) &_gammaFrailty_drofz_cpp, 1},
    {"_gammaFrailty_rmultinom_wrapper", (DL_FUNC) &_gammaFrailty_rmultinom_wrapper, 4},
    {"_gammaFrailty_sampleJ", (DL_FUNC) &_gammaFrailty_sampleJ, 4},
    {"_gammaFrailty_sampleH", (DL_FUNC) &_gammaFrailty_sampleH, 5},
    {"_gammaFrailty_sampleVar", (DL_FUNC) &_gammaFrailty_sampleVar, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_gammaFrailty(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
