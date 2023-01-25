// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// collapsed_efd_cpp
Rcpp::List collapsed_efd_cpp(arma::mat data, arma::colvec alpha, arma::colvec beta, arma::colvec tau, arma::colvec p, int K, int thin, int niter, double warmup, arma::colvec keep_index, int verbose);
RcppExport SEXP _EFLDA_collapsed_efd_cpp(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP pSEXP, SEXP KSEXP, SEXP thinSEXP, SEXP niterSEXP, SEXP warmupSEXP, SEXP keep_indexSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_efd_cpp(data, alpha, beta, tau, p, K, thin, niter, warmup, keep_index, verbose));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_lda_cpp
Rcpp::List collapsed_lda_cpp(arma::mat data, arma::colvec alpha, arma::colvec beta, int K, int thin, int niter, double warmup, arma::colvec keep_index, int verbose);
RcppExport SEXP _EFLDA_collapsed_lda_cpp(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP thinSEXP, SEXP niterSEXP, SEXP warmupSEXP, SEXP keep_indexSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_lda_cpp(data, alpha, beta, K, thin, niter, warmup, keep_index, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EFLDA_collapsed_efd_cpp", (DL_FUNC) &_EFLDA_collapsed_efd_cpp, 11},
    {"_EFLDA_collapsed_lda_cpp", (DL_FUNC) &_EFLDA_collapsed_lda_cpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_EFLDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
