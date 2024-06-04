// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// whichOne
int whichOne(IntegerVector x);
RcppExport SEXP _EFLDA_whichOne(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(whichOne(x));
    return rcpp_result_gen;
END_RCPP
}
// whichIndex
int whichIndex(NumericVector x, double val);
RcppExport SEXP _EFLDA_whichIndex(SEXP xSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(whichIndex(x, val));
    return rcpp_result_gen;
END_RCPP
}
// whichMax
int whichMax(arma::colvec vec);
RcppExport SEXP _EFLDA_whichMax(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(whichMax(vec));
    return rcpp_result_gen;
END_RCPP
}
// oneSampleMultinom
IntegerVector oneSampleMultinom(NumericVector probs);
RcppExport SEXP _EFLDA_oneSampleMultinom(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(oneSampleMultinom(probs));
    return rcpp_result_gen;
END_RCPP
}
// whichMultinom
int whichMultinom(NumericVector probs, int d, int w);
RcppExport SEXP _EFLDA_whichMultinom(SEXP probsSEXP, SEXP dSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(whichMultinom(probs, d, w));
    return rcpp_result_gen;
END_RCPP
}
// isInVector
bool isInVector(int value, NumericVector vec);
RcppExport SEXP _EFLDA_isInVector(SEXP valueSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(isInVector(value, vec));
    return rcpp_result_gen;
END_RCPP
}
// fDir
double fDir(arma::colvec x, arma::colvec alpha, arma::colvec tau, int position);
RcppExport SEXP _EFLDA_fDir(SEXP xSEXP, SEXP alphaSEXP, SEXP tauSEXP, SEXP positionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type position(positionSEXP);
    rcpp_result_gen = Rcpp::wrap(fDir(x, alpha, tau, position));
    return rcpp_result_gen;
END_RCPP
}
// cluster_allocation_cpp
arma::mat cluster_allocation_cpp(arma::cube theta_post, int D, int n_post, int K, arma::colvec alpha, arma::colvec tau, arma::colvec p);
RcppExport SEXP _EFLDA_cluster_allocation_cpp(SEXP theta_postSEXP, SEXP DSEXP, SEXP n_postSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP tauSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type theta_post(theta_postSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type n_post(n_postSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster_allocation_cpp(theta_post, D, n_post, K, alpha, tau, p));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_lda_cpp
Rcpp::List collapsed_lda_cpp(NumericMatrix& data, arma::colvec& alpha, arma::colvec& beta, int K, int niter, NumericVector& keep_index, Rcpp::List& z_init, int verbose, int nupd);
RcppExport SEXP _EFLDA_collapsed_lda_cpp(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP niterSEXP, SEXP keep_indexSEXP, SEXP z_initSEXP, SEXP verboseSEXP, SEXP nupdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type nupd(nupdSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_lda_cpp(data, alpha, beta, K, niter, keep_index, z_init, verbose, nupd));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_efd_cpp
Rcpp::List collapsed_efd_cpp(NumericMatrix data, arma::colvec& alpha, arma::colvec& beta, arma::colvec& tau, arma::colvec& p, int K, int niter, NumericVector& keep_index, Rcpp::List z_init, int verbose, int nupd);
RcppExport SEXP _EFLDA_collapsed_efd_cpp(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP pSEXP, SEXP KSEXP, SEXP niterSEXP, SEXP keep_indexSEXP, SEXP z_initSEXP, SEXP verboseSEXP, SEXP nupdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type nupd(nupdSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_efd_cpp(data, alpha, beta, tau, p, K, niter, keep_index, z_init, verbose, nupd));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_lda_cpp_pred
Rcpp::List collapsed_lda_cpp_pred(NumericMatrix& data, arma::colvec& alpha, arma::colvec& beta, arma::mat phi_post_mean, int K, int niter, NumericVector& keep_index, Rcpp::List& z_init, int verbose, int nupd);
RcppExport SEXP _EFLDA_collapsed_lda_cpp_pred(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP phi_post_meanSEXP, SEXP KSEXP, SEXP niterSEXP, SEXP keep_indexSEXP, SEXP z_initSEXP, SEXP verboseSEXP, SEXP nupdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_post_mean(phi_post_meanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type nupd(nupdSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_lda_cpp_pred(data, alpha, beta, phi_post_mean, K, niter, keep_index, z_init, verbose, nupd));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_efd_cpp_pred
Rcpp::List collapsed_efd_cpp_pred(NumericMatrix data, arma::colvec& alpha, arma::colvec& beta, arma::colvec& tau, arma::colvec& p, arma::mat phi_post_mean, int K, int niter, NumericVector& keep_index, Rcpp::List z_init, int verbose, int nupd);
RcppExport SEXP _EFLDA_collapsed_efd_cpp_pred(SEXP dataSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP pSEXP, SEXP phi_post_meanSEXP, SEXP KSEXP, SEXP niterSEXP, SEXP keep_indexSEXP, SEXP z_initSEXP, SEXP verboseSEXP, SEXP nupdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_post_mean(phi_post_meanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type keep_index(keep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type z_init(z_initSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type nupd(nupdSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_efd_cpp_pred(data, alpha, beta, tau, p, phi_post_mean, K, niter, keep_index, z_init, verbose, nupd));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EFLDA_whichOne", (DL_FUNC) &_EFLDA_whichOne, 1},
    {"_EFLDA_whichIndex", (DL_FUNC) &_EFLDA_whichIndex, 2},
    {"_EFLDA_whichMax", (DL_FUNC) &_EFLDA_whichMax, 1},
    {"_EFLDA_oneSampleMultinom", (DL_FUNC) &_EFLDA_oneSampleMultinom, 1},
    {"_EFLDA_whichMultinom", (DL_FUNC) &_EFLDA_whichMultinom, 3},
    {"_EFLDA_isInVector", (DL_FUNC) &_EFLDA_isInVector, 2},
    {"_EFLDA_fDir", (DL_FUNC) &_EFLDA_fDir, 4},
    {"_EFLDA_cluster_allocation_cpp", (DL_FUNC) &_EFLDA_cluster_allocation_cpp, 7},
    {"_EFLDA_collapsed_lda_cpp", (DL_FUNC) &_EFLDA_collapsed_lda_cpp, 9},
    {"_EFLDA_collapsed_efd_cpp", (DL_FUNC) &_EFLDA_collapsed_efd_cpp, 11},
    {"_EFLDA_collapsed_lda_cpp_pred", (DL_FUNC) &_EFLDA_collapsed_lda_cpp_pred, 10},
    {"_EFLDA_collapsed_efd_cpp_pred", (DL_FUNC) &_EFLDA_collapsed_efd_cpp_pred, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_EFLDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
