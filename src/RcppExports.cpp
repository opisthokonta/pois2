// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dpois2
NumericVector dpois2(IntegerVector& x1, IntegerVector& x2, NumericVector& lambda1, NumericVector& lambda2, NumericVector& lambda3);
RcppExport SEXP _pois2_dpois2(SEXP x1SEXP, SEXP x2SEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP lambda3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda3(lambda3SEXP);
    rcpp_result_gen = Rcpp::wrap(dpois2(x1, x2, lambda1, lambda2, lambda3));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pois2_dpois2", (DL_FUNC) &_pois2_dpois2, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_pois2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}