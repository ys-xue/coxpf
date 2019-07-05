// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lad_shape
arma::vec lad_shape(const arma::mat& profileMat, const float bandwidth, const int gridprecision, const arma::vec& timevec);
RcppExport SEXP _coxpf_lad_shape(SEXP profileMatSEXP, SEXP bandwidthSEXP, SEXP gridprecisionSEXP, SEXP timevecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type profileMat(profileMatSEXP);
    Rcpp::traits::input_parameter< const float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< const int >::type gridprecision(gridprecisionSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type timevec(timevecSEXP);
    rcpp_result_gen = Rcpp::wrap(lad_shape(profileMat, bandwidth, gridprecision, timevec));
    return rcpp_result_gen;
END_RCPP
}
// lad_variability
arma::vec lad_variability(const arma::mat& cprofileMat, const float bandwidth, const int gridprecision, const arma::vec& timevec);
RcppExport SEXP _coxpf_lad_variability(SEXP cprofileMatSEXP, SEXP bandwidthSEXP, SEXP gridprecisionSEXP, SEXP timevecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cprofileMat(cprofileMatSEXP);
    Rcpp::traits::input_parameter< const float >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< const int >::type gridprecision(gridprecisionSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type timevec(timevecSEXP);
    rcpp_result_gen = Rcpp::wrap(lad_variability(cprofileMat, bandwidth, gridprecision, timevec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_coxpf_lad_shape", (DL_FUNC) &_coxpf_lad_shape, 4},
    {"_coxpf_lad_variability", (DL_FUNC) &_coxpf_lad_variability, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_coxpf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
