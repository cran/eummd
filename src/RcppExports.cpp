// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mmd_lap_Rcpp
Rcpp::List mmd_lap_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::NumericVector beta_);
RcppExport SEXP _eummd_mmd_lap_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP beta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    rcpp_result_gen = Rcpp::wrap(mmd_lap_Rcpp(X, Y, nX_, dX_, nY_, dY_, beta_));
    return rcpp_result_gen;
END_RCPP
}
// mmd_gau_Rcpp
Rcpp::List mmd_gau_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::NumericVector beta_);
RcppExport SEXP _eummd_mmd_gau_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP beta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    rcpp_result_gen = Rcpp::wrap(mmd_gau_Rcpp(X, Y, nX_, dX_, nY_, dY_, beta_));
    return rcpp_result_gen;
END_RCPP
}
// fast_median_diff_Rcpp
Rcpp::NumericVector fast_median_diff_Rcpp(Rcpp::NumericVector X_);
RcppExport SEXP _eummd_fast_median_diff_Rcpp(SEXP X_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X_(X_SEXP);
    rcpp_result_gen = Rcpp::wrap(fast_median_diff_Rcpp(X_));
    return rcpp_result_gen;
END_RCPP
}
// naive_median_diff_Rcpp
Rcpp::NumericVector naive_median_diff_Rcpp(Rcpp::NumericVector Z_, Rcpp::IntegerVector nZ_, Rcpp::IntegerVector dZ_, Rcpp::IntegerVector kmethod_);
RcppExport SEXP _eummd_naive_median_diff_Rcpp(SEXP Z_SEXP, SEXP nZ_SEXP, SEXP dZ_SEXP, SEXP kmethod_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Z_(Z_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nZ_(nZ_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dZ_(dZ_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type kmethod_(kmethod_SEXP);
    rcpp_result_gen = Rcpp::wrap(naive_median_diff_Rcpp(Z_, nZ_, dZ_, kmethod_));
    return rcpp_result_gen;
END_RCPP
}
// mmd_lap_pval_Rcpp
Rcpp::List mmd_lap_pval_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector seednum_, Rcpp::NumericVector beta_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_);
RcppExport SEXP _eummd_mmd_lap_pval_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numperm_SEXP, SEXP seednum_SEXP, SEXP beta_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    rcpp_result_gen = Rcpp::wrap(mmd_lap_pval_Rcpp(X, Y, nX_, dX_, nY_, dY_, numperm_, seednum_, beta_, twosided_, boundedminpval_));
    return rcpp_result_gen;
END_RCPP
}
// mmd_gau_pval_Rcpp
Rcpp::List mmd_gau_pval_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector seednum_, Rcpp::NumericVector beta_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_);
RcppExport SEXP _eummd_mmd_gau_pval_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numperm_SEXP, SEXP seednum_SEXP, SEXP beta_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    rcpp_result_gen = Rcpp::wrap(mmd_gau_pval_Rcpp(X, Y, nX_, dX_, nY_, dY_, numperm_, seednum_, beta_, twosided_, boundedminpval_));
    return rcpp_result_gen;
END_RCPP
}
// eummd_Rcpp
Rcpp::List eummd_Rcpp(Rcpp::NumericVector X_, Rcpp::NumericVector Y_, Rcpp::NumericVector beta_);
RcppExport SEXP _eummd_eummd_Rcpp(SEXP X_SEXP, SEXP Y_SEXP, SEXP beta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    rcpp_result_gen = Rcpp::wrap(eummd_Rcpp(X_, Y_, beta_));
    return rcpp_result_gen;
END_RCPP
}
// eummd_pval_Rcpp
Rcpp::List eummd_pval_Rcpp(Rcpp::NumericVector X_, Rcpp::NumericVector Y_, Rcpp::NumericVector beta_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector seednum_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_);
RcppExport SEXP _eummd_eummd_pval_Rcpp(SEXP X_SEXP, SEXP Y_SEXP, SEXP beta_SEXP, SEXP numperm_SEXP, SEXP seednum_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    rcpp_result_gen = Rcpp::wrap(eummd_pval_Rcpp(X_, Y_, beta_, numperm_, seednum_, twosided_, boundedminpval_));
    return rcpp_result_gen;
END_RCPP
}
// meammd_proj_Rcpp
double meammd_proj_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numproj_, Rcpp::IntegerVector seednum_, Rcpp::NumericVector beta_);
RcppExport SEXP _eummd_meammd_proj_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numproj_SEXP, SEXP seednum_SEXP, SEXP beta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numproj_(numproj_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    rcpp_result_gen = Rcpp::wrap(meammd_proj_Rcpp(X, Y, nX_, dX_, nY_, dY_, numproj_, seednum_, beta_));
    return rcpp_result_gen;
END_RCPP
}
// meammd_proj_pval_Rcpp
Rcpp::List meammd_proj_pval_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector numproj_, Rcpp::IntegerVector seednum_, Rcpp::NumericVector beta_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_, Rcpp::IntegerVector faster_);
RcppExport SEXP _eummd_meammd_proj_pval_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numperm_SEXP, SEXP numproj_SEXP, SEXP seednum_SEXP, SEXP beta_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP, SEXP faster_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numproj_(numproj_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type faster_(faster_SEXP);
    rcpp_result_gen = Rcpp::wrap(meammd_proj_pval_Rcpp(X, Y, nX_, dX_, nY_, dY_, numperm_, numproj_, seednum_, beta_, twosided_, boundedminpval_, faster_));
    return rcpp_result_gen;
END_RCPP
}
// meammd_dist_pval_Rcpp
double meammd_dist_pval_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector seednum_, Rcpp::NumericVector beta_, Rcpp::NumericVector pmethod_, Rcpp::NumericVector nmethod_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_);
RcppExport SEXP _eummd_meammd_dist_pval_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numperm_SEXP, SEXP seednum_SEXP, SEXP beta_SEXP, SEXP pmethod_SEXP, SEXP nmethod_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pmethod_(pmethod_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nmethod_(nmethod_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    rcpp_result_gen = Rcpp::wrap(meammd_dist_pval_Rcpp(X, Y, nX_, dX_, nY_, dY_, numperm_, seednum_, beta_, pmethod_, nmethod_, twosided_, boundedminpval_));
    return rcpp_result_gen;
END_RCPP
}
// energydist_Rcpp
Rcpp::List energydist_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_);
RcppExport SEXP _eummd_energydist_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    rcpp_result_gen = Rcpp::wrap(energydist_Rcpp(X, Y, nX_, dX_, nY_, dY_));
    return rcpp_result_gen;
END_RCPP
}
// energydist_pval_Rcpp
Rcpp::List energydist_pval_Rcpp(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::IntegerVector nX_, Rcpp::IntegerVector dX_, Rcpp::IntegerVector nY_, Rcpp::IntegerVector dY_, Rcpp::IntegerVector numperm_, Rcpp::IntegerVector seednum_, Rcpp::IntegerVector twosided_, Rcpp::IntegerVector boundedminpval_);
RcppExport SEXP _eummd_energydist_pval_Rcpp(SEXP XSEXP, SEXP YSEXP, SEXP nX_SEXP, SEXP dX_SEXP, SEXP nY_SEXP, SEXP dY_SEXP, SEXP numperm_SEXP, SEXP seednum_SEXP, SEXP twosided_SEXP, SEXP boundedminpval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nX_(nX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dX_(dX_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nY_(nY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type dY_(dY_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numperm_(numperm_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seednum_(seednum_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type twosided_(twosided_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boundedminpval_(boundedminpval_SEXP);
    rcpp_result_gen = Rcpp::wrap(energydist_pval_Rcpp(X, Y, nX_, dX_, nY_, dY_, numperm_, seednum_, twosided_, boundedminpval_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eummd_mmd_lap_Rcpp", (DL_FUNC) &_eummd_mmd_lap_Rcpp, 7},
    {"_eummd_mmd_gau_Rcpp", (DL_FUNC) &_eummd_mmd_gau_Rcpp, 7},
    {"_eummd_fast_median_diff_Rcpp", (DL_FUNC) &_eummd_fast_median_diff_Rcpp, 1},
    {"_eummd_naive_median_diff_Rcpp", (DL_FUNC) &_eummd_naive_median_diff_Rcpp, 4},
    {"_eummd_mmd_lap_pval_Rcpp", (DL_FUNC) &_eummd_mmd_lap_pval_Rcpp, 11},
    {"_eummd_mmd_gau_pval_Rcpp", (DL_FUNC) &_eummd_mmd_gau_pval_Rcpp, 11},
    {"_eummd_eummd_Rcpp", (DL_FUNC) &_eummd_eummd_Rcpp, 3},
    {"_eummd_eummd_pval_Rcpp", (DL_FUNC) &_eummd_eummd_pval_Rcpp, 7},
    {"_eummd_meammd_proj_Rcpp", (DL_FUNC) &_eummd_meammd_proj_Rcpp, 9},
    {"_eummd_meammd_proj_pval_Rcpp", (DL_FUNC) &_eummd_meammd_proj_pval_Rcpp, 13},
    {"_eummd_meammd_dist_pval_Rcpp", (DL_FUNC) &_eummd_meammd_dist_pval_Rcpp, 13},
    {"_eummd_energydist_Rcpp", (DL_FUNC) &_eummd_energydist_Rcpp, 6},
    {"_eummd_energydist_pval_Rcpp", (DL_FUNC) &_eummd_energydist_pval_Rcpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_eummd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
