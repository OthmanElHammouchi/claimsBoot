// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// glmBoot
Rcpp::NumericVector glmBoot(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::String boot_type, Rcpp::String opt, int seed);
RcppExport SEXP _claimsBoot_glmBoot(SEXP triangleSEXP, SEXP n_bootSEXP, SEXP boot_typeSEXP, SEXP optSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type triangle(triangleSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type boot_type(boot_typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type opt(optSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(glmBoot(triangle, n_boot, boot_type, opt, seed));
    return rcpp_result_gen;
END_RCPP
}
// glmSim
Rcpp::DataFrame glmSim(Rcpp::NumericMatrix triangle, Rcpp::String sim_type, int n_boot, Rcpp::NumericVector factors, Rcpp::CharacterVector boot_types, Rcpp::String sim_dist, bool show_progress, int seed);
RcppExport SEXP _claimsBoot_glmSim(SEXP triangleSEXP, SEXP sim_typeSEXP, SEXP n_bootSEXP, SEXP factorsSEXP, SEXP boot_typesSEXP, SEXP sim_distSEXP, SEXP show_progressSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type triangle(triangleSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sim_type(sim_typeSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type boot_types(boot_typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sim_dist(sim_distSEXP);
    Rcpp::traits::input_parameter< bool >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(glmSim(triangle, sim_type, n_boot, factors, boot_types, sim_dist, show_progress, seed));
    return rcpp_result_gen;
END_RCPP
}
// mackBoot
Rcpp::NumericVector mackBoot(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::String boot_type, Rcpp::String opt, bool cond, int seed);
RcppExport SEXP _claimsBoot_mackBoot(SEXP triangleSEXP, SEXP n_bootSEXP, SEXP boot_typeSEXP, SEXP optSEXP, SEXP condSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type triangle(triangleSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type boot_type(boot_typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type opt(optSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(mackBoot(triangle, n_boot, boot_type, opt, cond, seed));
    return rcpp_result_gen;
END_RCPP
}
// mackSim
Rcpp::DataFrame mackSim(Rcpp::NumericMatrix triangle, Rcpp::String sim_type, int n_boot, Rcpp::NumericVector mean_factors, Rcpp::NumericVector sd_factors, Rcpp::CharacterVector boot_types, Rcpp::String sim_dist, bool show_progress, int seed);
RcppExport SEXP _claimsBoot_mackSim(SEXP triangleSEXP, SEXP sim_typeSEXP, SEXP n_bootSEXP, SEXP mean_factorsSEXP, SEXP sd_factorsSEXP, SEXP boot_typesSEXP, SEXP sim_distSEXP, SEXP show_progressSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type triangle(triangleSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sim_type(sim_typeSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mean_factors(mean_factorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sd_factors(sd_factorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type boot_types(boot_typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sim_dist(sim_distSEXP);
    Rcpp::traits::input_parameter< bool >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(mackSim(triangle, sim_type, n_boot, mean_factors, sd_factors, boot_types, sim_dist, show_progress, seed));
    return rcpp_result_gen;
END_RCPP
}
// test_pois
Rcpp::NumericVector test_pois(int n, double lambda);
RcppExport SEXP _claimsBoot_test_pois(SEXP nSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(test_pois(n, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_claimsBoot_glmBoot", (DL_FUNC) &_claimsBoot_glmBoot, 5},
    {"_claimsBoot_glmSim", (DL_FUNC) &_claimsBoot_glmSim, 8},
    {"_claimsBoot_mackBoot", (DL_FUNC) &_claimsBoot_mackBoot, 6},
    {"_claimsBoot_mackSim", (DL_FUNC) &_claimsBoot_mackSim, 9},
    {"_claimsBoot_test_pois", (DL_FUNC) &_claimsBoot_test_pois, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_claimsBoot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
