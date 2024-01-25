// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/epidemics.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// model_default_cpp_internal
Rcpp::List model_default_cpp_internal(const Eigen::MatrixXd& initial_state, const double& transmissibility, const double& infectiousness_rate, const double& recovery_rate, const Eigen::MatrixXd& contact_matrix, const Rcpp::NumericVector& npi_time_begin, const Rcpp::NumericVector& npi_time_end, const Rcpp::NumericMatrix& npi_cr, const Eigen::MatrixXd& vax_time_begin, const Eigen::MatrixXd& vax_time_end, const Eigen::MatrixXd& vax_nu, const Rcpp::List& rate_interventions, const Rcpp::List& time_dependence, const double& time_end, const double& increment);
RcppExport SEXP _epidemics_model_default_cpp_internal(SEXP initial_stateSEXP, SEXP transmissibilitySEXP, SEXP infectiousness_rateSEXP, SEXP recovery_rateSEXP, SEXP contact_matrixSEXP, SEXP npi_time_beginSEXP, SEXP npi_time_endSEXP, SEXP npi_crSEXP, SEXP vax_time_beginSEXP, SEXP vax_time_endSEXP, SEXP vax_nuSEXP, SEXP rate_interventionsSEXP, SEXP time_dependenceSEXP, SEXP time_endSEXP, SEXP incrementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type initial_state(initial_stateSEXP);
    Rcpp::traits::input_parameter< const double& >::type transmissibility(transmissibilitySEXP);
    Rcpp::traits::input_parameter< const double& >::type infectiousness_rate(infectiousness_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type recovery_rate(recovery_rateSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type contact_matrix(contact_matrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type npi_time_begin(npi_time_beginSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type npi_time_end(npi_time_endSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type npi_cr(npi_crSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_time_begin(vax_time_beginSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_time_end(vax_time_endSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_nu(vax_nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type rate_interventions(rate_interventionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type time_dependence(time_dependenceSEXP);
    Rcpp::traits::input_parameter< const double& >::type time_end(time_endSEXP);
    Rcpp::traits::input_parameter< const double& >::type increment(incrementSEXP);
    rcpp_result_gen = Rcpp::wrap(model_default_cpp_internal(initial_state, transmissibility, infectiousness_rate, recovery_rate, contact_matrix, npi_time_begin, npi_time_end, npi_cr, vax_time_begin, vax_time_end, vax_nu, rate_interventions, time_dependence, time_end, increment));
    return rcpp_result_gen;
END_RCPP
}
// model_diphtheria_cpp_internal
Rcpp::List model_diphtheria_cpp_internal(const Eigen::MatrixXd& initial_state, const double& transmissibility, const double& infectiousness_rate, const double& recovery_rate, const double& reporting_rate, const double& prop_hosp, const double& hosp_entry_rate, const double& hosp_exit_rate, const Rcpp::List& rate_interventions, const Rcpp::List& time_dependence, const Rcpp::NumericVector& pop_change_times, const Rcpp::List& pop_change_values, const double& time_end, const double& increment);
RcppExport SEXP _epidemics_model_diphtheria_cpp_internal(SEXP initial_stateSEXP, SEXP transmissibilitySEXP, SEXP infectiousness_rateSEXP, SEXP recovery_rateSEXP, SEXP reporting_rateSEXP, SEXP prop_hospSEXP, SEXP hosp_entry_rateSEXP, SEXP hosp_exit_rateSEXP, SEXP rate_interventionsSEXP, SEXP time_dependenceSEXP, SEXP pop_change_timesSEXP, SEXP pop_change_valuesSEXP, SEXP time_endSEXP, SEXP incrementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type initial_state(initial_stateSEXP);
    Rcpp::traits::input_parameter< const double& >::type transmissibility(transmissibilitySEXP);
    Rcpp::traits::input_parameter< const double& >::type infectiousness_rate(infectiousness_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type recovery_rate(recovery_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type reporting_rate(reporting_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type prop_hosp(prop_hospSEXP);
    Rcpp::traits::input_parameter< const double& >::type hosp_entry_rate(hosp_entry_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type hosp_exit_rate(hosp_exit_rateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type rate_interventions(rate_interventionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type time_dependence(time_dependenceSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type pop_change_times(pop_change_timesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pop_change_values(pop_change_valuesSEXP);
    Rcpp::traits::input_parameter< const double& >::type time_end(time_endSEXP);
    Rcpp::traits::input_parameter< const double& >::type increment(incrementSEXP);
    rcpp_result_gen = Rcpp::wrap(model_diphtheria_cpp_internal(initial_state, transmissibility, infectiousness_rate, recovery_rate, reporting_rate, prop_hosp, hosp_entry_rate, hosp_exit_rate, rate_interventions, time_dependence, pop_change_times, pop_change_values, time_end, increment));
    return rcpp_result_gen;
END_RCPP
}
// model_vacamole_cpp_internal
Rcpp::List model_vacamole_cpp_internal(const Eigen::MatrixXd& initial_state, const double& transmissibility, const double& transmissibility_vax, const double& infectiousness_rate, const double& mortality_rate, const double& mortality_rate_vax, const double& hospitalisation_rate, const double& hospitalisation_rate_vax, const double& recovery_rate, const Eigen::MatrixXd& contact_matrix, const Rcpp::NumericVector& npi_time_begin, const Rcpp::NumericVector& npi_time_end, const Rcpp::NumericMatrix& npi_cr, const Eigen::MatrixXd& vax_time_begin, const Eigen::MatrixXd& vax_time_end, const Eigen::MatrixXd& vax_nu, const Rcpp::List& rate_interventions, const Rcpp::List& time_dependence, const double& time_end, const double& increment);
RcppExport SEXP _epidemics_model_vacamole_cpp_internal(SEXP initial_stateSEXP, SEXP transmissibilitySEXP, SEXP transmissibility_vaxSEXP, SEXP infectiousness_rateSEXP, SEXP mortality_rateSEXP, SEXP mortality_rate_vaxSEXP, SEXP hospitalisation_rateSEXP, SEXP hospitalisation_rate_vaxSEXP, SEXP recovery_rateSEXP, SEXP contact_matrixSEXP, SEXP npi_time_beginSEXP, SEXP npi_time_endSEXP, SEXP npi_crSEXP, SEXP vax_time_beginSEXP, SEXP vax_time_endSEXP, SEXP vax_nuSEXP, SEXP rate_interventionsSEXP, SEXP time_dependenceSEXP, SEXP time_endSEXP, SEXP incrementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type initial_state(initial_stateSEXP);
    Rcpp::traits::input_parameter< const double& >::type transmissibility(transmissibilitySEXP);
    Rcpp::traits::input_parameter< const double& >::type transmissibility_vax(transmissibility_vaxSEXP);
    Rcpp::traits::input_parameter< const double& >::type infectiousness_rate(infectiousness_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mortality_rate(mortality_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mortality_rate_vax(mortality_rate_vaxSEXP);
    Rcpp::traits::input_parameter< const double& >::type hospitalisation_rate(hospitalisation_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type hospitalisation_rate_vax(hospitalisation_rate_vaxSEXP);
    Rcpp::traits::input_parameter< const double& >::type recovery_rate(recovery_rateSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type contact_matrix(contact_matrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type npi_time_begin(npi_time_beginSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type npi_time_end(npi_time_endSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type npi_cr(npi_crSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_time_begin(vax_time_beginSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_time_end(vax_time_endSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type vax_nu(vax_nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type rate_interventions(rate_interventionsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type time_dependence(time_dependenceSEXP);
    Rcpp::traits::input_parameter< const double& >::type time_end(time_endSEXP);
    Rcpp::traits::input_parameter< const double& >::type increment(incrementSEXP);
    rcpp_result_gen = Rcpp::wrap(model_vacamole_cpp_internal(initial_state, transmissibility, transmissibility_vax, infectiousness_rate, mortality_rate, mortality_rate_vax, hospitalisation_rate, hospitalisation_rate_vax, recovery_rate, contact_matrix, npi_time_begin, npi_time_end, npi_cr, vax_time_begin, vax_time_end, vax_nu, rate_interventions, time_dependence, time_end, increment));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epidemics_model_default_cpp_internal", (DL_FUNC) &_epidemics_model_default_cpp_internal, 15},
    {"_epidemics_model_diphtheria_cpp_internal", (DL_FUNC) &_epidemics_model_diphtheria_cpp_internal, 14},
    {"_epidemics_model_vacamole_cpp_internal", (DL_FUNC) &_epidemics_model_vacamole_cpp_internal, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_epidemics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
