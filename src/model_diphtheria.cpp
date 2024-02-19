// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <boost/numeric/odeint.hpp>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

//' @title Run an SEIHR ODE model for diphtheria using a Boost solver
//'
//' @description A compartmental model for diphtheria with parameters to help
//' account for case reporting rate, delays in seeking hospitalisation, and the
//' time spent in the hospitalised compartment.
//'
//' This function is intended to only be called internally from
//' [model_diphtheria_cpp()].
//'
//' @param initial_state A matrix for the initial state of the compartments.
//' @param transmissibility The transmission rate \eqn{\beta}.
//' @param infectiousness_rate The rate of transition from exposed to infectious
//' \eqn{\alpha}.
//' @param recovery_rate The recovery rate \eqn{\gamma}.
//' @param reporting_rate The recovery rate \eqn{\r}.
//' @param prop_hosp The proportion of individuals hospitalised \eqn{\eta}.
//' @param hosp_entry_rate The rate at which individuals are hospitalised,
//' represented as 1 / time to hospitalisation \eqn{\tau_1}.
//' @param hosp_exit_rate The rate at which individuals are discharged from
//' hospital, represented as 1 / time to discharge \eqn{\tau_2}.
//' @param rate_interventions A named list of `<rate_intervention>` objects.
//' @param time_dependence A named list of functions for parameter time
//' dependence.
//' @param pop_change_times A numeric vector of times and which the population
//' of susceptibles changes.
//' @param pop_change_values An Rcpp List of numeric vectors giving the value of
//' changes to each demographic group at each change in population.
//' @param time_end The end time of the simulation.
//' @param increment The time increment of the simulation.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix.
//' The second list element is a vector of timesteps.
//' @keywords internal
// [[Rcpp::export(name=".model_diphtheria_cpp")]]
Rcpp::List model_diphtheria_internal(
    const Eigen::MatrixXd &initial_state, const double &transmissibility,
    const double &infectiousness_rate, const double &recovery_rate,
    const double &reporting_rate, const double &prop_hosp,
    const double &hosp_entry_rate, const double &hosp_exit_rate,
    const Rcpp::List &rate_interventions, const Rcpp::List &time_dependence,
    const Rcpp::NumericVector &pop_change_times,
    const Rcpp::List &pop_change_values,
    const double &time_end = 100.0,  // double required by boost solver
    const double &increment = 1.0) {
  // initial conditions from input
  odetools::state_type x = initial_state;

  // create a map of the model parameters
  std::unordered_map<std::string, double> model_params{
      {"transmissibility", transmissibility},
      {"infectiousness_rate", infectiousness_rate},
      {"recovery_rate", recovery_rate},
      {"reporting_rate", reporting_rate},
      {"prop_hosp", prop_hosp},
      {"hosp_entry_rate", hosp_entry_rate},
      {"hosp_exit_rate", hosp_exit_rate}};

  // create a map of the rate interventions
  const std::unordered_map<std::string, intervention::rate_intervention>
      rate_interventions_cpp =
          intervention::rate_intervention_cpp(rate_interventions);

  // prepare population changes if any
  const population::population_change pop_change(
      pop_change_times, pop_change_values, initial_state.rows());

  // create a diphtheria epidemic with parameters
  epidemics::epidemic_diphtheria this_model(
      model_params, rate_interventions_cpp, time_dependence, pop_change);

  // prepare storage containers for the observer
  std::vector<odetools::state_type> x_vec;  // is a vector of MatrixXd
  std::vector<double> times;

  // a controlled stepper for constant step sizes
  boost::numeric::odeint::runge_kutta4<
      odetools::state_type, double, odetools::state_type, double,
      boost::numeric::odeint::vector_space_algebra>
      stepper;

  // run the function without assignment
  boost::numeric::odeint::integrate_const(stepper, this_model, x, 0.0, time_end,
                                          increment,
                                          odetools::observer(x_vec, times));

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x_vec),
                            Rcpp::Named("time") = Rcpp::wrap(times));
  // return process_data(data);
}
