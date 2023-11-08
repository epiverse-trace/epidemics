// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <boost/numeric/odeint.hpp>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

//' @title Run an age-structured SEIR-V epidemic ODE model using a Boost solver
//'
//' @description A compartmental model with an optional non-pharmaceutical
//' intervention and an optional vaccination regime.
//'
//' This function is intended to only be called internally from
//' [model_default_cpp()].
//'
//' Allows heterogeneity in social contact patterns, and variable sizes of
//' demographic groups.
//' Also allows for group-specific initial proportions in each model
//' compartment, as well as group-specific vaccination start dates and
//' vaccination rates, and also group-specific effects of implementing a
//' non-pharmaceutical intervention.
//' The model only allows for single, population-wide rates of
//' transition between the 'susceptible' and 'exposed' compartments, between the
//' 'exposed' and 'infectious' compartments, and in the recovery rate.
//'
//' @param initial_state An `Eigen::MatrixXd` holding the initial state of
//' each demographic-compartmental combination. Rows must represent demographic
//' groups, while columns represent compartments in the order S, E, I, R, V.
//' @param contact_matrix An `Eigen::MatrixXd` holding the population
//' contact matrix.
//' @param beta The transmission rate \eqn{\beta}.
//' @param alpha The rate of transition from exposed to infectious \eqn{\alpha}.
//' @param gamma The recovery rate \eqn{\gamma}.
//' @param time_end The maximum time, defaults to 200.0.
//' @param intervention A non-pharmaceutical intervention applied during the
//' course of the epidemic, with a start and end time, and age-specific effect
//' on contacts. See [intervention()].
//' @param vaccination A vaccination regime followed during the
//' course of the epidemic, with a start and end time, and age-specific effect
//' on the transition of individuals from susceptible to vaccinated.
//' See [vaccination()].
//' @param increment The increment time, defaults to 0.1.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix (see [population()]).
//' The second list element is a vector of timesteps.
//' @keywords internal
// [[Rcpp::export(name=".model_default_cpp")]]
Rcpp::List model_default_cpp_internal(
    const Eigen::MatrixXd &initial_state, const double &transmissibility,
    const double &infectiousness_rate, const double &recovery_rate,
    const Eigen::MatrixXd &contact_matrix,
    const Rcpp::NumericVector &npi_time_begin,
    const Rcpp::NumericVector &npi_time_end, const Rcpp::NumericMatrix &npi_cr,
    const Eigen::MatrixXd &vax_time_begin, const Eigen::MatrixXd &vax_time_end,
    const Eigen::MatrixXd &vax_nu, const Rcpp::List &rate_interventions,
    const Rcpp::List &time_dependence,
    const double &time_end = 100.0,  // double required by boost solver
    const double &increment = 1.0) {
  // initial conditions from input
  odetools::state_type x = initial_state;

  // create a map of the infection parameters
  std::unordered_map<std::string, double> infection_params{
      {"transmissibility", transmissibility},
      {"infectiousness_rate", infectiousness_rate},
      {"recovery_rate", recovery_rate}};

  // create a map of the rate interventions
  const std::unordered_map<std::string, intervention::rate_intervention>
      rate_interventions_cpp =
          intervention::rate_intervention_cpp(rate_interventions);

  // create a default epidemic with parameters
  epidemics::epidemic_default this_model(
      infection_params, contact_matrix, npi_time_begin, npi_time_end, npi_cr,
      vax_time_begin, vax_time_end, vax_nu, rate_interventions_cpp,
      time_dependence);

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
