// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <boost/numeric/odeint.hpp>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

//' @title Run the RIVM Vacamole model
//'
//' @description Vacamole is a deterministic, compartmental epidemic model built
//' by Kylie Ainslie and others at RIVM, the Dutch Public Health Institute for
//' the Covid-19 pandemic, with a focus on scenario modelling for
//' hospitalisation and vaccination.
//' Model code: https://github.com/kylieainslie/vacamole
//' Manuscript describing the model and its application:
//' https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090
//' This function is intended to only be called internally from
//' [model_vacamole_cpp()].
//'
//' @param initial_state A matrix for the initial state of the compartments.
//' @param transmission_rate The transmission rate \eqn{\beta} at which
//' unvaccinated and partially vaccinated individuals are infected by the
//' disease.
//' @param transmission_rate_vax The transmission rate \eqn{\beta_V} at which
//' individuals who have received two vaccine doses are infected by the disease.
//' @param infectiousness_rate The rate of transition from exposed to infectious
//' \eqn{\alpha}.
//' This is common to fully susceptible, partially vaccinated, and fully
//' vaccinated individuals (where fully vaccinated represents two doses).
//' @param mortality_rate The mortality rate of fully susceptible and partially
//' vaccinated and unprotected individuals.
//' @param mortality_rate_vax The mortality rate of individuals who are
//' protected by vaccination.
//' @param hospitalisation_rate The hospitalisation rate of fully susceptible
//' and partially vaccinated and unprotected individuals.
//' @param hospitalisation_rate_vax The hospitalisation rate of individuals who
//' are protected by vaccination.
//' @param recovery_rate The recovery rate \eqn{\gamma}.
//' @param contact_matrix The population contact matrix.
//' @param npi_time_begin The start time of any non-pharmaceutical interventions
//' .
//' @param npi_time_end The end time of any non-pharmaceutical interventions.
//' @param npi_cr The reduction in contacts from any non-pharmaceutical
//' interventions.
//' @param vax_time_begin The start time of any vaccination campaigns.
//' @param vax_time_end The end time of any vaccination campaigns.
//' @param vax_nu The vaccination rate of any vaccination campaigns.
//' @param rate_interventions A named list of `<rate_intervention>` objects.
//' @param time_dependence A named list of functions for parameter time
//' dependence.
//' @param time_end The end time of the simulation.
//' @param increment The time increment of the simulation.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix (see [population()]).
//' The second list element is a vector of timesteps.
//' @keywords internal
// [[Rcpp::export(name=".model_vacamole_cpp")]]
Rcpp::List model_vacamole_internal(
    const Eigen::MatrixXd &initial_state, const double &transmission_rate,
    const double &transmission_rate_vax, const double &infectiousness_rate,
    const double &mortality_rate, const double &mortality_rate_vax,
    const double &hospitalisation_rate, const double &hospitalisation_rate_vax,
    const double &recovery_rate, const Eigen::MatrixXd &contact_matrix,
    const Rcpp::NumericVector &npi_time_begin,
    const Rcpp::NumericVector &npi_time_end, const Rcpp::NumericMatrix &npi_cr,
    const Eigen::MatrixXd &vax_time_begin, const Eigen::MatrixXd &vax_time_end,
    const Eigen::MatrixXd &vax_nu, const Rcpp::List &rate_interventions,
    const Rcpp::List &time_dependence,
    const double &time_end = 100.0,  // double required by boost solver
    const double &increment = 1.0) {
  // initial conditions from input
  odetools::state_type x = initial_state;

  // create a map of the model parameters
  std::unordered_map<std::string, double> model_params{
      {"transmission_rate", transmission_rate},
      {"infectiousness_rate", infectiousness_rate},
      {"mortality_rate", mortality_rate},
      {"hospitalisation_rate", hospitalisation_rate},
      {"recovery_rate", recovery_rate}};

  // create a map of the rate interventions
  const std::unordered_map<std::string, intervention::rate_intervention>
      rate_interventions_cpp =
          intervention::rate_intervention_cpp(rate_interventions);

  // create a default epidemic with parameters
  epidemics::epidemic_vacamole this_model(
      model_params, contact_matrix, npi_time_begin, npi_time_end, npi_cr,
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
