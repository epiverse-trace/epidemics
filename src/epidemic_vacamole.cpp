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
//'
//' @details The original model has 8 conceptual compartments - four
//' epidemiological compartments (SEIR), three hospitalisation compartments
//' (H, ICU, ICU2H), and death - see the manuscript in Eurosurveillance.
//' Only infected individuals can enter the hospitalisation or death
//' compartments.
//' Vacamole was implemented as a stand-alone R package, and some versions have
//' been used to generate scenarios for the ECDC Covid-19 Scenario Hub.
//'
//' Individuals from the susceptible compartment may be vaccinated partially
//' or fully (assuming a two dose regimen), with only the second dose reducing
//' their probability of being infected, and of being hospitalised or dying.
//'
//' @param population An object of the `population` class, which holds a
//' population contact matrix, a demography vector, and the initial conditions
//' of each demographic group. See [population()].
//' @param beta The transmission rate \eqn{\beta} at which unvaccinated and
//' partially vaccinated individuals are infected by the disease.
//' @param beta_v The transmission rate \eqn{\beta_V} at which individuals who
//' have received two vaccine doses are infected by the disease.
//' @param alpha The rate of transition from exposed to infectious \eqn{\alpha}.
//' This is common to fully susceptible, partially vaccinated, and fully
//' vaccinated individuals (where fully vaccinated represents two doses).
//' @param omega The mortality rate of fully susceptible and partially
//' vaccinated and unprotected individuals.
//' @param omega_v The mortality rate of individuals who are protected by
//' vaccination.
//' @param eta The hospitalisation rate of fully susceptible and partially
//' vaccinated and unprotected individuals.
//' @param eta_v The hospitalisation rate of individuals who are protected by
//' vaccination.
//' @param gamma The recovery rate \eqn{\gamma}.
//' @param time_end The maximum time. See [epidemic()] for default value.
//' @param intervention A non-pharmaceutical intervention applied during the
//' course of the epidemic, with a start and end time, and age-specific effect
//' on contacts. See [intervention()].
//' @param vaccination A vaccination regime followed during the
//' course of the epidemic, with a group- and dose-specific start and end time,
//' and age-specific rates of delivery of first and second doses.
//' See [vaccination()].
//' @param increment The increment time. See [epidemic()] for default value.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix (see [population()]).
//' The second list element is a vector of timesteps.
//' @keywords internal
// [[Rcpp::export(name=".epidemic_vacamole_cpp")]]
Rcpp::List epidemic_vacamole_cpp(
    const Rcpp::List &population, const double &beta, const double &beta_v,
    const double &alpha, const double &omega, const double &omega_v,
    const double &eta, const double &eta_v, const double &gamma,
    const Rcpp::List &intervention, const Rcpp::List &vaccination,
    const double &time_end,  // double required by boost solver
    const double &increment) {
  // initial conditions from input
  odetools::state_type x = odetools::initial_state_from_pop(population);

  // create a default epidemic with parameters
  epidemics::epidemic_vacamole this_model(beta, beta_v, alpha, omega, omega_v,
                                          eta, eta_v, gamma, population,
                                          intervention, vaccination);

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
