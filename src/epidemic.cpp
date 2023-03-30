// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <boost/numeric/odeint.hpp>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

//' @title Run an age-structured SEIR epidemic model
//'
//' @param population An object of the `population` class, which holds a
//' population contact matrix, a demography vector, and the initial conditions
//' of each demographic group. See [population()].
//' @param beta The transmission rate \eqn{\beta}.
//' @param alpha The rate of transition from exposed to infectious \eqn{\alpha}.
//' @param gamma The recovery rate \eqn{\gamma}.
//' @param time_end The maximum time, defaults to 200.0.
//' @param intervention A non-pharamaceutical intervention applied during the
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
// [[Rcpp::export(name=".epidemic_default_cpp")]]
Rcpp::List epidemic_default_cpp(
    const Rcpp::List &population, const Eigen::ArrayXd &beta,
    const Eigen::ArrayXd &alpha, const Eigen::ArrayXd &gamma,
    const Rcpp::List &intervention, const Rcpp::List &vaccination,
    const double &time_end = 200.0,  // double required by boost solver
    const double &increment = 0.1) {
  // initial conditions from input
  odetools::state_type x = odetools::initial_state_from_pop(population);

  // create a default epidemic with parameters
  epidemics::epidemic_default this_model(beta, alpha, gamma, population,
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
