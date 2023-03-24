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
//' @param init The initial conditions, represented as a matrix, in which the
//' rows \eqn{i} represent age or demographic groups, and columns \eqn{j}
//' represent the proportions of each age group in the epidemiological
//' compartment. Compartments are arranged in order: "S", "E", "I", and "R".
//' Multiple (N) age groups are currently supported, thus `init` must be
//' a \eqn{N \times 4} matrix.
//' @param contact_matrix A matrix whose elements \eqn{i,j} give the number of
//' contacts between individuals of groups \eqn{i} and \eqn{j}.
//' @param beta The transmission rate \eqn{\beta}.
//' @param alpha The rate of transition from exposed to infectious \eqn{\alpha}.
//' @param gamma The recovery rate \eqn{\gamma}.
//' @param time_end The maximum time, defaults to 200.0.
//' @param increment The increment time, defaults to 0.1.
//' @export
// [[Rcpp::export]]
Rcpp::List epidemic_default_cpp(
    const Eigen::MatrixXd &init, const Eigen::MatrixXd &contact_matrix,
    const float &beta, const float &alpha, const float &gamma,
    const double &time_end = 200.0,  // double required by boost solver
    const double &increment = 0.1) {
  // initial conditions from input
  odetools::state_type x = init;

  // create a default epidemic with parameters
  epidemics::epidemic_default this_model(beta, alpha, gamma, contact_matrix);

  //[ integrate_observ
  std::vector<odetools::state_type> x_vec;  // is a vector of double vectors
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
