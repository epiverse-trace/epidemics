// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_EPIDEMIC_DEFAULT_H_
#define INST_INCLUDE_EPIDEMIC_DEFAULT_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
#include "intervention.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/* The rhs of x' = f(x) defined as a struct with an operator */
struct epidemic_default {
  const Eigen::VectorXd beta, alpha, gamma;
  const Eigen::MatrixXd contact_matrix;
  const Rcpp::List intervention;
  // npi, interv, pop
  epidemic_default(const Eigen::VectorXd beta, const Eigen::VectorXd alpha,
                   const Eigen::VectorXd gamma,
                   const Eigen::MatrixXd contact_matrix,
                   const Rcpp::List intervention)
      : beta(beta),
        alpha(alpha),
        gamma(gamma),
        contact_matrix(contact_matrix),
        intervention(intervention) {}

  void operator()(const odetools::state_type& x,
                  odetools::state_type& dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of x
    dxdt.resize(x.rows(), x.cols());

    // modify contact matrix if time is within intervention timespan
    Eigen::MatrixXd cm = contact_matrix;
    if (intervention::is_intervention_active(t, intervention)) {
      cm = intervention::intervention_on_cm(contact_matrix, intervention);
    }

    // compartmental equations
    dxdt.col(0) = -(beta * x.col(0)) * (cm * x.col(2));  // -beta*S*contacts*I
    dxdt.col(1) = (beta * x.col(0) * (cm * x.col(2))) -
                  (alpha * x.col(1));  // beta*S*contacts*I - alpha*E
    dxdt.col(2) = (alpha * x.col(1)) - (gamma * x.col(2));  // alpha*E - gamma*I
    dxdt.col(3) = gamma * x.col(2);                         // gamma*I
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_DEFAULT_H_
