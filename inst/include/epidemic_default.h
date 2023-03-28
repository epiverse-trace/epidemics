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
  const Eigen::ArrayXd beta, alpha, gamma;
  const Eigen::MatrixXd contact_matrix;
  const Rcpp::List intervention;
  // npi, interv, pop
  epidemic_default(const Eigen::ArrayXd beta, const Eigen::ArrayXd alpha,
                   const Eigen::ArrayXd gamma,
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

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for correct group-specific coefficient multiplications

    // compartmental transitions without accounting for contacts
    Eigen::VectorXd sToE = beta * x.col(0).array() * x.col(2).array();
    Eigen::VectorXd eToI = alpha * x.col(1).array();
    Eigen::VectorXd iToR = gamma * x.col(2).array();

    // compartmental changes accounting for contacts (for dS and dE)
    dxdt.col(0) = -cm * sToE;                          // -contacts * β*S*I
    dxdt.col(1) = (cm * sToE).array() - eToI.array();  // β*S*contacts*I - α*E
    dxdt.col(2) = eToI.array() - iToR.array();         // α*E - γ*I
    dxdt.col(3) = iToR;                                // γ*I
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_DEFAULT_H_
