// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_EPIDEMIC_DEFAULT_H_
#define INST_INCLUDE_EPIDEMIC_DEFAULT_H_

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
#include "intervention.h"
#include "vaccination.h"
#include "population.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/* The rhs of x' = f(x) defined as a struct with an operator */
struct epidemic_default {
  const double beta, alpha, gamma;
  const Eigen::MatrixXd nu;
  const Rcpp::List population, intervention, vaccination;
  const Eigen::MatrixXd contact_matrix;
  // npi, interv, pop
  epidemic_default(const double beta, const double alpha, const double gamma,
                   const Rcpp::List population, const Rcpp::List intervention,
                   const Rcpp::List vaccination)
      : beta(beta),
        alpha(alpha),
        gamma(gamma),
        nu(vaccination::get_nu(vaccination)),
        population(population),
        intervention(intervention),
        vaccination(vaccination),
        contact_matrix(population::get_contact_matrix(population)) {}

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

    // get current vaccination rate
    Eigen::MatrixXd current_nu = vaccination::current_nu(nu, vaccination, t);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations

    // compartmental transitions without accounting for contacts
    Eigen::ArrayXd sToE = beta * x.col(0).array() * (cm * x.col(2)).array();
    Eigen::ArrayXd eToI = alpha * x.col(1).array();
    Eigen::ArrayXd iToR = gamma * x.col(2).array();
    Eigen::ArrayXd sToV = current_nu.col(0).array() * x.col(0).array();

    // compartmental changes accounting for contacts (for dS and dE)
    dxdt.col(0) = -sToE - sToV;  // -β*S*contacts*I - ν*S
    dxdt.col(1) = sToE - eToI;   // β*S*contacts*I - α*E
    dxdt.col(2) = eToI - iToR;   // α*E - γ*I
    dxdt.col(3) = iToR;          // γ*I
    dxdt.col(4) = sToV;          // ν*S
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_DEFAULT_H_
