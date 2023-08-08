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
  const Eigen::MatrixXd contact_matrix;
  Eigen::MatrixXd cm_temp;
  // related to interventions
  const Rcpp::NumericVector npi_time_begin, npi_time_end;
  const Rcpp::NumericMatrix npi_cr;
  // related to vaccination
  const Eigen::MatrixXd vax_time_begin, vax_time_end, vax_nu;
  Eigen::MatrixXd vax_nu_current;

  // npi, interv, pop
  epidemic_default(const double beta, const double alpha, const double gamma,
                   const Eigen::MatrixXd contact_matrix,
                   const Rcpp::NumericVector npi_time_begin,
                   const Rcpp::NumericVector npi_time_end,
                   const Rcpp::NumericMatrix npi_cr,
                   const Eigen::MatrixXd vax_time_begin,
                   const Eigen::MatrixXd vax_time_end,
                   const Eigen::MatrixXd vax_nu)
      : beta(beta),
        alpha(alpha),
        gamma(gamma),
        contact_matrix(contact_matrix),
        cm_temp(contact_matrix),
        npi_time_begin(npi_time_begin),
        npi_time_end(npi_time_end),
        npi_cr(npi_cr),
        vax_time_begin(vax_time_begin),
        vax_time_end(vax_time_end),
        vax_nu(vax_nu),
        vax_nu_current(vax_nu) {}

  void operator()(const odetools::state_type& x,
                  odetools::state_type& dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of x
    dxdt.resize(x.rows(), x.cols());

    // modify contact matrix if time is within intervention timespan
    cm_temp = intervention::intervention_on_cm(
        t, contact_matrix, npi_time_begin, npi_time_end, npi_cr);

    // get current vaccination rate
    vax_nu_current =
        vaccination::current_nu(t, vax_nu, vax_time_begin, vax_time_end);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations

    // compartmental transitions without accounting for contacts
    Eigen::ArrayXd sToE =
        beta * x.col(0).array() * (cm_temp * x.col(2)).array();
    Eigen::ArrayXd eToI = alpha * x.col(1).array();
    Eigen::ArrayXd iToR = gamma * x.col(2).array();
    Eigen::ArrayXd sToV = vax_nu_current.col(0).array() * x.col(0).array();

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
