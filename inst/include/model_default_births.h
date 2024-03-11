// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_MODEL_DEFAULT_H_
#define INST_INCLUDE_MODEL_DEFAULT_H_

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <unordered_map>
#include <string>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
#include "intervention.h"
#include "vaccination.h"
#include "population.h"
#include "time_dependence.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/* The rhs of x' = f(x) defined as a struct with an operator */

/// @brief Struct containing the default epidemic ODE system
struct epidemic_default {
  // two maps for the model parameters, one for dynamic modification
  const std::unordered_map<std::string, double> model_params;
  std::unordered_map<std::string, double> model_params_temp;
  const Eigen::MatrixXd contact_matrix;
  Eigen::MatrixXd cm_temp;
  // related to interventions
  const Rcpp::NumericVector npi_time_begin, npi_time_end;
  const Rcpp::NumericMatrix npi_cr;
  // related to vaccination
  const Eigen::MatrixXd vax_time_begin, vax_time_end, vax_nu;
  Eigen::MatrixXd vax_nu_current;
  const std::unordered_map<std::string, intervention::rate_intervention>
      interventions;
  const Rcpp::List time_dependence;

  /// @brief Constructor for the default epidemic struct
  /// @param model_params An unordered map of string-double pairs, with the
  /// model parameters as keys, and parameter values as values. The
  /// model parameters are:
  /// - transmissibility The transmission rate
  /// - infectiousness_rate The rate at which individuals become infectious
  /// - recovery_rate The recovery rate
  /// @param contact_matrix The population contact matrix
  /// @param npi_time_begin The intervention start times
  /// @param npi_time_end The intervention end times
  /// @param npi_cr The intervention contact reduction
  /// @param vax_time_begin The age- and dose-specific vaccination start time
  /// @param vax_time_end The age- and dose-specific vaccination end time
  /// @param vax_nu The age- and dose-specific vaccination rate
  /// @param interventions An unordered map of string-intervention pairs. The
  /// keys must refer to parameters in `model_params`. The `intervention`
  /// struct is defined in `inst/include/intervention.h`.
  /// @param time_dependence An Rcpp List with named elements, where each name
  /// is a model parameter (see above), and each element is a function with
  /// the first two arguments being the current simulation time, and x, a value
  /// that is dependent on time (x is supposed to be a model parameter).
  epidemic_default(
      const std::unordered_map<std::string, double>& model_params,
      const Eigen::MatrixXd contact_matrix,
      const Rcpp::NumericVector npi_time_begin,
      const Rcpp::NumericVector npi_time_end, const Rcpp::NumericMatrix npi_cr,
      const Eigen::MatrixXd vax_time_begin, const Eigen::MatrixXd vax_time_end,
      const Eigen::MatrixXd vax_nu,
      const std::unordered_map<std::string, intervention::rate_intervention>&
          interventions,
      const Rcpp::List& time_dependence)
      : model_params(model_params),
        model_params_temp(model_params),
        contact_matrix(contact_matrix),
        cm_temp(contact_matrix),
        npi_time_begin(npi_time_begin),
        npi_time_end(npi_time_end),
        npi_cr(npi_cr),
        vax_time_begin(vax_time_begin),
        vax_time_end(vax_time_end),
        vax_nu(vax_nu),
        vax_nu_current(vax_nu),
        interventions(interventions),
        time_dependence(time_dependence) {}

  /// @brief Operator for the default model
  /// @param x The initial state of the population - rows represent age groups
  /// while columns represent compartments
  /// @param dxdt An object of the same type as `x` to hold the current state of
  /// the system
  /// @param t The simulation time
  void operator()(const odetools::state_type& x,
                  odetools::state_type& dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of x
    dxdt.resize(x.rows(), x.cols());

    // modify contact matrix if time is within intervention timespan
    cm_temp = intervention::intervention_on_cm(
        t, contact_matrix, npi_time_begin, npi_time_end, npi_cr);

    // apply time dependence
    model_params_temp = time_dependence::apply_time_dependence(t, model_params,
                                                               time_dependence);

    // rate interventions
    model_params_temp = intervention::intervention_on_rates(
        t, model_params_temp, interventions);

    // get current vaccination rate
    vax_nu_current =
        vaccination::current_nu(t, vax_nu, vax_time_begin, vax_time_end);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations

    // compartmental transitions without accounting for contacts
    Eigen::ArrayXd sToE = model_params_temp["transmissibility"] *
                          x.col(0).array() * (cm_temp * x.col(2)).array();
    Eigen::ArrayXd eToI =
        model_params_temp["infectiousness_rate"] * x.col(1).array();
    Eigen::ArrayXd iToR = model_params_temp["recovery_rate"] * x.col(2).array();
    Eigen::ArrayXd sToV = vax_nu_current.col(0).array() * x.col(0).array();

    // compartmental changes accounting for contacts (for dS and dE)
    // β: transmissibility; ν: vaccination rate; σ: infectiousness rate
    // γ: recovery rate
    dxdt.col(0) = -sToE - sToV;  // -β*S*contacts*I - ν*S
    dxdt.col(1) = sToE - eToI;   // β*S*contacts*I - σ*E
    dxdt.col(2) = eToI - iToR;   // σ*E - γ*I
    dxdt.col(3) = iToR;          // γ*I
    dxdt.col(4) = sToV;          // ν*S
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_MODEL_DEFAULT_H_
