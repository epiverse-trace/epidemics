// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_EPIDEMIC_DIPHTHERIA_H_
#define INST_INCLUDE_EPIDEMIC_DIPHTHERIA_H_

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

/// @brief Struct containing the diphtheria epidemic ODE system
struct epidemic_diphtheria {
  // two maps for the model parameters, one for dynamic modification
  const std::unordered_map<std::string, double> model_params;
  std::unordered_map<std::string, double> model_params_temp;
  const std::unordered_map<std::string, intervention::rate_intervention>
      interventions;
  const Rcpp::List time_dependence;

  /// @brief Constructor for the diphtheria epidemic struct
  /// @param model_params An unordered map of string-double pairs, with the
  /// model parameters as keys, and parameter values as values. The
  /// model parameters are:
  /// - transmissibility The transmission rate
  /// - infectiousness_rate The rate at which individuals become infectious
  /// - recovery_rate The recovery rate
  /// - reporting_rate
  /// - prop_hosp
  /// - hosp_entry_rate
  /// - hosp_exit_rate
  /// @param interventions An unordered map of string-intervention pairs. The
  /// keys must refer to parameters in `model_params`. The `intervention`
  /// struct is defined in `inst/include/intervention.h`.
  /// @param time_dependence An Rcpp List with named elements, where each name
  /// is a model parameter (see above), and each element is a function with
  /// the first two arguments being the current simulation time, and x, a value
  /// that is dependent on time (x is supposed to be a model parameter).
  epidemic_diphtheria(
      const std::unordered_map<std::string, double>& model_params,
      const std::unordered_map<std::string, intervention::rate_intervention>&
          interventions,
      const Rcpp::List& time_dependence)
      : model_params(model_params),
        model_params_temp(model_params),
        interventions(interventions),
        time_dependence(time_dependence) {}

  /// @brief Operator for the diphtheria model
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

    // apply time dependence
    model_params_temp = time_dependence::apply_time_dependence(t, model_params,
                                                               time_dependence);

    // rate interventions
    model_params_temp = intervention::intervention_on_rates(
        t, model_params_temp, interventions);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for broadcasting
    // columns are as follows
    // 0|1|2|3|4
    // S|E|I|H|R

    const double total_infections = (x.col(2).array()).sum();

    // compartmental transitions without accounting for stratified contacts
    // NOTE: division by population size - this is typically included in the
    // contact matrix scaling in other models
    Eigen::ArrayXd sToE = model_params_temp["transmissibility"] *
                          total_infections * x.col(0) / x.sum();
    Eigen::ArrayXd eToI =
        model_params_temp["infectiousness_rate"] * x.col(1).array();
    // hospitalised individuals = a function of reporting, p(hospitalised),
    // time to hospitalisation, and time to discharge from hospital
    Eigen::ArrayXd iToH =
        model_params_temp["prop_hosp"] * model_params_temp["reporting_rate"] *
        model_params_temp["hosp_entry_rate"] * x.col(2).array();
    Eigen::ArrayXd iToR = model_params_temp["recovery_rate"] * x.col(2).array();
    Eigen::ArrayXd hToR =
        model_params_temp["hosp_exit_rate"] * x.col(3).array();

    // compartmental changes; note that there are no contacts
    // β: transmissibility; σ: infectiousness rate; γ: recovery rate
    // ν: reporting rate; τ1: 1 / time to hospitalisation;
    // τ2: 1 / time to discharge from hospital; η: prop. hospitalised
    dxdt.col(0) = -sToE;        // -β*S*I
    dxdt.col(1) = sToE - eToI;  // β*S*I - σ*E
    dxdt.col(2) = eToI - iToR;  // σ*E - γ*I
    dxdt.col(3) = iToH - hToR;  // τ1*η*ν*I - τ2*H
    dxdt.col(4) = iToR + hToR;  // γ*I + τ2*H
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_DIPHTHERIA_H_
