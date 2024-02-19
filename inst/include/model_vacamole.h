// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
/*
    Adapted from the Vacamole model developed by Kylie Ainslie et al. at the
    National Institute for Public Health and the Environment (RIVM),
    The Netherlands, during the Covid-19 pandemic.
    Source code for the original Vacamole model was developed at:
    https://github.com/kylieainslie/vacamole
    and the model is described in Eurosurveillance at:
    https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090
*/
#ifndef INST_INCLUDE_MODEL_VACAMOLE_H_
#define INST_INCLUDE_MODEL_VACAMOLE_H_

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <string>
#include <unordered_map>

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

/// @brief Struct containing the Vacamole epidemic ODE system
struct epidemic_vacamole {
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

  /// @brief Constructor for the Vacamole epidemic struct
  /// @param model_params An unordered map of string-double pairs, with the
  /// model parameters as keys, and parameter values as values. The
  /// model parameters are:
  /// - transmissibility The transmission rate for un-or-single vaccinated
  /// individuals
  /// - transmissibility_vax The transmission rate for double-vaccinated
  /// individuals
  /// - infectiousness_rate The rate at which individuals become infectious
  /// - mortality_rate The mortality rate of un-or-single-vaccinated individuals
  /// - mortality_rate_vax The mortality rate of double-vaccinated individuals
  /// - hospitalisation_rate The hospitalisation rate of un-or-single-vaccinated
  /// individuals
  /// - hospitalisation_rate_vax The hospitalisation rate of double-vaccinated
  /// individuals
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
  epidemic_vacamole(
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
    // columns are as follows
    // 0| 1| 2|3| 4|5| 6|7| 8|9|10
    // S|V1|V2|E|EV|I|IV|H|HV|D|R

    // compartmental transitions without accounting for contacts
    // Susceptible (unvaccinated) to exposed
    Eigen::ArrayXd sToE = model_params_temp["transmissibility"] *
                          x.col(0).array() *
                          (cm_temp * (x.col(5) + x.col(6))).array();

    // Susceptible to vaccinated with one dose
    Eigen::ArrayXd sToV1 = vax_nu_current.col(0).array() * x.col(0).array();
    // Vaccinated one dose to vaccinated with two doses
    Eigen::ArrayXd v1ToV2 = vax_nu_current.col(1).array() * x.col(1).array();

    // Vaccinated one dose to exposed - same as susceptible to exposed
    Eigen::ArrayXd v1ToE = model_params_temp["transmissibility"] *
                           x.col(1).array() *
                           (cm_temp * (x.col(5) + x.col(6))).array();
    // Vaccinated two doses to exposed - uses different transmissibility
    Eigen::ArrayXd v2ToEv = model_params_temp["transmissibility_vax"] *
                            x.col(2).array() *
                            (cm_temp * (x.col(5) + x.col(6))).array();

    // Exposed unvaccinated or not protected to infectious
    Eigen::ArrayXd eToI =
        model_params_temp["infectiousness_rate"] * x.col(3).array();
    // Exposed vaccinated to infectious
    Eigen::ArrayXd evToIv =
        model_params_temp["infectiousness_rate"] * x.col(4).array();

    // Infectious to hospitalised
    Eigen::ArrayXd iToH =
        model_params_temp["hospitalisation_rate"] * x.col(5).array();
    // Vaccinated infectious to hospitalised
    Eigen::ArrayXd ivToHv =
        model_params_temp["hospitalisation_rate_vax"] * x.col(6).array();

    // Infectious to dead
    Eigen::ArrayXd iToD =
        model_params_temp["mortality_rate"] * x.col(5).array();
    // Infectious vaccinated to dead
    Eigen::ArrayXd ivToD =
        model_params_temp["mortality_rate_vax"] * x.col(6).array();

    // Hospitalised to dead
    Eigen::ArrayXd hToD =
        model_params_temp["mortality_rate"] * x.col(7).array();
    // Hospitalised vaccinated to dead
    Eigen::ArrayXd hvToD =
        model_params_temp["mortality_rate_vax"] * x.col(8).array();

    // Infectious to recovered
    Eigen::ArrayXd iToR = model_params_temp["recovery_rate"] * x.col(5).array();
    // Infectious vaccinated to recovered
    Eigen::ArrayXd ivToR =
        model_params_temp["recovery_rate"] * x.col(6).array();

    // Hospitalised to recovered
    Eigen::ArrayXd hToR = model_params_temp["recovery_rate"] * x.col(7).array();
    // Hospitalised vaccinated to recovered
    Eigen::ArrayXd hvToR =
        model_params_temp["recovery_rate"] * x.col(8).array();

    // compartmental changes accounting for contacts
    // β: transmissibility; βv: transmissibility for doubly vaccinated;
    // ν1, ν2: vaccination rate first, second dose;
    // σ: infectiousness rate; γ: recovery rate;
    // η: hospitalisation rate; η_v: hosp. rate doubly vaccinated
    // ω: mortality rate; ω_v: mort. rate doubly vaccinated
    dxdt.col(0) = -sToE - sToV1;           // -β*S*contacts*(I+Iv) - ν1*S
    dxdt.col(1) = sToV1 - v1ToE - v1ToV2;  // ν1*S -β*V1*contacts*(I+Iv) - ν2*V
    dxdt.col(2) = v1ToV2 - v2ToEv;         // ν2*V - βv*V2*contacts*(I+Iv)
    dxdt.col(3) = sToE + v1ToE - eToI;     // β*(S+V1)*contacts*(I+Iv) - σE
    dxdt.col(4) = v2ToEv - evToIv;         // βv*V2*contacts*(I+Iv) - αE
    dxdt.col(5) = eToI - iToH - iToR - iToD;        // σE - I(γ + η + ω)
    dxdt.col(6) = evToIv - ivToHv - ivToR - ivToD;  // σEv - Iv(γ + η_v + ω_v)
    dxdt.col(7) = iToH - hToR - hToD;               // ηI - γH - ωH
    dxdt.col(8) = ivToHv - hvToR - hvToD;           // ηIv - γHv - ωHv
    dxdt.col(9) = iToD + ivToD + hToD + hvToD;      // all deaths
    dxdt.col(10) = iToR + ivToR + hToR + hvToR;     // all recoveries
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_MODEL_VACAMOLE_H_
