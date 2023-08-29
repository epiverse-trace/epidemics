// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_TIME_DEPENDENCE_H_
#define INST_INCLUDE_TIME_DEPENDENCE_H_

// clang-format off
#include <Rcpp.h>
#include <unordered_map>
#include <string>
// clang-format on

/// @brief A namespace for helpers related to time dependence functions
namespace time_dependence {

/// @brief Apply time-dependence functions to model parameters
/// @param t The time; this argument is passed to the functions in
/// `time_dependence`
/// @param infection_params The infection parameters as key-value pairs
/// @param time_dependence The time-dependence functions as a list of 
/// Rcpp functions whose first argument is time, and whose second argument is
/// the parameter.
/// @return A map of the same size and with the same keys as `infection_params`.
inline const std::unordered_map<std::string, double> apply_time_dependence(
    const double &t,
    const std::unordered_map<std::string, double> &infection_params,
    const Rcpp::List &time_dependence) {
  // make copy of infection params
  std::unordered_map<std::string, double> params_temp = infection_params;

  // get time-dependence names - these are the target parameters
  Rcpp::CharacterVector time_dep_targets = time_dependence.names();

  Rcpp::Rcout << "t = " << t << "\n";

  // loop over rate_interventions and check
  for (size_t i = 0; i < time_dependence.size(); i++) {
    std::string name = Rcpp::as<std::string>(time_dep_targets(i));

    // get function output
    Rcpp::Function f = time_dependence[i];
    double t_mod_param = Rcpp::as<double>(f(t, params_temp.at(name)));

    Rcpp::Rcout << "t_mod_param = " << t_mod_param << "\n";

    params_temp.at(name) = t_mod_param;
  }

  return params_temp;
}

}  // namespace time_dependence

#endif  // INST_INCLUDE_TIME_DEPENDENCE_H_
