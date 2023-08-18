// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_INFECTION_H_
#define INST_INCLUDE_INFECTION_H_

// clang-format off
#include <Rcpp.h>
#include <map>
// clang-format on

/// @brief A namespace for helpers related to the infection class
namespace infection {

/// @brief Convert an `Rcpp::List` of infection parameters into a `std::map`
/// @param infection_parameters
/// @return
inline std::map<std::string, double> make_infection_param_map(
    Rcpp::List &infection_parameters) {
  std::map<std::string, double> result;

  Rcpp::CharacterVector keys = infection_parameters.names();
  Rcpp::NumericVector values(infection_parameters);

  for (int i = 0; i < infection_parameters.size(); ++i) {
    result[std::string(keys[i])] = values[i];
  }

  return result;
}

}  // namespace infection

#endif  // INST_INCLUDE_INFECTION_H_
