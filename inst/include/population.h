
// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_POPULATION_H_
#define INST_INCLUDE_POPULATION_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>
// clang-format on

// add to namespace population
namespace population {

/// @brief Get the contact matrix from a `population` object
/// @param population A `population` class object, handled here as an Rcpp List.
/// @return An Eigen::Matrix object corresponding to the population contact
/// matrix.
inline Eigen::MatrixXd get_contact_matrix(const Rcpp::List &population) {
  return Rcpp::as<Eigen::MatrixXd>(population["contact_matrix"]);
}

/// @brief Get the total population size
/// @param population A `population` class object, handled as an Rcpp::List.
/// @return An integer value for the total population size.
inline int get_population_size(const Rcpp::List &population) {
  Rcpp::NumericVector demography_vector = population["demography_vector"];
  return static_cast<int>(Rcpp::sum(demography_vector));
}

/// @brief Get the population initial conditions
/// @param population A `population` class object, handled as an Rcpp::List.
/// @return An Rcpp::NumericMatrix of the initial conditions per age group. See
/// the documentation in R for the `population` class.
inline Rcpp::NumericMatrix get_initial_conditions(
    const Rcpp::List &population) {
  Rcpp::NumericMatrix initial_conditions = population["initial_conditions"];
  return initial_conditions;
}

}  // namespace population

#endif  // INST_INCLUDE_POPULATION_H_
