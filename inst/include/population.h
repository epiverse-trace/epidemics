
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

}  // namespace population

#endif  // INST_INCLUDE_POPULATION_H_
