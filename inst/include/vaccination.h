// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_VACCINATION_H_
#define INST_INCLUDE_VACCINATION_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>
// clang-format on

// add to namespace vaccination
namespace vaccination {

/// @brief Get the group specific vaccination rate from a `vaccination` object
/// @param vaccination A `vaccination` object in the form of an Rcpp List
/// @return An Eigen Array of group-specific vaccination rates, with dimensions
/// matching those the `vaccination` class member `nu`, in which each element
/// `i,j` gives the group- and dose-specific vaccination rate, in row `i` and
/// column `j` respectively.
inline Eigen::MatrixXd get_nu(const Rcpp::List &vaccination) {
  // convert to Eigen array, then resize to match the nu matrix
  Eigen::MatrixXd nu(Rcpp::as<Eigen::MatrixXd>(vaccination["nu"]));
  return nu;
}

/// @brief Get the vaccination rate at a time in the vaccination regime
/// @param nu An Eigen Matrix of the group- and dose-specific vaccination rates
/// @param vaccination An Rcpp List object that holds matrices giving the group-
/// and dose-specific vaccination start and end times.
/// @param t A double value giving the current time, which is compared against
/// the vaccination start and end times.
/// @return An Eigen Matrix of the same dimensions as `nu`, with values modified
/// to 0.0 if vaccination for that group and dose is not active.
inline Eigen::MatrixXd current_nu(const Eigen::MatrixXd &nu,
                                  const Rcpp::List &vaccination,
                                  const double &t) {
  // create Eigen 2D array from R matrix passed in an list (class intervention)
  // resize the array to have the correct dimensions as default is 1D
  Eigen::MatrixXd time_begin(
      Rcpp::as<Eigen::MatrixXd>(vaccination["time_begin"]));
  Eigen::MatrixXd time_end(Rcpp::as<Eigen::MatrixXd>(vaccination["time_end"]));

  Eigen::MatrixXd nu_active = nu;
  nu_active.fill(0.0);

  // iterate over a potentially 2D matrix
  for (size_t i = 0; i < nu.rows(); i++) {
    for (size_t j = 0; j < nu.cols(); j++) {
      nu_active(i, j) = (t >= time_begin(i, j) && t <= time_end(i, j))
                            ? nu(i, j)
                            : nu_active(i, j);
    }
  }

  return nu_active;
}

}  // namespace vaccination

#endif  // INST_INCLUDE_VACCINATION_H_
