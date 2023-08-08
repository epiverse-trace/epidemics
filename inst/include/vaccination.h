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
/// @param t A double value giving the current time, which is compared against
/// the vaccination start and end times.
/// @param vax_nu An Eigen Matrix of the group- and dose-specific vaccination
/// rates.
/// @param vax_time_begin An Eigen Matrix of the group- and dose-specific
/// vaccination start times.
/// @param vax_time_end An Eigen Matrix of the group- and dose-specific
/// vaccination end times.
/// @return An Eigen Matrix of the same dimensions as `vax_nu`, with values
/// modified to 0.0 if vaccination for that group and dose is not active.
inline Eigen::MatrixXd current_nu(const double &t,
                                  const Eigen::MatrixXd &vax_nu,
                                  const Eigen::MatrixXd &vax_time_begin,
                                  const Eigen::MatrixXd &vax_time_end) {
  // create empty array
  Eigen::MatrixXd nu_active = vax_nu;
  nu_active.fill(0.0);

  // iterate over a potentially 2D matrix
  for (size_t i = 0; i < vax_nu.rows(); i++) {
    for (size_t j = 0; j < vax_nu.cols(); j++) {
      nu_active(i, j) = (t >= vax_time_begin(i, j) && t <= vax_time_end(i, j))
                            ? vax_nu(i, j)
                            : nu_active(i, j);
    }
  }

  return nu_active;
}

}  // namespace vaccination

#endif  // INST_INCLUDE_VACCINATION_H_
