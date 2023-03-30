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
/// @return An Eigen Array of group-specific vaccination rates.
inline Eigen::ArrayXd get_nu(const Rcpp::List &vaccination) {
  return Rcpp::as<Eigen::ArrayXd>(vaccination["nu"]);
}

/// @brief
/// @param nu
/// @param vaccination
/// @param t
/// @return
inline Eigen::ArrayXd current_nu(const Eigen::ArrayXd &nu,
                                 const Rcpp::List &vaccination,
                                 const double &t) {
  // create Eigen 1D array from R matrix passed in an list (class intervention)
  Eigen::ArrayXd time_begin(
      Rcpp::as<Eigen::ArrayXd>(vaccination["time_begin"]));
  Eigen::ArrayXd time_end(Rcpp::as<Eigen::ArrayXd>(vaccination["time_end"]));

  Eigen::ArrayXd nu_active(nu.size());
  nu_active.fill(0.0);

  for (size_t i = 0; i < nu.size(); i++) {
    nu_active[i] =
        (t >= time_begin[i] && t <= time_end[i]) ? nu[i] : nu_active[i];
  }

  return nu_active;
}

}  // namespace vaccination

#endif  // INST_INCLUDE_VACCINATION_H_
