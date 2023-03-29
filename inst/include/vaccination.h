// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_VACCINATION_H_
#define INST_INCLUDE_VACCINATION_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <iostream>
#include <vector>
// clang-format on

// add to namespace intervention
namespace vaccination {

inline Eigen::ArrayXd current_nu(const Eigen::ArrayXd &nu,
                                 const Eigen::ArrayXd &t_vax_begin,
                                 const Eigen::ArrayXd &t_vax_end,
                                 const double &t) {
  Eigen::ArrayXd nu_active(nu.size());
  nu_active.fill(0.0);

  for (size_t i = 0; i < nu.size(); i++) {
    nu_active[i] =
        (t >= t_vax_begin[i] && t <= t_vax_end[i]) ? nu[i] : nu_active[i];
  }

  return nu_active;
}

}  // namespace vaccination

#endif  // INST_INCLUDE_VACCINATION_H_
