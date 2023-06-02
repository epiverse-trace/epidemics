// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_INTERVENTION_H_
#define INST_INCLUDE_INTERVENTION_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>
// clang-format on

// add to namespace ode
namespace intervention {

inline Eigen::MatrixXd intervention_on_cm(const Eigen::MatrixXd &cm,
                                          const Rcpp::List &intervention) {
  // create Eigen 1D array from R matrix passed in an list (class intervention)
  Eigen::ArrayXd contact_reduction(
      Rcpp::as<Eigen::ArrayXd>(intervention["contact_reduction"]));

  // modify the contact matrix as cm_mod = cm * (1 - intervention)
  // for a percentage reduction in contacts
  Eigen::MatrixXd modified_cm =
      cm.array().colwise() * (1.0 - contact_reduction);
  // transpose for rowwise array multiplication, as Eigen is col-major
  return modified_cm;
}

inline bool is_intervention_active(const double &t,
                                   const Rcpp::List &intervention) {
  const double time_begin = intervention["time_begin"];
  const double time_end = intervention["time_end"];
  // input checking on intervention
  bool intervention_active = t >= time_begin && t <= time_end;
  return intervention_active;
}

}  // namespace intervention

#endif  // INST_INCLUDE_INTERVENTION_H_
