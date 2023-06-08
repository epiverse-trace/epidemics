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

inline Eigen::ArrayXd cumulative_intervention(const double &t,
                                              const Rcpp::List &intervention) {
  // copy the list elements
  // TODO(all): reduce copying
  Rcpp::NumericVector time_begin = intervention["time_begin"];
  Rcpp::NumericVector time_end = intervention["time_end"];
  Rcpp::NumericMatrix cr = intervention["contact_reduction"];

  Rcpp::NumericVector eff_con_red(cr.size());

  // iterate over the Rcpp vector time_begin, a member of the list-like
  // intervention-class
  for (size_t i = 0; i < time_begin.size(); i++) {
    if (t >= time_begin(i) && t <= time_end(i)) {
      eff_con_red += cr(Rcpp::_, i);
    }
  }

  Eigen::ArrayXd effective_contact_reduction(
      Rcpp::as<Eigen::ArrayXd>(eff_con_red));

  return effective_contact_reduction;
}

inline Eigen::MatrixXd intervention_on_cm(const double &t,
                                          const Eigen::MatrixXd &cm,
                                          const Rcpp::List &intervention) {
  // create Eigen 1D array from R matrix passed in an list (class intervention)
  Eigen::ArrayXd contact_reduction = cumulative_intervention(t, intervention);

  // modify the contact matrix as cm_mod = cm * (1 - intervention)
  // for a percentage reduction in contacts
  Eigen::MatrixXd modified_cm =
      cm.array().colwise() * (1.0 - contact_reduction);
  // transpose for rowwise array multiplication, as Eigen is col-major
  return modified_cm;
}

}  // namespace intervention

#endif  // INST_INCLUDE_INTERVENTION_H_
