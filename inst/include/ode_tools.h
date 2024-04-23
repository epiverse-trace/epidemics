// Copyright 2024 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_ODE_TOOLS_H_
#define INST_INCLUDE_ODE_TOOLS_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <vector>
// clang-format on

// add to namespace ode
namespace odetools {

//[ rhs_function
/* The type of container used to hold the state vector */
typedef Eigen::MatrixXd state_type;  // flexible

//[ integrate_observer
struct observer {
  std::vector<state_type> &m_states;
  std::vector<double> &m_times;

  observer(std::vector<state_type> &states,  // NOLINT
           std::vector<double> &times)       // NOLINT
      : m_states(states), m_times(times) {}

  void operator()(const state_type &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

/// @brief Prepare a state_type object from population initial conditions
/// @param population An object of the R S3 class `population`, which is
/// handled as an Rcpp List.
/// @return A `state_type` Eigen::MatriXd that corresponds to population initial
/// conditions.
inline state_type initial_state_from_pop(const Rcpp::List &population) {
  // TODO(pratikunterwegs): input checking here
  return Eigen::MatrixXd(
      Rcpp::as<Eigen::MatrixXd>(population["initial_conditions"]));
}
//]

}  // namespace odetools

#endif  // INST_INCLUDE_ODE_TOOLS_H_
