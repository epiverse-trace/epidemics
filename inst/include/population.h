
// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_POPULATION_H_
#define INST_INCLUDE_POPULATION_H_

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <vector>
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

/// @brief Hold parameters for changes to the model population.
/// This is special functionality for the diphtheria model applied to
/// humanitarian camps.
struct population_change {
  const Rcpp::NumericVector times;
  const Rcpp::List value;  // note these are actually absolute values
  const int n_demo_groups;

  /// @brief Constructor for the population change struct
  /// @param times The times at which the population changes. Changes apply to
  /// the SUSCEPTIBLES compartment only; the assumption is that the camp is the
  /// locus of the outbreak.
  /// @param value The ABSOLUTE change in the population size; an Rcpp List
  /// of Rcpp NumericVectors, each of the same length as the number of
  /// demographic groups.
  population_change(const Rcpp::NumericVector &times, const Rcpp::List &value,
                    const int &n_demo_groups)
      : times(times), value(value), n_demo_groups(n_demo_groups) {}

  // member function for population change at time t
  /// @brief Calculate the population change at time t
  /// @param t The current timestep
  /// @return An Eigen Array of population size changes.
  const Eigen::ArrayXd get_population_change(const double &t);
};

/// @brief Calculate the population change at time t
/// @param t The current timestep
/// @return An Eigen Array of population size changes, which defaults to zeros
/// if no population change is scheduled
inline const Eigen::ArrayXd population_change::get_population_change(
    const double &t) {
  // empty value for timepoints of no population change
  Eigen::ArrayXd pop_change(n_demo_groups);
  pop_change.fill(0.0);

  // crudely check for a match between current time and scheduled pop change
  for (size_t i = 0; i < times.size(); i++) {
    if (t > times[i] && t < (times[i] + 1.0)) {
      pop_change = Rcpp::as<Eigen::ArrayXd>(value[i]);
      return pop_change;
    }
  }

  return pop_change;
}

}  // namespace population

#endif  // INST_INCLUDE_POPULATION_H_
