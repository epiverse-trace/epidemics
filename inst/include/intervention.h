// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_INTERVENTION_H_
#define INST_INCLUDE_INTERVENTION_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>
// clang-format on

// add to namespace ode
namespace intervention {

/// @brief A struct for rate interventions
struct rate_intervention {
  std::vector<double> time_begin;
  std::vector<double> time_end;
  std::vector<double> reduction;
  int n_interventions;

  /// @brief Constructor for the rate_intervention struct
  /// @param time_begin The start time of the intervention
  /// @param time_end The end time of the intervention
  /// @param reduction The proportional reduction in the parameter affected
  rate_intervention(const std::vector<double> &time_begin,
                    const std::vector<double> &time_end,
                    const std::vector<double> &reduction)
      : time_begin(time_begin),
        time_end(time_end),
        reduction(reduction),
        n_interventions(time_begin.size()) {}

  /// @brief Constructor for the rate_intervention struct
  rate_intervention()
      : time_begin({0.0}),
        time_end({0.0}),
        reduction({0.0}),
        n_interventions(0) {}
};

inline std::unordered_map<std::string, rate_intervention> rate_intervention_cpp(
    const Rcpp::List &rate_interventions) {
  // to hold output
  std::unordered_map<std::string, rate_intervention> result;

  // get intervention names - these are the targets
  Rcpp::CharacterVector intervention_targets = rate_interventions.names();

  // add intervention objects to map
  for (size_t i = 0; i < rate_interventions.size(); i++) {
    // make a copy of the list element, which is also a list
    Rcpp::List intervention_element = rate_interventions[i];

    // make copies of the intervention members
    std::vector<double> time_begin =
        Rcpp::as<std::vector<double> >(intervention_element["time_begin"]);
    std::vector<double> time_end =
        Rcpp::as<std::vector<double> >(intervention_element["time_end"]);
    std::vector<double> reduction =
        Rcpp::as<std::vector<double> >(intervention_element["reduction"]);

    // create new intervention object from list element members
    rate_intervention intervention_temp(time_begin, time_end, reduction);

    // add to map
    std::string name = Rcpp::as<std::string>(intervention_targets[i]);
    result[name] = intervention_temp;
  }

  return result;
}

/// @brief Get the cumulative effect of interventions
/// @param t The current simulation time.
/// @param time_begin The time for each intervention to begin.
/// @param time_end The time for each intervention to end.
/// @param cr A matrix with the demographic group and intervention specific
/// effect. When two interventions are simultaneously active, their cumulative
/// effect is additive, that is, the two interventions' effects on social
/// contacts are added together.
/// @return An Eigen Array with as many elements as there are demographic groups
/// in the population. This is the number of rows of `cr`.
/// The array gives the current cumulative effect of interventions on the
/// corresponding age group.
inline Eigen::ArrayXd cumulative_contacts_intervention(
    const double &t, const Rcpp::NumericVector &time_begin,
    const Rcpp::NumericVector &time_end, const Rcpp::NumericMatrix &cr) {
  // a vector with as elements as the number of rows, i.e., age groups
  Rcpp::NumericVector eff_con_red(cr.nrow());

  // iterate over the Rcpp vector time_begin, a member of the list-like
  // intervention-class
  for (size_t i = 0; i < time_begin.size(); i++) {
    if (t >= time_begin(i) && t <= time_end(i)) {
      eff_con_red += cr(Rcpp::_, i);
    }
  }

  // correct for contact reductions greater than 1.0
  for (size_t i = 0; i < eff_con_red.size(); i++) {
    eff_con_red[i] = std::min(1.0, eff_con_red[i]);
  }

  // transform to an Eigen Array and return
  Eigen::ArrayXd effective_contact_reduction(
      Rcpp::as<Eigen::ArrayXd>(eff_con_red));

  return effective_contact_reduction;
}

/// @brief Get the contact matrix modified by interventions
/// @param t The current simulation time.
/// @param cm The population contact matrix.
/// @param time_begin The time for each intervention to begin.
/// @param time_end The time for each intervention to end.
/// @param cr A matrix with the demographic group and intervention specific
/// effect. When two interventions are simultaneously active, their cumulative
/// effect is additive, that is, the two interventions' effects on social
/// contacts are added together.
/// @return An Eigen Matrix of the same dimensions as `cm`.
inline Eigen::MatrixXd intervention_on_cm(const double &t,
                                          const Eigen::MatrixXd &cm,
                                          const Rcpp::NumericVector &time_begin,
                                          const Rcpp::NumericVector &time_end,
                                          const Rcpp::NumericMatrix &cr) {
  // create Eigen 1D array from R matrix passed in an list (class intervention)
  Eigen::ArrayXd contact_scaling =
      1.0 - cumulative_contacts_intervention(t, time_begin, time_end, cr);

  // modify the contact matrix by multiplying each i-th row and column
  // with the i-th element of the contact reduction array
  // first multiply the rows. THIS REQUIRES COLWISE as Eigen is col-major!
  Eigen::MatrixXd modified_cm = cm.array().colwise() * contact_scaling;
  // then multiply the columns. THIS USES ROWWISE with array transpose.
  modified_cm = modified_cm.array().rowwise() * contact_scaling.transpose();
  return modified_cm;
}

/// @brief Get the cumulative effect of interventions
/// @param t The current simulation time.
/// @param rate_intervention A `rate_intervention` struct with the members:
/// - time_begin The time for each intervention to begin.
/// - time_end The time for each intervention to end.
/// - reduction A vector with the intervention-specific effects on a model
/// rate parameter. When two interventions are simultaneously active, their
/// cumulative effect is additive, that is, the two interventions' effects on
/// a model rate parameter are added together.
/// @return An single double value for the proportional reduction in a model
/// rate parameter. Values > 1.0 are capped at 1.0.
inline double cumulative_rate_intervention(
    const double &t, const rate_intervention &rate_interv) {
  // a vector with as elements as the number of rows, i.e., age groups
  double effect = 0.0;

  // iterate over the rate interventions
  for (size_t i = 0; i < rate_interv.n_interventions; i++) {
    if (t >= rate_interv.time_begin[i] && t <= rate_interv.time_end[i]) {
      effect += rate_interv.reduction[i];
    }
  }

  // correct for reductions greater than 1.0
  effect = effect > 1.0 ? 1.0 : effect;

  return effect;
}

/// @brief Apply rate_interventions on the rate parameters
/// @param t The current simulation time
/// @param infection_params A map of the infection rate parameters
/// @param rate_interventions A map of the rate_interventions
inline std::unordered_map<std::string, double> intervention_on_rates(
    const double &t,
    const std::unordered_map<std::string, double> &infection_params,
    const std::unordered_map<std::string, rate_intervention>
        &rate_interventions) {
  // make copy of infection params
  std::unordered_map<std::string, double> params_temp = infection_params;

  // loop over rate_interventions and check
  for (const auto &pair : rate_interventions) {
    // assign the cumulative effect of each rate intervention to the rate params
    params_temp.at(pair.first) *=
        (1.0 - cumulative_rate_intervention(t, pair.second));
  }

  return params_temp;
}

}  // namespace intervention

#endif  // INST_INCLUDE_INTERVENTION_H_
