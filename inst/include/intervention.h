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
// clang-format on

// add to namespace ode
namespace intervention {

/// @brief
struct intervention {
  double time_begin;
  double time_end;
  double reduction;

  /// @brief A struct that holds intervention parameters
  /// @param time_begin The start time of the intervention
  /// @param time_end The end time of the intervention
  /// @param reduction The proportional reduction in the parameter affected
  intervention(double time_begin, double time_end, double reduction)
      : time_begin(time_begin), time_end(time_end), reduction(reduction) {}

  intervention() : time_begin(0.0), time_end(0.0), reduction(0.0) {}
};

inline std::unordered_map<std::string, intervention> translate_interventions(
    Rcpp::List &interventions) {
  // to hold output
  std::unordered_map<std::string, intervention> result;

  Rcpp::CharacterVector intervention_targets = interventions.names();

  for (size_t i = 0; i < interventions.length(); i++) {
    intervention this_intervention(interventions[i]["time_begin"],
                                   interventions[i]["time_end"],
                                   interventions[i]["reduction)"]);

    std::string name = Rcpp::as<std::string>(intervention_targets[i]);

    result[name] = this_intervention;
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
inline Eigen::ArrayXd cumulative_intervention(
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
  Eigen::ArrayXd contact_reduction =
      cumulative_intervention(t, time_begin, time_end, cr);

  // modify the contact matrix as cm_mod = cm * (1 - intervention)
  // for a percentage reduction in contacts
  Eigen::MatrixXd modified_cm =
      cm.array().colwise() * (1.0 - contact_reduction);
  // transpose for rowwise array multiplication, as Eigen is col-major
  return modified_cm;
}

inline void apply_interventions(
    const double &t, std::unordered_map<std::string, double> &infection_params,
    const std::unordered_map<std::string, intervention> &interventions) {
  // loop over interventions and check
  for (const auto &pair : interventions) {
    intervention temp = pair.second;
    double effect = std::abs(t - temp.time_begin) < 1e-6
                        ? temp.reduction
                        : 0.0;

    infection_params[pair.first] *= effect;
  }
}

}  // namespace intervention

#endif  // INST_INCLUDE_INTERVENTION_H_
