// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
/*
    Adapted from a stochastic discrete-time model for ebola virus disease
    developed by Wayne Getz and Eric Dougherty (Getz and Dougherty 2017),
    with code adapted by Hạ Minh Lâm for the Epirecipes Cookbook.
    Source code for the model in Epirecipes is at:
    http://epirecip.es/epicookbook/chapters/erlang/r
    and the model is described in the Journal of Biological Dynamics at:
    https://doi.org/10.1080/17513758.2017.1401677
*/
#ifndef INST_INCLUDE_EPIDEMIC_EBOLA_H_
#define INST_INCLUDE_EPIDEMIC_EBOLA_H_

// clang-format off
#include <vector>
#include <utility>
#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
#include "intervention.h"
#include "vaccination.h"
#include "population.h"
#include "helpers.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/// @brief Run a stochastic SEIR model of ebola with Erlang passage times
/// @param initial_conditions A vector representing the number of individuals in
/// each compartment.
/// @param population_size The population size
/// @param beta The transmission rate.
/// @param shape_E The shape of the Erlang distribution of passage times through
/// the exposed compartment
/// @param rate_E The rate of the Erlang distribution of passage times through
/// the exposed compartment
/// @param shape_I The shape of the Erlang distribution of passage times through
/// the infectious compartment
/// @param rate_I The rate of the Erlang distribution of passage times through
/// the infectious compartment
/// @param max_time The time at which the simulation ends
/// @return An Rcpp::List with the states and times
inline Rcpp::List epidemic_ebola(const Rcpp::IntegerVector &initial_conditions,
                                 const int &population_size, const double &beta,
                                 const int &shape_E, const double &rate_E,
                                 const int &shape_I, const double &rate_I,
                                 const int &max_time) {
  // copy conditions
  Rcpp::IntegerVector current_conditions = initial_conditions;

  // exposed and infectious rates --- must be Rcpp vectors
  Rcpp::NumericVector exposed_rates =
      helpers::prob_discrete_erlang(shape_E, rate_E);
  Rcpp::NumericVector infectious_rates =
      helpers::prob_discrete_erlang(shape_I, rate_I);

  // numbers of exposed and infectious blocks
  const int n_exposed_blocks = exposed_rates.size();
  const int n_infectious_blocks = infectious_rates.size();

  // vectors to store the numbers of individuals in each block, past and
  // current
  std::vector<int> exposed_blocks_past(n_exposed_blocks);
  std::vector<int> infectious_blocks_past(n_infectious_blocks);

  // fill vectors with initial values
  // exposed is expected to be the second element -- hardcoded
  // infectious is expected to be the third element -- hardcoded
  // TODO(all): replace all position-based access with name-based access via map
  exposed_blocks_past.back() = initial_conditions[1];
  infectious_blocks_past.back() = initial_conditions[2];

  // vec-of-vecs matrix for data storage --- four columns, 1 per compartment
  Rcpp::IntegerMatrix data_matrix(max_time + 1, 4);

  // assign initial conditions at time = 0
  data_matrix(0, Rcpp::_) = initial_conditions;

  // run the simulation from time 1 to max time (inclusive of max time)
  for (size_t time = 1; time <= max_time; time++) {
    // vectors for current values --- hold zeros
    std::vector<int> exposed_blocks_current(n_exposed_blocks);
    std::vector<int> infectious_blocks_current(n_infectious_blocks);

    // get current probability of exposure
    const double beta_now = beta * static_cast<double>(current_conditions[2]) /
                            static_cast<double>(population_size);
    const double prob_exposure = 1.0 - std::exp(-beta_now);

    // get new exposures, infectious, and recovered
    const int new_exposed = Rcpp::rbinom(
        1.0, static_cast<double>(current_conditions[0]), prob_exposure)[0];
    const int new_infectious = exposed_blocks_past[0];    // first index
    const int new_recovered = infectious_blocks_past[0];  // first index

    // handle movement across exposed blocks
    if (new_exposed > 0) {
      exposed_blocks_current = Rcpp::as<std::vector<int> >(
          helpers::rmultinom_1(new_exposed, exposed_rates, n_exposed_blocks));
    }
    // number of individuals in each block
    exposed_blocks_current[n_exposed_blocks - 1] += 0;
    for (size_t i = 0; i < n_exposed_blocks - 1; i++) {
      exposed_blocks_current[i] += exposed_blocks_past[i + 1];
    }

    // handle movement across infectious blocks
    if (new_infectious > 0) {
      infectious_blocks_current =
          Rcpp::as<std::vector<int> >(helpers::rmultinom_1(
              new_infectious, infectious_rates, n_infectious_blocks));
    }
    // number of individuals in each block
    infectious_blocks_current.back() += 0;
    for (size_t i = 0; i < n_infectious_blocks - 1; i++) {
      infectious_blocks_current[i] += infectious_blocks_past[i + 1];
    }

    // update current conditions
    current_conditions[0] = current_conditions[0] - new_exposed;
    current_conditions[1] = std::accumulate(exposed_blocks_current.begin(),
                                            exposed_blocks_current.end(), 0);
    current_conditions[2] = std::accumulate(infectious_blocks_current.begin(),
                                            infectious_blocks_current.end(), 0);
    current_conditions[3] = current_conditions[3] + new_recovered;

    // log data
    data_matrix(time, Rcpp::_) = current_conditions;

    // swap past vectors vectors
    exposed_blocks_past = exposed_blocks_current;
    infectious_blocks_past = infectious_blocks_current;
  }

  return Rcpp::List::create(Rcpp::Named("x") = data_matrix,
                            Rcpp::Named("time") = Rcpp::seq(0, max_time));
}

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_EBOLA_H_
