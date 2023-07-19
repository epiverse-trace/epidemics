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

inline Rcpp::List epidemic_ebola(const double &beta, const int &shape_E,
                                 const double &rate_E, const int &shape_I,
                                 const double &rate_I, const int &max_time,
                                 const Rcpp::List population) {
  // get population size and initial conditions
  const int population_size = population::get_population_size(population);
  Rcpp::NumericVector initial_conditions =
      Rcpp::NumericVector(population::get_initial_conditions(population));

  // copy conditions
  Rcpp::NumericVector current_conditions = initial_conditions;

  // exposed and infectious rates
  Rcpp::NumericVector exposed_rates =
      helpers::prob_discrete_erlang(shape_E, rate_E);
  Rcpp::NumericVector infectious_rates =
      helpers::prob_discrete_erlang(shape_I, rate_I);

  // numbers of exposed and infectious blocks
  const int n_exposed_blocks = exposed_rates.length();
  const int n_infectious_blocks = infectious_rates.length();

  // vectors to store the numbers of individuals in each block, past and
  // current
  Rcpp::IntegerVector exposed_blocks_past(n_exposed_blocks);
  Rcpp::IntegerVector infectious_blocks_past(n_infectious_blocks);

  // fill vectors with initial values
  // exposed is expected to be the second element -- hardcoded
  // infectious is expeted to be the third element -- hardcoded
  // TODO(all): replace all position-based access with name-based access
  exposed_blocks_past(n_exposed_blocks - 1) = initial_conditions[1];
  infectious_blocks_past(n_infectious_blocks - 1) = initial_conditions[2];

  // matrix for data storage --- four columns for each compartment
  Rcpp::IntegerMatrix data_matrix(max_time + 1, 4L);
  // assign initial conditions
  data_matrix(0, Rcpp::_) = initial_conditions;

  // run the simulation from 1 to max time
  for (size_t time = 1; time <= max_time; time++) {
    // vectors for current values --- hold zeros
    Rcpp::IntegerVector exposed_blocks_current(n_exposed_blocks);
    Rcpp::IntegerVector infectious_blocks_current(n_infectious_blocks);

    // get current probability of exposure
    const double beta_now = beta * static_cast<double>(current_conditions[2]) /
                            static_cast<double>(population_size);
    const double prob_exposure = 1.0 - std::exp(-beta_now);

    // get new exposures, infectious, and recovered
    const int new_exposed = Rcpp::rbinom(
        1.0, static_cast<double>(current_conditions[0]), prob_exposure)[0];
    const int new_infectious = exposed_blocks_past(0);    // first index
    const int new_recovered = infectious_blocks_past(0);  // first index

    // handle movement across exposed blocks
    if (new_exposed > 0) {
      exposed_blocks_current = Rcpp::IntegerVector(Rcpp::transpose(
          helpers::rmultinom_vectorised(1, new_exposed, exposed_rates)));
    }
    // number of individuals in each block
    Rcpp::IntegerVector vals_exposed(n_exposed_blocks);
    vals_exposed[n_exposed_blocks - 1] = 0;
    for (size_t i = 0; i < n_exposed_blocks - 1; i++) {
      vals_exposed[i] = exposed_blocks_past[i + 1];
    }
    exposed_blocks_current += vals_exposed;

    // handle movement across infectious blocks
    if (new_infectious > 0) {
      infectious_blocks_current = Rcpp::IntegerVector(Rcpp::transpose(
          helpers::rmultinom_vectorised(1, new_infectious, infectious_rates)));
    }
    // number of individuals in each block
    Rcpp::IntegerVector vals_infectious(n_infectious_blocks);
    vals_infectious[n_infectious_blocks - 1] = 0;
    for (size_t i = 0; i < n_infectious_blocks - 1; i++) {
      vals_infectious[i] = infectious_blocks_past[i + 1];
    }
    infectious_blocks_current += vals_infectious;

    // log data
    data_matrix(time, 0) = current_conditions[0] - new_exposed;
    data_matrix(time, 1) = Rcpp::sum(exposed_blocks_current);
    data_matrix(time, 2) = Rcpp::sum(infectious_blocks_current);
    data_matrix(time, 3) = current_conditions[3] + new_recovered;

    // update current conditions
    current_conditions[0] = current_conditions[0] - new_exposed;
    current_conditions[1] = Rcpp::sum(exposed_blocks_current);
    current_conditions[2] = Rcpp::sum(infectious_blocks_current);
    current_conditions[3] = current_conditions[3] + new_recovered;

    // swap past vectors vectors
    exposed_blocks_past = exposed_blocks_current;
    infectious_blocks_past = infectious_blocks_current;
  }

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::transpose(data_matrix),
                            Rcpp::Named("time") = Rcpp::seq(0, max_time));
}

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_EBOLA_H_
