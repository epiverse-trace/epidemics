// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <epidemics.h>

#include <vector>

//' @title Run an SEIR model with Erlang passage times
//'
//' @description Wrapper function for an SEIR compartmental model with Erlang
//' passage times based on Getz and Dougherty (2017) J. Biological Dynamics,
//' and developed to model the West African ebola virus disease outbreak of 2014
//' .
//'
//' @param initial_conditions An integer vector of the initial conditions
//' corresponding to an SEIR model.
//' @param population_size A single integer for the population size.
//' @param beta The transmission rate \eqn{\beta}.
//' @param shape_E A single integer for the shape parameter of the Erlang
//' distribution of passage times through the exposed compartment.
//' @param rate_E A single double for the rate parameter of the Erlang
//' distribution of passage times through the exposed compartment.
//' @param shape_I A single integer for the shape parameter of the Erlang
//' distribution of passage times through the infectious compartment.
//' @param rate_I A single double for the rate parameter of the Erlang
//' distribution of passage times through the infectious compartment.
//' @param time_end A single integer for the maximum simulation time; this is
//' assumed to be in days.
//' @return A list with two elements, `x`, and integer matrix with as many rows
//' as the number of timesteps (given by `time_end`), and four columns, one for
//' each compartment, susceptible, exposed, infectious, recovered, in that order
//' ; and `times`, a vector of the simulation times, taken to be days.
//' This output is intended to be passed to [output_to_df()] to be converted
//' into a data.frame for further analysis.
//' @keywords internal
// [[Rcpp::export(name=".model_ebola_cpp")]]
Rcpp::List model_ebola_cpp_internal(const Rcpp::IntegerVector &initial_state,
                                       const int &population_size,
                                       const double &beta, const int &shape_E,
                                       const double &rate_E, const int &shape_I,
                                       const double &rate_I,
                                       const int &time_end) {
  return epidemics::epidemic_ebola(initial_state, population_size, beta,
                                   shape_E, rate_E, shape_I, rate_I, time_end);
}
