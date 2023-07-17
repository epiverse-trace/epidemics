// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <helpers.h>

//' @title Compute the discrete probability of the truncated Erlang distribution
//'
//' @description A helper function that gives the probability of discrete values
//' from an Erlang distribution with a given shape and rate. The number of
//' values returned correspond to the number of discrete values over which the
//' cumulative probability reaches 0.99.
//'
//' @param shape A single integer-like number for the shape of the Erlang
//' distribution.
//' @param rate A single number for the rate of the Erlang distribution.
//' @return A vector of variable length giving the probability of each integer
//' value for a cumulative probability of 0.99.
//' @export
// [[Rcpp::export(name=".prob_discrete_erlang")]]
Rcpp::NumericVector prob_discrete_erlang(const int &shape, const double &rate) {
  return helpers::prob_discrete_erlang(shape, rate);
}
