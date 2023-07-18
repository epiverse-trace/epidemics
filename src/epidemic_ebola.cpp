// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

//' @title Run an SEIR model with Erlang passage times
//'
//' @param population An object of the `population` class, which holds a
//' population contact matrix, a demography vector, and the initial conditions
//' of each demographic group. See [population()].
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
//' @return An integer matrix with as many rows as the number of timesteps
//' (given by `time_end`), and four columns, one for each compartment, SEIR,
//' in that order.
//' @export
// [[Rcpp::export(name=".epidemic_ebola_cpp")]]
Rcpp::IntegerMatrix epidemic_ebola_cpp(const Rcpp::List population,
                                       const double &beta, const int &shape_E,
                                       const double &rate_E, const int &shape_I,
                                       const double &rate_I,
                                       const int &time_end) {
  return epidemics::epidemic_ebola(beta, shape_E, rate_E, shape_I, rate_I,
                                   time_end, population);
}
