// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <epidemics.h>

//' @title Compute the discrete probability of the truncated Erlang distribution
//' @name prob_discrete_erlang
//' @rdname prob_discrete_erlang
//'
// [[Rcpp::export]]
Rcpp::NumericVector prob_discrete_erlang_cpp(const int &shape,
                                             const double &rate) {
  return helpers::prob_discrete_erlang(shape, rate);
}
