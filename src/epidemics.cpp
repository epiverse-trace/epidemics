// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <sir_stochastic.h>

// [[Rcpp::interfaces(r, cpp)]]

//' @title Stochastic, continuous time SIR epidemic simulation
//' @description An initial implementation of an SIR epidemic simulation.
//' @param parameters A named list of parameters to the simulation.
//'
// [[Rcpp::export()]]
Rcpp::List run_sir_stochastic(const Rcpp::List &parameters) {
  return epidemics::sir_stochastic(
      parameters["beta"], parameters["gamma"], parameters["N"],
      parameters["S0"], parameters["I0"], parameters["R0"], parameters["tf"]);
}
