// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <sir_stochastic.h>

// [[Rcpp::interfaces(r, cpp)]]

//' @title Stochastic, continuous time SIR epidemic simulation
//' @description An initial implementation of an SIR epidemic simulation.
//' @param parameters A named list of parameters to the simulation.
//'
// [[Rcpp::export(name = ".sir_stochastic")]]
Rcpp::List sir_stochastic(const Rcpp::List &parameters) {
  return sir_stochastic_(parameters["beta"], parameters["gamma"],
                         parameters["N"], parameters["S0"], parameters["I0"],
                         parameters["R0"], parameters["tf"]);
}
