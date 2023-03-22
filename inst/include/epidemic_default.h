// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_EPIDEMIC_DEFAULT_H_
#define INST_INCLUDE_EPIDEMIC_DEFAULT_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]

// clang-format off
#include <Rcpp.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/* The rhs of x' = f(x) defined as a struct with an operator */
struct epidemic_default {
  const float beta, alpha, gamma;
  // npi, interv, pop
  epidemic_default(float beta, float alpha, float gamma)
      : beta(beta), alpha(alpha), gamma(gamma) {}

  void operator()(odetools::state_type const& x, odetools::state_type& dxdt,
                  const double t) {
    dxdt[0] = -beta * x[0] * x[2];                    // -beta*S*I
    dxdt[1] = (beta * x[0] * x[2]) - (alpha * x[1]);  // beta*S*I - alpha*E
    dxdt[2] = (alpha * x[1]) - (gamma * x[2]);        // alpha*E - gamma*I
    dxdt[3] = gamma * x[2];                           // gamma*I
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_DEFAULT_H_
