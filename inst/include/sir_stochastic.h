// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef SRC_SIRSTOCHASTIC_H_
#define SRC_SIRSTOCHASTIC_H_

// clang-format off
#include <algorithm>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include <Rcpp.h>
// clang-format on

/// function for a stochastic, continuous time epidemic
Rcpp::DataFrame sir_stochastic(const float &beta, const float &gamma,
                               const float &N, const float &S0, const float &I0,
                               const float &R0, const float &tf) {
  float t = 0.0;
  float S = S0;
  float I = I0;
  float R = R0;
  std::vector<float> ta, Sa, Ia, Ra;
  do {
    ta.push_back(t);
    Sa.push_back(S);
    Ia.push_back(I);
    Ra.push_back(R);
    float pf1 = beta * S * I;
    float pf2 = gamma * I;
    float pf = pf1 + pf2;
    float dt = Rcpp::rexp(1, pf)[0];
    t += dt;
    float r = Rcpp::runif(1)[0];
    if (r < pf1 / pf) {
      S--;
      I++;
    } else {
      I--;
      R++;
    }
    if (I == 0) {
      break;
    }
  } while (t <= tf && (I > 0));
  return Rcpp::DataFrame::create(Rcpp::Named("time") = ta,
                                 Rcpp::Named("S") = Sa, Rcpp::Named("I") = Ia,
                                 Rcpp::Named("R") = Ra);
}

#endif  // SRC_SIRSTOCHASTIC_H_
