// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
/*
    Adapted from a function by Hạ Minh Lâm for the Epirecipes Cookbook, to
    replicate a stochastic discrete-time model for ebola virus disease
    developed by Wayne Getz and Eric Dougherty (Getz and Dougherty 2017).
    Source code for the model in Epirecipes is at:
    http://epirecip.es/epicookbook/chapters/erlang/r
    and the model is described in the Journal of Biological Dynamics at:
    https://doi.org/10.1080/17513758.2017.1401677
*/
#ifndef INST_INCLUDE_HELPERS_H_
#define INST_INCLUDE_HELPERS_H_

#include <Rcpp.h>

#include <cmath>
#include <vector>

namespace helpers {
inline Rcpp::NumericVector prob_discrete_erlang(const int &shape,
                                                const double &rate) {
  // prepare factorials vector, note loop starts at 1
  std::vector<double> factorials(shape + 1, 0.0);
  factorials[0] = 1.0;
  for (int i = 1; i <= shape; i++) factorials[i] = i * factorials[i - 1];

  double cumulative_prob = 0.0;
  int n_bin = 0;
  std::vector<double> vec_cumulative_probs{1.0};  // note vec of similar name

  while (cumulative_prob <= 0.99) {
    n_bin++;  // increment bin number
    double val = 0.0;

    for (int j = 0; j < shape; j++) {
      val += std::exp(-n_bin * rate) * std::pow((n_bin * rate), j) /
             factorials[j + 1];
    }
    vec_cumulative_probs.push_back(val);
    cumulative_prob = 1.0 - val;
  }

  // probabilities at each integer
  std::vector<double> density_prob(vec_cumulative_probs.size() - 1);
  for (size_t i = 0; i < density_prob.size(); i++) {
    density_prob[i] = vec_cumulative_probs[i] - vec_cumulative_probs[i + 1];
  }

  return Rcpp::wrap(density_prob);
}
}  // namespace helpers

#endif  // INST_INCLUDE_HELPERS_H_
