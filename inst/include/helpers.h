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

/// @brief Get the discrete probabilities of an Erlang distribution
/// @param shape An integer number for the shape parameter of the Erlang
/// distribution.
/// @param rate A single number for the rate parameter of the Erlang
/// distribution.
/// @return An Rcpp double (numeric) vector of variable length. The length is
/// determined by the value at which the cumulative probability is >= 0.99.
inline Rcpp::NumericVector prob_discrete_erlang(const int &shape,
                                                const double &rate) {
  // prepare factorials vector, note loop starts at 1
  // different length from R implementation as last value seemingly never used
  std::vector<double> factorials(shape, 0.0);
  factorials[0] = 1.0;
  for (int i = 1; i < shape; i++) factorials[i] = i * factorials[i - 1];

  // initial values for while loop
  double cumulative_prob = 0.0;
  int n_bin = 1;
  std::vector<double> vec_cumulative_probs{1.0};  // note vec of similar name

  // while loop filling probability values
  while (cumulative_prob <= 0.99) {  // hardoced to 0.99
    double val = 0.0;

    for (int j = 0; j < shape; j++) {
      val += (std::exp(static_cast<double>(-n_bin) * rate) *
              std::pow((static_cast<double>(n_bin) * rate),
                       static_cast<double>(j)) /
              static_cast<double>(factorials[j]));
      // different from R due to 0 indexing
    }
    cumulative_prob = 1.0 - val;
    vec_cumulative_probs.push_back(val);
    n_bin++;  // increment bin number
  }

  // probabilities at each integer
  std::vector<double> density_prob(vec_cumulative_probs.size() - 1);
  for (size_t i = 0; i < density_prob.size(); i++) {
    density_prob[i] = (vec_cumulative_probs[i] - vec_cumulative_probs[i + 1]) /
                      cumulative_prob;
  }

  return Rcpp::wrap(density_prob);
}

/// Replicating the stats::rmultinom() function vectorised across a range of
/// probabilities, taken from
/// https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/

/// @brief Draw from a multinomial distribution.
/// @param size A single integer number for the number of multinomial outcomes
/// to sample for each data set.
/// @param probs An Rcpp numeric vector of probabilities.
/// @param N A single integer for the number of simulated data sets to produce.
/// @return An integer vector of size N, giving draws from multinomial outcomes.
inline Rcpp::IntegerVector rmultinom_1(const int &size,
                                       Rcpp::NumericVector &probs,  // NOLINT
                                       const int &N) {
  Rcpp::IntegerVector outcome(N);
  R::rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

}  // namespace helpers

#endif  // INST_INCLUDE_HELPERS_H_
