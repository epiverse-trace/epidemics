// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
#ifndef INST_INCLUDE_SIR_STOCHASTIC_H_
#define INST_INCLUDE_SIR_STOCHASTIC_H_

// clang-format off
#include <vector>
#include <Rcpp.h>
// clang-format on

// Enable C++14 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp14)]]

/// @brief Function to normalise a vector by a maximum value.
/// @param vec A numeric vector.
/// @param N The maximum value of the numeric vector.
/// @return A vector scaled between 0 and 1 by the maximum allowed value.
std::vector<float> normalise_vec(const std::vector<float> &vec,
                                 const float &N) {
  std::vector<float> normalised_vec = vec;
  for (size_t i = 0; i < vec.size(); i++) {
    normalised_vec[i] = vec[i] / N;
  }
  return normalised_vec;
}

/// @brief Function to make a long-format dataframe of epidemic data.
/// @param S A vector of susceptibles at each time point.
/// @param I A vector of infected individuals at each time point.
/// @param R A vector of recovered individuals at each time point.
/// @param time A vector of time points.
/// @return A data frame with columns 'time', 'value', and 'variable'.
Rcpp::DataFrame tidy_epidemic_data(const std::vector<float> &S,
                                   const std::vector<float> &I,
                                   const std::vector<float> &R,
                                   const std::vector<float> &time) {
  // prepare vectors for time etc.
  std::vector<float> value, v_time;

  // prepare variables
  Rcpp::CharacterVector variable =
      Rcpp::rep(Rcpp::CharacterVector({"S", "I", "R"}), time.size());

  // prepare time and values
  for (size_t i = 0; i < time.size(); i++) {
    value.push_back(S[i]);
    value.push_back(I[i]);
    value.push_back(R[i]);
    for (size_t j = 0; j < 3; j++) {
      v_time.push_back(time[i]);
    }
  }

  return Rcpp::DataFrame::create(Rcpp::Named("time") = v_time,
                                 Rcpp::Named("variable") = variable,
                                 Rcpp::Named("value") = value);
}

/// @brief Function for a continuous-time, stochastic SIR epidemic.
/// @param beta Transmission rate \eqn{beta}.
/// @param gamma Recovery rate \eqn{gamma}.
/// @param N Total number of individuals in the population.
/// @param S0 The number of individuals initially susceptible to infection.
/// @param I0 The number of individuals initially infected by the pathogen.
/// @param R0 The number of individuals already recovered from infection and
/// which cannot be infected again.
/// @param tf The final timepoint in the model.
/// @return A data frame of each time point and the proportion of each
/// class at that time point, in long format.
Rcpp::DataFrame sir_stochastic(const float &beta, const float &gamma,
                               const float &N, const float &S0, const float &I0,
                               const float &R0, const float &tf) {
  float t = 0.0;
  float S = S0;
  float I = I0;
  float R = R0;
  // created as numeric although these are integers
  std::vector<float> vec_times, vec_susceptible, vec_infected, vec_recovered;

  // prepare pf, time increment, and r // Clarification required
  float pf1, pf2, pf, dt, r;

  do {
    // record time and SIR classes
    vec_times.push_back(t);
    vec_susceptible.push_back(S);
    vec_infected.push_back(I);
    vec_recovered.push_back(R);

    // calculate probabilities of new infections and recoveries
    pf1 = beta * S * I;
    pf2 = gamma * I;
    pf = pf1 + pf2;

    // draw time to next event
    dt = Rcpp::rexp(1, pf)[0];
    // increment time point
    t += dt;
    // draw from a uniform distribution between 0.0 and 1.0
    // to determine which event (infection or recovery) occurs
    r = Rcpp::runif(1)[0];
    if (r < pf1 / pf) {
      S--;
      I++;
    } else {
      I--;
      R++;
    }
    // simulation ends when there are no more infected individuals
    if (I == 0) {
      break;
    }
  } while (t <= tf && (I > 0));  // simulation runs while there is time, and
                                 // individuals are infected

  // create output data
  return tidy_epidemic_data(normalise_vec(vec_susceptible, N),
                            normalise_vec(vec_infected, N),
                            normalise_vec(vec_recovered, N), vec_times);
}

#endif  // INST_INCLUDE_SIR_STOCHASTIC_H_
