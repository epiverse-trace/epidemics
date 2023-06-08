// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
/*
    Adapted from the Vacamole model developed by Kylie Ainslie et al. at the
    National Institute for Public Health and the Environment (RIVM),
    The Netherlands, during the Covid-19 pandemic.
    Source code for the original Vacamole model was developed at:
    https://github.com/kylieainslie/vacamole
    and the model is described in Eurosurveillance at:
    https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090
*/
#ifndef INST_INCLUDE_EPIDEMIC_VACAMOLE_H_
#define INST_INCLUDE_EPIDEMIC_VACAMOLE_H_

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/numeric/odeint.hpp>
#include "ode_tools.h"
#include "intervention.h"
#include "vaccination.h"
#include "population.h"
// clang-format on

// add to namespace epidemics
namespace epidemics {

/* The rhs of x' = f(x) defined as a struct with an operator */
struct epidemic_vacamole {
  const double beta, beta_v, alpha, omega, omega_v, eta, eta_v, gamma;
  const Eigen::MatrixXd nu;
  const Rcpp::List population, intervention, vaccination;
  const Eigen::MatrixXd contact_matrix;
  // npi, interv, pop
  epidemic_vacamole(const double beta, const double beta_v, const double alpha,
                    const double omega, const double omega_v, const double eta,
                    const double eta_v, const double gamma,
                    const Rcpp::List population, const Rcpp::List intervention,
                    const Rcpp::List vaccination)
      : beta(beta),
        beta_v(beta_v),
        alpha(alpha),
        omega(omega),
        omega_v(omega_v),
        eta(eta),
        eta_v(eta_v),
        gamma(gamma),
        nu(vaccination::get_nu(vaccination)),
        population(population),
        intervention(intervention),
        vaccination(vaccination),
        contact_matrix(population::get_contact_matrix(population)) {}

  void operator()(const odetools::state_type& x,
                  odetools::state_type& dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of x
    dxdt.resize(x.rows(), x.cols());

    // modify contact matrix if time is within intervention timespan
    Eigen::MatrixXd cm =
        intervention::intervention_on_cm(t, contact_matrix, intervention);

    // get current vaccination rate
    Eigen::MatrixXd current_nu = vaccination::current_nu(nu, vaccination, t);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations
    // columns are as follows
    // 0| 1| 2|3| 4|5| 6|7| 8|9|10
    // S|V1|V2|E|EV|I|IV|H|HV|D|R

    // compartmental transitions without accounting for contacts
    // Susceptible to exposed
    Eigen::ArrayXd sToE =
        beta * x.col(0).array() * (cm * (x.col(5) + x.col(6))).array();

    // Susceptible to vaccinated with one dose
    Eigen::ArrayXd sToV1 = current_nu.col(0).array() * x.col(0).array();
    // Vaccinated one dose to vaccinated with two doses
    Eigen::ArrayXd v1ToV2 = current_nu.col(1).array() * x.col(1).array();

    // Vaccinated one dose to exposed - same as susceptible to exposed
    Eigen::ArrayXd v1ToE =
        beta * x.col(1).array() * (cm * (x.col(5) + x.col(6))).array();
    // Vaccinated two doses to exposed - uses different beta
    Eigen::ArrayXd v2ToEv =
        beta_v * x.col(2).array() * (cm * (x.col(5) + x.col(6))).array();

    // Exposed unvaccinated or not protected to infectious
    Eigen::ArrayXd eToI = alpha * x.col(3).array();
    // Exposed vaccinated to infectious
    Eigen::ArrayXd evToIv = alpha * x.col(4).array();

    // Infectious to hospitalised
    Eigen::ArrayXd iToH = eta * x.col(5).array();
    // Vaccinated infectious to hospitalised
    Eigen::ArrayXd ivToHv = eta_v * x.col(6).array();

    // Infectious to dead
    Eigen::ArrayXd iToD = omega * x.col(5).array();
    // Infectious vaccinated to dead
    Eigen::ArrayXd ivToD = omega_v * x.col(6).array();

    // Hospitalised to dead
    Eigen::ArrayXd hToD = omega * x.col(7).array();
    // Hospitalised vaccinated to dead
    Eigen::ArrayXd hvToD = omega_v * x.col(8).array();

    // Infectious to recovered
    Eigen::ArrayXd iToR = gamma * x.col(5).array();
    // Infectious vaccinated to recovered
    Eigen::ArrayXd ivToR = gamma * x.col(6).array();

    // Hospitalised to recovered
    Eigen::ArrayXd hToR = gamma * x.col(7).array();
    // Hospitalised vaccinated to recovered
    Eigen::ArrayXd hvToR = gamma * x.col(8).array();

    // compartmental changes accounting for contacts
    dxdt.col(0) = -sToE - sToV1;           // -β*S*contacts*(I+Iv) - ν1*S
    dxdt.col(1) = sToV1 - v1ToE - v1ToV2;  // ν1*S -β*V1*contacts*(I+Iv) - ν2*V
    dxdt.col(2) = v1ToV2 - v2ToEv;         // ν2*V - βv*V2*contacts*(I+Iv)
    dxdt.col(3) = sToE + v1ToE - eToI;     // β*(S+V1)*contacts*(I+Iv) - αE
    dxdt.col(4) = v2ToEv - evToIv;         // βv*V2*contacts*(I+Iv) - αE
    dxdt.col(5) = eToI - iToH - iToR - iToD;        // αE - I(γ + η + ω)
    dxdt.col(6) = evToIv - ivToHv - ivToR - ivToD;  // αEv - Iv(γ + η_v + ω_v)
    dxdt.col(7) = iToH - hToR - hToD;               // ηI - γH - ωH
    dxdt.col(8) = ivToHv - hvToR - hvToD;           // ηIv - γHv - ωHv
    dxdt.col(9) = iToD + ivToD + hToD + hvToD;      // all deaths
    dxdt.col(10) = iToR + ivToR + hToR + hvToR;     // all recoveries
  }
};

}  // namespace epidemics

#endif  // INST_INCLUDE_EPIDEMIC_VACAMOLE_H_
