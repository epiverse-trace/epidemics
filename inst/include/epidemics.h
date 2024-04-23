// Copyright 2024 'epidemics' authors. See repository licence in LICENSE.md.
// Written manually to allow header export
// copied from https://github.com/r-pkg-examples/rcpp-shared-cpp-functions

#ifndef epidemics_epidemics_H_
#define epidemics_epidemics_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]

// Include the Rcpp header
#include <Rcpp.h>

// paths relative to this file
#include "ode_tools.h"
#include "population.h"
#include "intervention.h"
#include "vaccination.h"
#include "time_dependence.h"
#include "model_default.h"
#include "model_vacamole.h"
#include "model_diphtheria.h"

#endif  // epidemics_epidemics_H_
