// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.
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
#include "epidemic_default.h"  // paths must be relative to this file

#endif  // epidemics_epidemics_H_
