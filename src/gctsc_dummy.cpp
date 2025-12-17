#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "gauss_utils.h"
#include "sampling_utils.h"
#include "lnNpr.h"

using namespace Rcpp;

// dummy function so CRAN recognizes the cpp file
extern "C" void gctsc_dummy() {}
