// ZINB_MLE.h

#ifndef ZINB_MLE_H
#define ZINB_MLE_H

#include <Rcpp.h>

Rcpp::List ZINB_MLE(Rcpp::NumericVector X, int iterMax = 15000, double EPSILON = 1e-7, bool verbose = false);

#endif // ZINB_MLE_H
