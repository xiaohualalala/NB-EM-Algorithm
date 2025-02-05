// MixNB_MLE.h

#ifndef MixNB_MLE_H
#define MixNB_MLE_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List mixNB_MLE(Rcpp::NumericVector X, int nbS, Rcpp::NumericVector avg = Rcpp::NumericVector(), Rcpp::NumericVector var = Rcpp::NumericVector(), int iterMax = 15000, double EPSILON = 1e-7, bool verbose = false);

#endif // MixNB_MLE_H