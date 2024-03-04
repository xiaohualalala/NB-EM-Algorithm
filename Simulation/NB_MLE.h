// NB_MLE.h

#ifndef NB_MLE_H
#define NB_MLE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>

Rcpp::List Infor_Matrix(Rcpp::NumericVector X, Rcpp::List Param);
Rcpp::List NB_MLE(Rcpp::NumericVector X, int iterMax = 15000, double EPSILON = 1e-7, bool verbose = false, bool se = false, bool rate = false);

#endif // NB_MLE_H
