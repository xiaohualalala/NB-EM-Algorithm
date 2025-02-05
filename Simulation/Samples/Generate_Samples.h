// Generate_Samples.h

#ifndef GENERATE_SAMPLES_H
#define GENERATE_SAMPLES_H

#include <Rcpp.h>

Rcpp::List SamplesNB(double prob, double size, int num_samples, int sample_size);
Rcpp::List SamplesZINB(double prob, double size, double zero, int num_samples, int sample_size);
Rcpp::List SamplesMixNB(Rcpp::NumericVector probs, Rcpp::NumericVector sizes, Rcpp::NumericVector weights, int num_samples, int sample_size);
Rcpp::List SamplesMixZINB(Rcpp::NumericVector probs, Rcpp::NumericVector sizes, Rcpp::NumericVector zero_probs, Rcpp::NumericVector weights, int num_samples, int sample_size);

#endif // GENERATE_SAMPLES_H
