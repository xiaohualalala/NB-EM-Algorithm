// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#include "../Samples/Generate_Samples.h"
#include "../../MLE/C/NB_MLE.h"
/*-------------------------------------------------------*
 * Function: SimulateNB
 * Purpose: Simulate and analyze samples from a negative binomial distribution.
 * Input:
 *   prob (double) - Probability parameter for the negative binomial distribution.
 *   size (double) - Size parameter for the negative binomial distribution.
 *   num_samples (int) - Number of samples to simulate.
 *   sample_size (int) - Size of each individual sample.
 * Output:
 *   List - Contains iterations, estimated probabilities, and estimated sizes.
 * Description:
 *   This function simulates samples from a negative binomial distribution,
 *   applies the NB_MLE function to each sample, and aggregates the results.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List SimulateNB(double prob, double size, int num_samples, int sample_size, int max_iter = 15000)
{
    List data = SamplesNB(prob, size, num_samples, sample_size);

    std::vector<int> iterations(num_samples);
    std::vector<double> prob_estimates(num_samples);
    std::vector<double> size_estimates(num_samples);

    int index = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector X = data[i];
        List results = NB_MLE(X);

        List stats = results["Statistics"];
        List parameter = results["Parameter_Estimates"];

        int num_iterations = stats["Iterations"];
        if (num_iterations < max_iter)
        {
            iterations[index] = num_iterations;
            prob_estimates[index] = parameter["eProb"];
            size_estimates[index] = parameter["eSize"];
            index++;
        }
    }

    iterations.resize(index);
    prob_estimates.resize(index);
    size_estimates.resize(index);

    return List::create(
        Named("iterations") = wrap(iterations),
        Named("prob_estimates") = wrap(prob_estimates),
        Named("size_estimates") = wrap(size_estimates));
}

/*-------------------------------------------------------*
 * Function: SimulateNBs
 * Purpose: Perform simulations across a range of probabilities and sizes,
 *          and summarize the results for negative binomial maximum likelihood estimation.
 * Input:
 *   prob_values (NumericVector) - Vector of probabilities to simulate.
 *   size_values (NumericVector) - Vector of sizes to simulate.
 *   num_samples (int) - Number of samples to generate for each combination of probability and size.
 *   sample_size (int) - Size of each sample to generate.
 * Output:
 *   DataFrame - A data frame containing the summary of results for each combination of probability and size.
 * Description:
 *   This function iterates over combinations of probabilities and sizes,
 *   performs simulations, applies NB_MLE, and summarizes the results in a data frame.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
DataFrame SimulateNBs(NumericVector prob_values, NumericVector size_values, int num_samples, int sample_size)
{
    int num_prob = prob_values.size();
    int num_size = size_values.size();

    std::vector<double> avg_prob, sd_prob, avg_size, sd_size, avg_t, sd_t;
    std::vector<double> probs, sizes;

    for (int i = 0; i < num_prob; ++i)
    {
        for (int j = 0; j < num_size; ++j)
        {
            List result = SimulateNB(prob_values[i], size_values[j], num_samples, sample_size);

            NumericVector prob_estimates = result["prob_estimates"];
            NumericVector size_estimates = result["size_estimates"];
            NumericVector iterations = result["iterations"];

            probs.push_back(prob_values[i]);
            sizes.push_back(size_values[j]);
            avg_prob.push_back(mean(prob_estimates));
            sd_prob.push_back(sd(prob_estimates));
            avg_size.push_back(mean(size_estimates));
            sd_size.push_back(sd(size_estimates));
            avg_t.push_back(mean(iterations));
            sd_t.push_back(sd(iterations));
        }
    }

    return DataFrame::create(
        Named("prob") = probs,
        Named("size") = sizes,
        Named("avg_prob") = avg_prob,
        Named("sd_prob") = sd_prob,
        Named("avg_size") = avg_size,
        Named("sd_size") = sd_size,
        Named("avg_t") = avg_t,
        Named("sd_t") = sd_t);
}

// [[Rcpp::export]]
List SimulateRate(double prob, double size, int num_samples, int sample_size, int max_iter = 15000)
{
    List data = SamplesNB(prob, size, num_samples, sample_size);

    std::vector<double> conv_rates(num_samples);

    int index = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector X = data[i];
        List results = NB_MLE(X);

        List stats = results["Statistics"];
        List parameter = results["Parameter_Estimates"];
        double eProb = parameter["eProb"];
        double eSize = parameter["eSize"];
        int iter = stats["Iterations"];

        if (iter < max_iter)
        {
            double alpha = -1.0 / log(eProb);
            double Sigma = 0.0;

            for (int j = 0; j < sample_size; j++)
            {
                Sigma += R::trigamma(eSize) - R::trigamma(eSize + X[j]);
            }

            double a1 = alpha * alpha * (1 - eProb) * (1 - eProb) / (alpha * (1 - eProb) - eProb);
            double b1 = (alpha * alpha * (1 - eProb) * (1 - eProb) - eProb) / (alpha * (1 - eProb) - eProb);
            double j = 1 - alpha * eSize * Sigma / (1.0 * sample_size);
            double b = (a1 + 1) * j + b1;
            double c = b1 * j;
            double r = (b + sqrt(b * b - 4 * c)) / 2;
            conv_rates[index] = r;
            
            index++;
        }
    }
    conv_rates.resize(index);

    return List::create(
        Named("conv_rates") = wrap(conv_rates));
}


// [[Rcpp::export]]
DataFrame SimulateRates(NumericVector prob_values, NumericVector size_values, int num_samples, int sample_size)
{
    int num_prob = prob_values.size();
    int num_size = size_values.size();

    std::vector<double> avg_rate, sd_rate;
    std::vector<double> probs, sizes;

    for (int i = 0; i < num_prob; ++i)
    {
        for (int j = 0; j < num_size; ++j)
        {
            List result = SimulateRate(prob_values[i], size_values[j], num_samples, sample_size);

            NumericVector rate = result["conv_rates"];

            probs.push_back(prob_values[i]);
            sizes.push_back(size_values[j]);
            avg_rate.push_back(mean(rate));
            sd_rate.push_back(sd(rate));
        }
    }

    return DataFrame::create(
        Named("prob") = probs,
        Named("size") = sizes,
        Named("avg_rate") = avg_rate,
        Named("sd_rate") = sd_rate);
}


// [[Rcpp::export]]
List SimulateError(double prob, double size, int num_samples, int sample_size, int max_iter = 15000)
{
    List data = SamplesNB(prob, size, num_samples, sample_size);

    std::vector<double> error_1(num_samples);
    std::vector<double> error_2(num_samples);

    int index = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector X = data[i];
        List results = NB_MLE(X);

        List stats = results["Statistics"];
        List parameter = results["Parameter_Estimates"];
        double eProb = parameter["eProb"];
        double eSize = parameter["eSize"];
        int iter = stats["Iterations"];

        if (iter < max_iter)
        {
            double alpha = -1.0 / log(eProb);
            double Sigma = 0.0;

            for (int j = 0; j < sample_size; j++)
            {
                Sigma += R::trigamma(eSize) - R::trigamma(eSize + X[j]);
            }

            arma::mat obeI(2, 2);

            obeI(0, 0) = eSize / (pow(eProb,2) * (1 - eProb)) * (sample_size * (1 - 2 * alpha * (1-eProb)) + pow(alpha,2) * eSize * (1-eProb) * Sigma);

            obeI(0, 1) = - alpha / eProb * (sample_size - alpha * eSize * Sigma);

            obeI(1, 0) = obeI(0, 1);

            obeI(1, 1) = pow(alpha, 2) * Sigma;

            arma::mat obeIinv(2, 2);
            bool invertible = arma::inv(obeIinv, obeI);

            error_1[index] = sqrt(obeIinv(0,0));
            error_2[index] = sqrt(obeIinv(1,1));
            index++;
        }
    }
    error_1.resize(index);
    error_2.resize(index);

    return List::create(
        Named("error_1") = wrap(error_1),
        Named("error_2") = wrap(error_2));
}

// [[Rcpp::export]]
DataFrame SimulateErrors(NumericVector prob_values, NumericVector size_values, int num_samples, int sample_size)
{
    int num_prob = prob_values.size();
    int num_size = size_values.size();

    std::vector<double> avg_1, avg_2, sd_1, sd_2;
    std::vector<double> probs, sizes;

    for (int i = 0; i < num_prob; ++i)
    {
        for (int j = 0; j < num_size; ++j)
        {
            List result = SimulateError(prob_values[i], size_values[j], num_samples, sample_size);

            NumericVector error_1 = result["error_1"];
            NumericVector error_2 = result["error_2"];

            probs.push_back(prob_values[i]);
            sizes.push_back(size_values[j]);
            avg_1.push_back(mean(error_1));
            avg_2.push_back(mean(error_2));
            sd_1.push_back(sd(error_1));
            sd_2.push_back(sd(error_2));
        }
    }

    return DataFrame::create(
        Named("prob") = probs,
        Named("size") = sizes,
        Named("avg_error_1") = avg_1,
        Named("sd_error_1") = sd_1,
        Named("avg_error_2") = avg_2,
        Named("sd_error_2") = sd_2);
}