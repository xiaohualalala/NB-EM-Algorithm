// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#include "../Samples/Generate_Samples.h"
#include "../../MLE/C/ZINB_MLE.h"
/*-------------------------------------------------------*
 * Function: SimulateZINB
 * Purpose: Simulate and analyze samples from a zero-inflated negative binomial distribution.
 * Input:
 *   prob (double) - Probability parameter for the zero-inflated negative binomial distribution.
 *   size (double) - Size parameter for the zero-inflated negative binomial distribution.
 *   zero (double) - Zero inflation parameter for the zero-inflated negative binomial distribution.
 *   num_samples (int) - Number of samples to simulate.
 *   sample_size (int) - Size of each individual sample.
 * Output:
 *   List - Contains iterations, estimated probabilities, estimated sizes and estimated zeros.
 * Description:
 *   This function simulates samples from a zero-inflated negative binomial distribution,
 *   applies the ZINB_MLE function to each sample, and aggregates the results.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List SimulateZINB(double prob, double size, double zero, int num_samples, int sample_size, int max_iter = 15000)
{
    List ZINB = SamplesZINB(prob, size, zero, num_samples, sample_size);
    List data = ZINB["samples"];
    List labels= ZINB["labels"];

    std::vector<int> iterations(num_samples);
    std::vector<double> prob_estimates(num_samples);
    std::vector<double> size_estimates(num_samples);
    std::vector<double> zero_estimates(num_samples);

    int index = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector X = data[i];
        List results = ZINB_MLE(X);

        List stats = results["Statistics"];
        List parameter = results["Parameter_Estimates"];

        int num_iterations = stats["Iterations"];
        if (num_iterations < max_iter) { 
            iterations[index] = num_iterations;
            prob_estimates[index] = parameter["eProb"];
            size_estimates[index] = parameter["eSize"];
            zero_estimates[index] = parameter["eZero"];
            index++;
        }
    }

    iterations.resize(index); 
    prob_estimates.resize(index);
    size_estimates.resize(index);
    zero_estimates.resize(index);

    return List::create(
        Named("iterations") = iterations,
        Named("prob_estimates") = prob_estimates,
        Named("size_estimates") = size_estimates,
        Named("zero_estimates") = zero_estimates,
        Named("labels") = labels);
}

/*-------------------------------------------------------*
 * Function: SimulateZINBs
 * Purpose: Perform simulations across a range of probabilities, sizes and zeros,
 *          and summarize the results for zero-inflated negative binomial maximum likelihood estimation.
 * Input:
 *   prob_values (NumericVector) - Vector of probabilities to simulate.
 *   size_values (NumericVector) - Vector of sizes to simulate.
 *   zero_values (NumericVector) - Vector of zeros to simulate.
 *   num_samples (int) - Number of samples to generate for each combination of probability and size.
 *   sample_size (int) - Size of each sample to generate.
 * Output:
 *   DataFrame - A data frame containing the summary of results for each combination of probability, size and zero.
 * Description:
 *   This function iterates over combinations of probabilities, sizes and zeros,
 *   performs simulations, applies ZINB_MLE, and summarizes the results in a data frame.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
DataFrame SimulateZINBs(NumericVector prob_values, NumericVector size_values, NumericVector zero_values, int num_samples, int sample_size)
{
    int num_prob = prob_values.size();
    int num_size = size_values.size();
    int num_zero = zero_values.size();

    std::vector<double> avg_prob, sd_prob, avg_size, sd_size, avg_t, sd_t, avg_zero, sd_zero;
    std::vector<double> probs, sizes, zeros;

    for (int i = 0; i < num_prob; ++i)
    {
        for (int j = 0; j < num_size; ++j)
        {
            for (int k = 0; k < num_zero; ++k)
            {
                List result = SimulateZINB(prob_values[i], size_values[j], zero_values[k], num_samples, sample_size);

                NumericVector prob_estimates = result["prob_estimates"];
                NumericVector size_estimates = result["size_estimates"];
                NumericVector zero_estimates = result["zero_estimates"];
                NumericVector iterations = result["iterations"];

                probs.push_back(prob_values[i]);
                sizes.push_back(size_values[j]);
                zeros.push_back(zero_values[k]);
                avg_prob.push_back(mean(prob_estimates));
                sd_prob.push_back(sd(prob_estimates));
                avg_size.push_back(mean(size_estimates));
                sd_size.push_back(sd(size_estimates));
                avg_zero.push_back(mean(zero_estimates));
                sd_zero.push_back(sd(zero_estimates));
                avg_t.push_back(mean(iterations));
                sd_t.push_back(sd(iterations));
            }
        }
    }

    return DataFrame::create(
        Named("zero") = zeros,
        Named("prob") = probs,
        Named("size") = sizes,
        Named("avg_zero") = avg_zero,
        Named("sd_zero") = sd_zero,
        Named("avg_prob") = avg_prob,
        Named("sd_prob") = sd_prob,
        Named("avg_size") = avg_size,
        Named("sd_size") = sd_size,
        Named("avg_t") = avg_t,
        Named("sd_t") = sd_t);
}