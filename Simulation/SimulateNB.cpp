// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#include "Generate_Samples.h"
#include "NB_MLE.h"
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
        if (num_iterations < max_iter) { 
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
