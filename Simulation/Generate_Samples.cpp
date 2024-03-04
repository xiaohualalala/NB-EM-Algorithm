#include <Rcpp.h>
using namespace Rcpp;

/*-------------------------------------------------------*
 * Function: SamplesNB
 * Purpose: Generate simulation data based on a negative binomial distribution.
 * Input:
 *   prob (double) - The probability of success in each trial.
 *   size (double) - The dispersion parameter (size) for the negative binomial distribution.
 *   num_samples (int) - The number of samples to generate.
 *   sample_size (int) - The size of each sample.
 * Output:
 *   List - A list of numeric vectors, each representing a sample.
 * Description:
 *   This function generates a list of samples from a negative binomial distribution.
 *   The 'prob' and 'size' parameters are used to define the characteristics of the distribution.
 *   Each sample is generated until its variance is greater than its mean,
 *   ensuring a distribution consistent with the negative binomial characteristics.
 * Examples:
 *   samples <- SamplesNB(0.5, 1, 10, 100)
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List SamplesNB(double prob, double size, int num_samples, int sample_size)
{
    // Initialize a list to store the samples
    List samples(num_samples);

    // Generate each sample
    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector sample(sample_size);
        while (true)
        {
            // Generate a sample using the negative binomial distribution
            for (int j = 0; j < sample_size; ++j)
            {
                sample[j] = R::rnbinom(size, prob);
            }
            // Generate a sample using the negative binomial distribution

            // Ensure the sample meets the negative binomial distribution's property
            if (var(sample) > mean(sample))
            {
                break;
            }
        }
        // Store the generated sample in the list
        samples[i] = sample;
    }

    return samples;
}

/*-------------------------------------------------------*
 * Function: SamplesZINB
 * Purpose: Generate samples from a zero-inflated negative binomial (ZINB) distribution.
 * Input:
 *   prob (double) - Probability of success in each trial for the negative binomial part.
 *   size (double) - Dispersion parameter (size) for the negative binomial part.
 *   zero (double) - Probability of generating a zero (inflation part).
 *   num_samples (int) - Number of samples to generate.
 *   sample_size (int) - Size of each sample.
 * Output:
 *   List - A list of numeric vectors, each representing a ZINB sample.
 * Description:
 *   This function generates a list of samples from a ZINB distribution.
 *   For each observation in a sample, there is a zero_prob chance of being zero
 *   (zero-inflated part), otherwise, it is generated from a negative binomial
 *   distribution defined by prob and size.
 *   The function is useful for simulating count data with excess zeros.
 * Examples:
 *   // Example usage in R
 *   samples <- SamplesZINB(0.5, 1.0, 0.2, 10, 100)
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List SamplesZINB(double prob, double size, double zero, int num_samples, int sample_size)
{
    List samples(num_samples);

    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector sample(sample_size);
        for (int j = 0; j < sample_size; ++j)
        {
            sample[j] = (R::runif(0, 1) < zero) ? 0 : R::rnbinom(size, prob);
        }
        samples[i] = sample;
    }

    return samples;
}

/*-------------------------------------------------------*
 * Function: SamplesMixNB
 * Purpose: Generate samples from a mixed negative binomial distribution.
 * Input:
 *   probs (NumericVector) - A vector of probabilities for success in each component's negative binomial distribution.
 *   sizes (NumericVector) - A vector of dispersion parameters (sizes) for each component's negative binomial distribution.
 *   weights (NumericVector) - A vector of weights representing the probability of each component being selected.
 *   num_samples (int) - The number of samples to generate.
 *   sample_size (int) - The size of each sample.
 * Output:
 *   List - A list of numeric vectors, each representing a sample drawn from the mixed negative binomial distribution.
 * Description:
 *   This function generates a list of samples where each sample is drawn from a mixed negative binomial distribution.
 *   The mixed distribution is defined by a set of component negative binomial distributions, each with its own probability
 *   and size parameter. The 'weights' vector defines the probability of each component being selected for each trial in a sample.
 * Examples:
 *   // Example usage in R
 *   probs <- c(0.5, 0.3)
 *   sizes <- c(1, 2)
 *   weights <- c(0.6, 0.4)
 *   samples <- SamplesMixNB(probs, sizes, weights, 10, 100)
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List SamplesMixNB(NumericVector probs, NumericVector sizes, NumericVector weights, int num_samples, int sample_size)
{
    int num_components = weights.size();
    List samples(num_samples);
    List labels(num_samples);

    for (int i = 0; i < num_samples; ++i)
    {
        NumericVector sample(sample_size);
        NumericVector label(sample_size);
        for (int j = 0; j < sample_size; ++j)
        {
            int component = Rcpp::sample(num_components, 1, true, weights)[0] - 1; // Select component
            label[j] = component + 1;
            sample[j] = R::rnbinom(sizes[component], probs[component]);
        }
        samples[i] = sample;
        labels[i] = label;
    }

    return List::create(Named("samples") = samples, Named("labels") = labels);
}

/*-------------------------------------------------------*
* Function: SamplesMixZINB
* Purpose: Generate samples from a mixed zero-inflated negative binomial distribution.
* Input:
*   probs (NumericVector) - Vector of probabilities for each component in the mixture.
*   sizes (NumericVector) - Vector of sizes (dispersion parameters) for each component in the mixture.
*   zero (NumericVector) - Vector of zero inflation probabilities for each component.
*   weights (NumericVector) - Vector of weights for each component, determining the mixture proportion.
*   num_samples (int) - Number of samples to be generated.
*   sample_size (int) - Size of each individual sample.
* Output:
*   List - A list where each element is a numeric vector representing a sample from the mixed distribution.
* Description:
*   This function generates samples from a mixed zero-inflated negative binomial distribution.
*   Each component of the mixture has its own size, probability, and zero inflation probability parameters.
*   The mixture of different components is determined by the specified weights.
* Examples:
*   // Example usage in R
*   probs <- c(0.3, 0.6)
*   sizes <- c(1, 2)
*   zero <- c(0.2, 0.4)
*   weights <- c(0.5, 0.5)
*   mixZINB_samples <- SamplesMixZINB(probs, sizes, zero_probs, weights, 10, 100)
*-------------------------------------------------------*/
// [[Rcpp::export]]
List SamplesMixZINB(NumericVector probs, NumericVector sizes, NumericVector zero_probs, NumericVector weights, int num_samples, int sample_size) {
    int num_components = weights.size();
    List samples(num_samples);

    for (int i = 0; i < num_samples; ++i) {
        NumericVector sample(sample_size);
        for (int j = 0; j < sample_size; ++j) {
            int component = Rcpp::sample(num_components, 1, true, weights)[0] - 1;
            sample[j] = (R::runif(0, 1) < zero_probs[component]) ? 0 : R::rnbinom(sizes[component], probs[component]);
        }
        samples[i] = sample;
    }

    return samples;
}