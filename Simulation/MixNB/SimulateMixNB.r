source("../../MLE/R/MixNB_MLE.r")

# Define the density function for a mixed Negative Binomial distribution
mixNB_density <- function(x, sizes, probs, weights) {
    # Initialize the density vector
    density <- numeric(length(x))
    # Compute the mixed density
    for (i in 1:length(x)) {
        for (j in 1:length(weights)) {
            # Add the weighted density of each Negative Binomial component
            density[i] <- density[i] + weights[j] * dnbinom(x[i], size = sizes[j], prob = probs[j])
        }
    }
    return(density)
}

# Function to analyze mixed Negative Binomial distribution data
SimulateMixNB_r <- function(probs, sizes, weights, sample_size) {
    # Generate samples from the mixed Negative Binomial distribution
    MixNB <- SamplesMixNB(probs, sizes, weights, 1, sample_size)
    observations <- unlist(MixNB$samples)
    true_labels <- unlist(MixNB$labels)
    K <- length(weights) # Number of components in the mixture
    N <- length(observations)
    predicted_labels <- numeric(N)

    # Estimate parameters using Maximum Likelihood Estimation (MLE)
    estimation <- mixNB_MLE(observations, K)

    # Predict labels
    for (i in 1:N) {
        # Calculate the probability of each observation belonging to each component
        probability <- sapply(1:K, function(k) {
            dnbinom(observations[i], size = estimation$eSize[k], prob = estimation$eProb[k])
        })
        # Assign label based on highest probability
        predicted_labels[i] <- which.max(probability)
    }

    # Calculate accuracy
    accuracy <- sum(true_labels == predicted_labels) / N

    return(list(
        observations = observations,
        estimation = estimation,
        true_labels = true_labels,
        predicted_labels = predicted_labels,
        accuracy = accuracy
    ))
}