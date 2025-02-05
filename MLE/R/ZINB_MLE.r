#' Maximum Likelihood Estimation for Zero-Inflated Negative Binomial Distribution
#'
#' This function performs maximum likelihood estimation for a Zero-Inflated Negative Binomial (ZINB) distribution.
#' It iteratively estimates the zero inflation coefficient, size parameter, and probability parameter of the distribution.
#'
#' @param X A numeric vector of observations.
#' @param iterMax An integer specifying the maximum number of iterations for the algorithm. Defaults to 15000.
#' @param EPSILON A convergence threshold for the algorithm. If the sum of squared changes in the parameters is less than this value, the algorithm stops. Defaults to 1e-7.
#' @param verbose Logical; if TRUE, prints the number of iterations and computational time.
#' @return A list containing the estimated parameters: Zero inflation coefficient, size, probability, mean, and variance.
#' @export
#' @examples
#' # Example usage:
#' library(extraDistr)
#' data <- rzinb(n = 1000, size = 10, prob = 0.5, pi = 0.2)
#' resZINB <- ZINB_MLE(X = data)
#'
ZINB_MLE <- function(X, intial_values = NULL, iterMax = 15000, EPSILON = 1e-7, verbose = TRUE) {
    start.time <- Sys.time() # Start timing the computation

    nbT <- length(X) # Number of observations
    if(is.null(intial_values)){
       # Initialization of parameters
       eZero <- sum(X == 0) / length(X) # Initial estimate for zero inflation coefficient
       nonZero <- X[X > 0]
       esAvg <- mean(nonZero) # Mean of observations
       esVar <- var(nonZero) # Variance of observations
       # Initial estimates for size and probability parameters
       eSize <- esAvg * esAvg / (esVar - esAvg)
       eProb <- esAvg / esVar
    } else {
       eZero <- intial_values[1]
       eSize <- intial_values[2]
       eProb <- intial_values[3]
       esAvg <- eSize * (1 - eProb) / eProb
       esVar <- esAvg / eProb
    }


    # Iteration for parameter estimation
    for (iter in 1:iterMax) {
        prevParam <- c(eProb = eProb, eSize = eSize, eZero = eZero)

        # E-Step: Expectation calculation
        beta <- as.vector(1 - 1 / (1 - eProb) - 1 / log(eProb))
        E1Vec <- numeric(nbT)
        E2Vec <- numeric(nbT)
        Dvec <- numeric(nbT)

        for (i in 1:nbT) {
            E1Vec[i] <- ifelse(X[i] == 0, eZero, 0) # Probability of excess zeros
            E2Vec[i] <- (1 - eZero) * dnbinom(X[i], eSize, eProb) # Probability from Negative Binomial
            Dvec[i] <- eSize * (digamma(eSize + X[i]) - digamma(eSize))
        }

        # Normalization of probabilities
        E1Vec <- E1Vec / (E1Vec + E2Vec)
        E2Vec <- E2Vec / (E1Vec + E2Vec)

        # M-Step: Maximization to update parameters
        eZero <- sum(E1Vec) / sum(E1Vec + E2Vec) # Update zero inflation coefficient
        eProb <- sum(E2Vec * (beta * Dvec)) / sum(E2Vec * (X + beta * Dvec - Dvec)) # Update probability parameter
        eSize <- -sum(E2Vec * Dvec) / sum(E2Vec) / log(eProb) # Update size parameter
        esAvg <- (1 - eZero) * eSize * (1 - eProb) / eProb # Update mean estimate
        esVar <- esAvg / eProb + eZero / (1 - eZero) * esAvg * esAvg # Update variance estimate

        newParam <- c(eProb = eProb, eSize = eSize, eZero = eZero) # Combine new parameters

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    # Prepare and return results
    parameter <- list(eZero = eZero, eSize = eSize, eProb = eProb, esAvg = esAvg, esVar = esVar)
    res <- list()
    res$Parameter <- parameter

    # Create a data frame for the estimates
    res_df <- data.frame(
        Size = eSize,
        Probability = eProb,
        Zero = eZero,
        Mean = esAvg,
        Variance = esVar
    )

    # Verbose output
    if (verbose) {
        stop.time <- Sys.time()
        minutes <- difftime(stop.time, start.time, units = "min")
        cat("Computational iterations:", iter, "\n")
        cat("Computational time:", minutes, "minutes \n")
    }

    # Output estimated results
    cat("\nEstimate Results:")
    print(knitr::kable(res_df, align = "c"))

    return(res)
}
