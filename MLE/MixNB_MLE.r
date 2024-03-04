#' Perform inference of mixture of negative binomial models using closed form algorithm
#'
#' This function performs inference for a mixture of negative binomial models using a closed form algorithm.
#' It estimates the parameters of each negative binomial component in the mixture model.
#'
#' @param X A vector of observations.
#' @param nbS An integer specifying the number of states (components) in the mixture model.
#' @param avg An optional vector of initial mean values for each component in the mixture model.
#'            If not provided, k-means algorithm is used to initialize.
#' @param var An optional vector of initial variance values for each component in the mixture model.
#'            If not provided, k-means algorithm is used to initialize.
#' @param iterMax An optional integer specifying the maximum number of iterations for the algorithm.
#'                Defaults to 15000.
#' @param EPSILON An optional value for the convergence threshold used in the stopping criteria.
#'                Defaults to 1e-7.
#' @param verbose Logical; if TRUE, prints the number of iterations and computational time.
#' @return A list containing the estimated parameters: Size, Probability, Mean, Variance, and Weight of each component.
#' @export
#' @examples
#' # Example usage:
#' data(toyexample)
#' resmixNB <- mixNB_MLE(X = toydata, nbS = 3)
#'
mixNB_MLE <- function(X, nbS, avg = NULL, var = NULL, iterMax = 15000, EPSILON = 1e-7, verbose = TRUE) {
    start.time <- Sys.time() # Start timing the computation

    nbT <- length(X) # Number of samples

    # Initialization of parameters
    eWeight <- rep(1 / nbS, nbS) # Initial weight of each component
    esAvg <- numeric(nbS) # Initial mean of each component
    esVar <- numeric(nbS) # Initial variance of each component

    # Use k-means for initial values if not provided
    if (is.null(c(avg, var))) {
        res_kmeans <- kmeans(X, nbS)
        for (i in 1:nbS) {
            esAvg[i] <- res_kmeans$centers[i]
            esVar[i] <- res_kmeans$withinss[i] / (res_kmeans$size[i] - 1)
        }
    } else {
        esAvg <- avg
        esVar <- var
    }

    if (any(esAvg >= esVar)) {
        stop("Invalid data: Variance should be greater than mean for each cluster.")
    }

    # Initial size and probability calculations
    eSize <- esAvg * esAvg / (esVar - esAvg)
    eProb <- esAvg / esVar

    ones_nbT <- matrix(1, nrow = nbT, ncol = 1)
    ones_nbS <- matrix(1, nrow = nbS, ncol = 1)
    ones_nbSS <- matrix(1, nrow = nbS, ncol = nbS)

    # Iteration for parameter estimation
    for (iter in 1:iterMax) {
        prevParam <- cbind(esAvg = esAvg, eSize = eSize, eWeight = eWeight)

        # E-Step: Expectation calculation
        Pmat <- matrix(0, nrow = nbS, ncol = nbT)
        Dmat <- matrix(0, nrow = nbS, ncol = nbT)
        Xmat <- ones_nbS %*% t(X)

        for (i in 1:nbS) {
            for (j in 1:nbT) {
                Dmat[i, j] <- eSize[i] * (digamma(eSize[i] + X[j]) - digamma(eSize[i]))
                Pmat[i, j] <- eWeight[i] * dnbinom(X[j], size = eSize[i], prob = eProb[i])
            }
        }

        Tmat <- Pmat / (ones_nbSS %*% Pmat)
        beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        Bmat <- beta %*% t(ones_nbT)

        # M-Step: Maximization to update parameters
        eWeight <- (Tmat / nbT) %*% ones_nbT # new value of weight

        denominator <- (Tmat * Dmat * Bmat) %*% ones_nbT
        numerator <- (Tmat * (Xmat + Dmat * Bmat - Dmat)) %*% ones_nbT
        eProb <- denominator / numerator # new value of prob

        temp <- ((Tmat * Dmat) %*% ones_nbT) / (Tmat %*% ones_nbT)
        eSize <- -temp / log(eProb) # new value of size

        esAvg <- eSize * (1 - eProb) / eProb # new value of mean
        esVar <- esAvg / eProb # new value of variance

        if (any((is.na(esAvg))) || any((is.na(eSize))) || any((is.na(eWeight)))) {
            warning("mixNB_MLE - NA is appeared!")
            esAvg <- prevParam[, 1]
            eSize <- prevParam[, 2]
            eWeight <- prevParam[, 3]
            break
        }
        
        # reorder estimated parameters
        ordAvg <- order(esAvg)
        esAvg <- esAvg[ordAvg]
        eSize <- eSize[ordAvg]
        eWeight <- eWeight[ordAvg]
        eProb <- eProb[ordAvg]
        esVar <- esVar[ordAvg]
        newParam <- cbind(esAvg = esAvg, eSize = eSize, eWeight = eWeight) # update iterative parameters

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    # Prepare and return results
    parameter <- list(eWeight = as.vector(eWeight), eSize = as.vector(eSize), eProb = as.vector(eProb), esAvg = as.vector(esAvg), esVar = as.vector(esVar))

    res_df <- data.frame(
        Component = seq_along(eWeight),
        Weight = eWeight,
        Size = eSize,
        Probability = eProb,
        Mean = esAvg,
        Variance = esVar
    )

    # Stop timing the computation
    stop.time <- Sys.time()
    minutes <- difftime(stop.time, start.time, units = "min")

    # Verbose output
    if (verbose) {
        cat("Computational iterations:", iter, "\n")
        cat("Computational time:", minutes, "minutes \n")

        # estimated results
        cat("\nEstimate Results:")
        print(knitr::kable(res_df, align = "c", format = "simple"))
    }

    return(parameter)
}
