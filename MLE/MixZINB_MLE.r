#' Perform inference of mixture of zero inflated negative binomial models using closed form algorithm
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @param zero a vector of initial value of zero inflation coefficient (optional).
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm (optional).
#' @param EPSILON a value for the threshold used for the stopping criteria for the mixNB_MLE algorithm (optional).
#' @return A list with the inferred parameters: Size, Probability, Mean, Variance, Weight, Zero Inflation Coefficient.
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of mixture of negative binomial models using closed form algorithm
#' resmixZINB <- mixZINB_MLE(X = toydata, nbS = 3)
#'
mixZINB_MLE <- function(X, nbS, avg = NULL, var = NULL, zero = NULL, iterMax = 15000, EPSILON = 1e-5, verbose = TRUE) {
    start.time <- Sys.time() # Start timing the computation

    nbT <- length(X) # Number of samples

    # Initialization of parameters
    eWeight <- rep(1 / nbS, nbS) # Initial weight of each component
    eZero <- numeric(nbS) 
    esAvg <- numeric(nbS)
    esVar <- numeric(nbS)

    if (is.null(zero)) {
        eZero <- rep(sum(X == 0) / nbT, nbS) # Intial value of zero inflation coefficient 
    } else {
        eZero <- zero
    }

    # Use k-means for initial values if not provided
    if (is.null(c(avg, var))) {
        nonZeroX <- X[X > 0]
        res_kmeans <- kmeans(nonZeroX, nbS)
        for (i in 1:nbS) {
            esAvg[i] <- res_kmeans$centers[i]
            esVar[i] <- res_kmeans$withinss[i] / (res_kmeans$size[i] - 1)
        }
    } else {
        esAvg <- avg
        esVar <- var
    }

    if (any(esAvg * (1 - eZero) >= esVar)) {
        stop("Invalid data: Variance should be greater than mean for each cluster.")
    }

    # Initial size and probability calculations
    eSize <- esAvg * esAvg / (esVar / (1 - eZero) - esAvg)
    eProb <- (1 - eZero) * esAvg / esVar

    ones_nbT <- matrix(1, nrow = nbT, ncol = 1)
    ones_nbS <- matrix(1, nrow = nbS, ncol = 1)
    ones_nbSS <- matrix(1, nrow = nbS, ncol = nbS)

    # Iteration for parameter estimation
    for (iter in 1:iterMax) {
        prevParam <- cbind(esAvg = esAvg, eSize = eSize, eWeight = eWeight, eZero = eZero)

        # E-Step: Expectation calculation
        Pmat <- matrix(0, nrow = nbS, ncol = nbT)
        Dmat <- matrix(0, nrow = nbS, ncol = nbT)
        E1mat <- matrix(0, nrow = nbS, ncol = nbT)
        E2mat <- matrix(0, nrow = nbS, ncol = nbT)
        Xmat <- ones_nbS %*% t(X)

        for (i in 1:nbS) {
            for (j in 1:nbT) {
                Dmat[i, j] <- eSize[i] * (digamma(eSize[i] + X[j]) - digamma(eSize[i]))
                E1mat[i, j] <- ifelse(X[j] == 0, eZero[i], 0)
                E2mat[i, j] <- (1 - eZero[i]) * dnbinom(X[j], size = eSize[i], prob = eProb[i])
                Pmat[i, j] <- eWeight[i] * (E1mat[i, j] + E2mat[i, j])
            }
        }

        Tmat <- Pmat / (ones_nbSS %*% Pmat)
        beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        Bmat <- beta %*% t(ones_nbT)
        E1mat <- E1mat / (E1mat + E2mat) 
        E2mat <- E2mat / (E1mat + E2mat)

        # M-Step: Maximization to update parameters
        eWeight <- (Tmat / nbT) %*% ones_nbT # new value of weight

        eZero <- ((Tmat * E1mat) %*% ones_nbT) / ((Tmat * (E1mat + E2mat)) %*% ones_nbT) # new value of zero inflation coefficients

        denominator <- (Tmat * E2mat * Dmat * Bmat) %*% ones_nbT
        numerator <- (Tmat * E2mat * (Xmat + Dmat * Bmat - Dmat)) %*% ones_nbT
        eProb <- denominator / numerator # new value of prob

        temp <- ((Tmat * E2mat * Dmat) %*% ones_nbT) / ((Tmat * E2mat) %*% ones_nbT)
        eSize <- -temp / log(eProb) # new value of size

        esAvg <- (1 - eZero) * eSize * (1 - eProb) / eProb # new value of mean
        esVar <- esAvg / eProb + eZero / (1 - eZero) * esAvg * esAvg # new value of variance

        if (any((is.na(esAvg))) || any((is.na(eSize))) || any((is.na(eWeight))) || any((is.na(eZero)))) {
            warning("mixZINB_MLE - NA is appeared!")
            esAvg <- prevParam[, 1]
            eSize <- prevParam[, 2]
            eWeight <- prevParam[, 3]
            eZero <- prevParam[, 4]
            break
        }

        newParam <- cbind(esAvg = esAvg, eSize = eSize, eWeight = eWeight, eZero = eZero) # update iterative parameters

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    # Prepare and return results
    parameter <- list(eWeight = as.vector(eWeight), eZero = as.vector(eZero), eSize = as.vector(eSize), eProb = as.vector(eProb), esAvg = as.vector(esAvg), esVar = as.vector(esVar))

    res_df <- data.frame(
        Component = seq_along(eWeight),
        Weight = eWeight,
        Size = eSize,
        Probability = eProb,
        Zero = eZero,
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
    }

    # Output estimated results
    cat("\nEstimate Results:")
    print(knitr::kable(res_df, align = "c"))

    return(parameter)
}