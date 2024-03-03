#' Zero-Inflated Negative Binomial Regression Using EM and IWLS
#'
#' This function performs Zero-Inflated Negative Binomial (ZINB) regression using a combination
#' of the Expectation-Maximization (EM) algorithm and Iteratively Reweighted Least Squares (IWLS).
#' It estimates regression coefficients for ZINB distributed data.
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param iterMax Maximum number of iterations for the EM algorithm. Defaults to 15000.
#' @param EPSILON Convergence threshold for the algorithm. If the sum of squared changes 
#'                in the parameters is less than this value, the algorithm stops. Defaults to 1e-7.
#' @return A list containing the estimated regression coefficients, fitted values, residuals, convergence status, number of iterations, and the matched call.
#' @export
#' @examples
#' # Example usage:
#' Y <- rzinb(100, size = 2, prob = 0.5, pi = 0.2)
#' X1 <- rnorm(100)
#' X2 <- rnorm(100)
#' data <- data.frame(Y, X1, X2)
#' result <- ZINB_EM_IWLS(Y ~ X1 + X2, data)
#'
ZINB_EM_IWLS <- function(formula, data, iterMax = 15000, EPSILON = 1e-7) {
    # Creating model matrix and response vector
    X <- model.matrix(formula, data) # Design matrix 
    Y <- model.response(model.frame(formula, data)) # Response variable

    nbT <- length(Y) # Number of samples
    nbC <- ncol(X) # Number of regression coefficients

    # Initialization of parameters
    esCoef <- numeric(nbC) # Initial regression coefficients
    eZero <- sum(Y == 0) / length(Y) # Initial zero inflation coefficient
    nonZeroY <- Y[Y > 0] 
    esAvg <- mean(nonZeroY) 
    esVar <- var(nonZeroY) 

    # Initial estimates for size and probability parameters
    eSize <- esAvg * esAvg / (esVar / (1 - eZero) - esAvg)
    eProb <- (1 - eZero) * esAvg / esVar

    # Iteration starts
    for (iter in 1:iterMax) {
        prevParam <- c(eProb = eProb, eSize = eSize, eZero = eZero, esCoef = esCoef)

        # E-Step: Expectation calculation
        beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        E1Vec <- numeric(nbT)
        E2Vec <- numeric(nbT)
        Dvec <- numeric(nbT)
        # value of mean
        Mvec <- numeric(nbT)
        # elements of Z vector
        Zvec <- numeric(nbT)
        # diagnoal elements of W matrix
        Wmat <- matrix(0, nrow = nbT, ncol = nbT)

        temp <- -eProb * log(eProb) / (1 - eProb)
        for (i in 1:nbT) {
            E1Vec[i] <- ifelse(Y[i] == 0, eZero, 0) # Probability of excess zeros
            E2Vec[i] <- (1 - eZero) * dnbinom(Y[i], eSize, eProb) # Probability from Negative Binomial
            Dvec[i] <- eSize * (digamma(eSize + Y[i]) - digamma(eSize))
        }

        # Normalization of probabilities
        E1Vec <- E1Vec / (E1Vec + E2Vec)
        E2Vec <- E2Vec / (E1Vec + E2Vec)

        for (i in 1:nbT) {
            Mvec[i] <- exp(X[i, ] %*% esCoef)
            Zvec[i] <- log(Mvec[i]) + Dvec[i] / (temp * Mvec[i]) - 1
            Wmat[i, i] <- temp * E2Vec[i] * Mvec[i]
        }

        # M-Step: Maximization to update parameters
        eZero <- sum(E1Vec) / sum(E1Vec + E2Vec) # Update zero inflation coefficient
        eProb <- sum(E2Vec * (beta * Dvec)) / sum(E2Vec * (Y + beta * Dvec - Dvec)) # Update probability parameter
        eSize <- -sum(E2Vec * Dvec) / sum(E2Vec) / log(eProb) # Update size parameter

        if (is.na(eProb) || is.na(eSize) || is.na(eZero)) {
            warning("ZINB_EM_IWLS - NA is appeared!")
            eProb <- prevParam["eProb"]
            eSize <- prevParam["eSize"]
            eZero <- prevParam["eZero"]
            break
        }

        esAvg <- (1 - eZero) * eSize * (1 - eProb) / eProb # Update mean estimate
        esVar <- esAvg / eProb + eZero / (1 - eZero) * esAvg * esAvg # Update variance estimate

        product <- t(X) %*% Wmat %*% X
        if (det(product) <= 0) {
            warning("ZINB_EM_IWLS - Invertible Matrix!")
            break
        }

        esCoef <- solve(product) %*% t(X) %*% Wmat %*% Zvec

        newParam <- c(eProb = eProb, eSize = eSize, eZero = eZero, esCoef = esCoef) 

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    fit <- list()
    fit$coefficients <- esCoef
    fit$fitted.values <- Mvec
    fit$residuals <- esAvg - Mvec
    fit$convergence <- (iter < iterMax)
    fit$iter <- iter
    fit$call <- match.call()

    class(fit) <- "ZINB_EM_IWLS"
    return(fit)
}
