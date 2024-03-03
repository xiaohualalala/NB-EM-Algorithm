#' Perform Inference for Mixture of Negative Binomial Models Using EM and IWLS Algorithms
#'
#' This function performs inference for a mixture of negative binomial models using a combination
#' of the Expectation-Maximization (EM) algorithm and Iteratively Reweighted Least Squares (IWLS).
#' It estimates the regression coefficients and other parameters for each component in the mixture model.
#'
#' @param formula An object of class \code{"formula"}: a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param nbS An integer specifying the number of states (components) in the mixture model.
#' @param avg Optional vector of initial mean values for each component in the mixture model.
#'            If not provided, k-means algorithm is used to initialize.
#' @param var Optional vector of initial variance values for each component in the mixture model.
#'            If not provided, k-means algorithm is used to initialize.
#' @param iterMax An integer specifying the maximum number of iterations for the EM algorithm.
#'                Defaults to 15000.
#' @param EPSILON A convergence threshold for the algorithm.
#'                If the sum of squared changes in the parameters is less than this value, the algorithm stops. Defaults to 1e-7.
#' @return A list containing the estimated regression coefficients, fitted values, residuals,
#'         convergence status, number of iterations, and the matched call.
#' @export
#' @examples
#' data(toyexample)
#' resMixNB <- MixNB_EM_IWLS(Y ~ X1 + X2, data = toydata, nbS = 3)
#'
MixNB_EM_IWLS <- function(formula, data, nbS, avg = NULL, var = NULL, iterMax = 15000, EPSILON = 1e-7) {
    # Creating model matrix and response vector
    X <- model.matrix(formula, data) # Design matrix 
    Y <- model.response(model.frame(formula, data)) # Response variable

    nbT <- length(Y) # Number of samples
    nbC <- ncol(X) # Number of regression coefficients

    # Initialization of parameters
    esCoef <- matrix(0, nrow = nbS, ncol = nbC) # Initial regression coefficients
    eWeight <- rep(1 / nbS, nbS) # Initial weights for each component
    esAvg <- numeric(nbS) # Initial means for each component
    esVar <- numeric(nbS) # Initial variances for each component

    # Use k-means for initial values if not provided
    if (is.null(c(avg, var))) {
        res_kmeans <- kmeans(X, nbS)
        for (i in 1:nbS) {
            esAvg[i] <- res_kmeans$centers[i]
            esVar[i] <- as.matrix(res_kmeans$withinss[i] / (res_kmeans$size[i] - 1))
        }
    } else {
        esAvg <- avg
        esVar <- var
    }

    if (any(esAvg >= esVar)) {
        stop("Invalid data: Variance should be greater than mean for each cluster.")
    }

    # Initial estimates for size and probability parameters
    eSize <- esAvg * esAvg / (esVar - esAvg)
    eProb <- esAvg / esVar

    ones_nbT <- matrix(1, nrow = nbT, ncol = 1)
    ones_nbS <- matrix(1, nrow = nbS, ncol = 1)
    ones_nbSS <- matrix(1, nrow = nbS, ncol = nbS)
    ones_nbST <- matrix(1, nrow = nbS, ncol = nbT)

    # Iteration starts
    for (iter in 1:iterMax) {
        prevParam <- cbind(eProb = eProb, eSize = eSize, eWeight = eWeight, esCoef = esCoef)

        # E-Step: Expectation calculation
        Pmat <- matrix(0, nrow = nbS, ncol = nbT)
        Dmat <- matrix(0, nrow = nbS, ncol = nbT)
        Mmat <- matrix(0, nrow = nbS, ncol = nbT)
        Wmat <- matrix(0, nrow = nbS, ncol = nbT)
        Zmat <- matrix(0, nrow = nbS, ncol = nbT)

        for (i in 1:nbS) {
            for (j in 1:nbT) {
                Pmat[i, j] <- eWeight[i] * dnbinom(Y[j], size = eSize[i], prob = eProb[i])
                Dmat[i, j] <- eSize[i] * (digamma(eSize[i] + Y[j]) - digamma(eSize[i]))
            }
        }

        Ymat <- ones_nbS %*% t(Y)
        Tmat <- Pmat / (ones_nbSS %*% Pmat)
        A <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        Bmat <- A %*% t(ones_nbT)

        Mmat <- exp(esCoef %*% t(X))
        temp1 <- (-eProb * log(eProb) / (1 - eProb)) %*% t(ones_nbT)
        Wmat <- temp1 * Tmat * Mmat
        Zmat <- log(Mmat) + Dmat / (temp1 * Mmat) - ones_nbST

        # M-Step: Maximization to update parameters
        eWeight <- (Tmat / nbT) %*% ones_nbT # new value of weight

        denominator <- (Bmat * Tmat * Dmat) %*% ones_nbT
        numerator <- (Tmat * (Ymat + Dmat * Bmat - Dmat)) %*% ones_nbT
        eProb <- denominator / numerator # new value of prob

        temp2 <- ((Tmat * Dmat) %*% ones_nbT) / (Tmat %*% ones_nbT)
        eSize <- -temp2 / log(eProb) # new value of size

        esAvg <- eSize * (1 - eProb) / eProb # new value of mean
        esVar <- esAvg / eProb # new value of variance

        # new value of regression coefficient
        for (i in 1:nbS) {
            product <- t(X) %*% diag(Wmat[i, ]) %*% X
            if (det(product) <= 0) {
                warning("NB_EM_IWLS - Invertible Matrix!")
                break
            }
            esCoef[i, ] <- t(solve(product) %*% t(X) %*% diag(Wmat[i, ]) %*% Zmat[i, ])
        }

        # update iterative parameters
        newParam <- c(eProb = eProb, eSize = eSize, eWeight = eWeight, esCoef = esCoef)

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    fit <- list()
    fit$coefficients <- data.frame(
        Component = 1:nbS,
        Intercept = esCoef[, 1],
        X1 = esCoef[, 2],
        X2 = esCoef[, 3]
    )
    fit$fitted.values <- Mmat
    fit$residuals <- esAvg %*% t(ones_nbT) - fit$fitted.values
    fit$convergence <- (iter < iterMax)
    fit$iter <- iter
    fit$call <- match.call()

    class(fit) <- "MixNB_EM_IWLS"
    return(fit)
}
