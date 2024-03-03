#' Negative Binomial Regression using EM and IWLS
#'
#' Performs the inference of negative binomial models using an Expectation-Maximization (EM) algorithm combined with Iteratively Reweighted Least Squares (IWLS) for parameter estimation.
#'
#' @param formula An object of class \code{"formula"}: a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param iterMax Maximum number of iterations for the EM algorithm. Defaults to 15000.
#' @param EPSILON The threshold for convergence. If the sum of squared changes in the parameters is less than this value, the algorithm stops. Defaults to 1e-7.
#'
#' @return A list with several components:
#' \itemize{
#'   \item \code{coefficients}: The estimated coefficients of the model.
#'   \item \code{fitted.values}: The fitted mean values for the negative binomial model.
#'   \item \code{residuals}: The residuals of the model.
#'   \item \code{convergence}: Logical value indicating whether the algorithm converged (TRUE) or reached the maximum number of iterations (FALSE).
#'   \item \code{iter}: The number of iterations performed.
#'   \item \code{call}: The matched call.
#' }
#' 
#' @examples
#' # Example usage
#' # n <- 100
#' # Y <- rnbinom(n, size = 2, prob = 0.5)
#' # X1 <- rnorm(n)
#' # X2 <- rnorm(n)
#' # data <- data.frame(Y, X1, X2)
#' # result <- NB_EM_IWLS(Y ~ X1 + X2, data)
#' # summary(result)
#'
#' @export
NB_EM_IWLS <- function(formula, data, iterMax = 15000, EPSILON = 1e-7) {
    # Creating model matrix and response vector
    X <- model.matrix(formula, data) # Design matrix 
    Y <- model.response(model.frame(formula, data)) # Response variable

    nbT <- length(Y) # Number of samples
    nbC <- ncol(X) # Number of regression coefficients

    # Initialization of parameters
    esCoef <- numeric(nbC)
    esAvg <- mean(Y)
    esVar <- var(Y)

    if (esAvg >= esVar) {
        stop("Invalid data: Variance should be greater than mean.")
    }

    # Initial estimates for size and probability parameters
    eSize <- esAvg * esAvg / (esVar - esAvg) 
    eProb <- esAvg / esVar 

    # Iteration starts
    for (iter in 1:iterMax) {
        prevParam <- c(eProb = eProb, eSize = eSize, esCoef = esCoef)

        # E-Step: Expectation calculation
        beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        # value of delta
        Dvec <- numeric(nbT)
        # value of mean
        Mvec <- numeric(nbT)
        # elements of Z vector
        Zvec <- numeric(nbT)
        # diagnoal elements of W matrix
        Wmat <- matrix(0, nrow = nbT, ncol = nbT)

        temp <- -eProb * log(eProb) / (1 - eProb)
        for (i in 1:nbT) {
            Dvec[i] <- eSize * (digamma(eSize + Y[i]) - digamma(eSize))
            Mvec[i] <- exp(X[i, ] %*% esCoef)
            Zvec[i] <- log(Mvec[i]) + Dvec[i] / (temp * Mvec[i]) - 1
            Wmat[i, i] <- temp * Mvec[i]
        }

        # M-Step: Maximization to update parameters
        eProb <- sum(beta * Dvec) / sum(Y + beta * Dvec - Dvec) # new value of prob
        eSize <- -mean(Dvec) / log(eProb) # new value of size

        if (is.na(eProb) || is.na(eSize)) {
            warning("NB_EM_IWLS - NA is appeared!")
            eProb <- prevParam["eProb"]
            eSize <- prevParam["eSize"]
            break
        }

        product <- t(X) %*% Wmat %*% X
        if (det(product) <= 0) {
            warning("NB_EM_IWLS - Invertible Matrix!")
            break
        }

        esCoef <- solve(product) %*% t(X) %*% Wmat %*% Zvec

        # update iterative parameters
        newParam <- c(eProb = eProb, eSize = eSize, esCoef = esCoef) 

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

    class(fit) <- "NB_EM_IWLS"
    return(fit)
}
