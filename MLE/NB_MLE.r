#' Maximum Likelihood Estimation for Negative Binomial Distribution
#'
#' Performs maximum likelihood estimation (MLE) for the parameters of a negative binomial distribution via closed form algorithm.
#'
#' @param X A numeric vector of observations.
#' @param iterMax Maximum number of iterations for the MLE algorithm. Defaults to 15000.
#' @param EPSILON Convergence threshold for the algorithm. If the sum of squared changes in the parameters is less than this value, the algorithm stops. Defaults to 1e-7.
#' @param verbose Logical; if TRUE, prints the number of iterations and computational time.
#' @param se Logical; if TRUE, calculates and returns the standard errors of the parameters.
#' @param rate Logical; if TRUE, calculates and returns the convergence rates of the parameters.
#'
#' @return A list with several components:
#' \itemize{
#'   \item \code{Parameter}: A list containing the estimated parameters (probability and size).
#'   \item \code{SE}: (Optional) A list containing the standard errors of the estimated parameters.
#'   \item \code{Rate}: (Optional) A list containing the convergence rates of the estimated parameters.
#' }
#'
#' @examples
#' # Example usage
#' # set.seed(123)
#' # X <- rnbinom(100, size = 2, prob = 0.5)
#' # result <- NB_MLE(X)
#' # print(result)
#'
#' @export
NB_MLE <- function(X, iterMax = 15000, EPSILON = 1e-5, verbose = TRUE, se = FALSE, rate = FALSE) {
    start.time <- Sys.time() # Start timing the computation

    nbT <- length(X) # Number of samples

    # Initialization of parameters
    esAvg <- mean(X) # Mean of the observations
    esVar <- var(X) # Variance of the observations

    # Check for valid input data
    if (esAvg >= esVar) {
        stop("Invalid data: Variance should be greater than mean.")
    }

    # Initial estimates for the size and probability parameters
    eSize <- esAvg * esAvg / (esVar - esAvg)
    eProb <- esAvg / esVar

    # Iteration starts
    for (iter in 1:iterMax) {
        prevParam <- c(esAvg = esAvg, eSize = eSize) # Store previous parameter values

        ## E-Step: Expectation calculation
        beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        Dvec <- numeric(nbT) # Initialize vector for storing intermediate values
        for (i in 1:nbT) {
            Dvec[i] <- eSize * (digamma(eSize + X[i]) - digamma(eSize))
        }

        ## M-Step: Maximization step to update the parameters
        eProb <- sum(beta * Dvec) / sum(X + beta * Dvec - Dvec)
        eSize <- -mean(Dvec) / log(eProb)

        # Update mean and variance estimates
        esAvg <- eSize * (1 - eProb) / eProb
        esVar <- esAvg / eProb

        # Store new parameters
        newParam <- c(esAvg = esAvg, eSize = eSize)

        # Check for NaN values
        if (is.na(esAvg) || is.na(eSize)) {
            warning("NB_MLE - NA is appeared!")
            esAvg <- prevParam["esAvg"]
            eSize <- prevParam["eSize"]
            break
        }

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }
    
    # Store the final parameters
    parameter <- list(eProb = eProb, eSize = eSize, esAvg = esAvg, esVar = esVar)

    # Prepare result list
    res <- list()
    res$Parameter <- parameter

    # Create a data frame for the estimates
    res_df <- data.frame(
        Parameter = c("Prob", "Size"),
        Estimate = c(eProb, eSize)
    )

    # Stop timing the computation
    stop.time <- Sys.time()
    minutes <- difftime(stop.time, start.time, units = "min")
    
    # Verbose output
    if (verbose) {
        cat("Computational iterations:", iter, "\n")
        cat("Computational time:", minutes, "minutes \n")
    }

    # Standard Error calculation
    if (se) {
        I <- Infor_Matrix_O(X, Param = parameter)
        if (det(I) <= 0) {
            warning("NB_MLE - Fisher Information Matrix is invertible!")
        } else {
            inverse <- solve(I)
            SE_eProb <- sqrt(inverse[1, 1])
            SE_eSize <- sqrt(inverse[2, 2])
            SE <- list(SE_eProb = SE_eProb, SE_eSize = SE_eSize)
            res$SE <- SE
            res_df$Standard_Error <- c(SE_eProb, SE_eSize)
        }
    }

    # Convergence Rate calculation
    if (rate) {
        J <- Conver_Matrix(X, Param = parameter)
        Rate_eProb <- J[1, 1]
        Rate_eSize <- J[2, 2]
        Rate_global <- max(eigen(J)$values)
        Rate <- list(Rate_eProb = Rate_eProb, Rate_eSize = Rate_eSize, Rate_global = Rate_global)
        res$Rate <- Rate
        res_df$Convergence_Rate <- c(Rate_eProb, Rate_eSize)
    }

    # Output estimated results
    cat("\nEstimate Results:")
    print(knitr::kable(res_df, align = "c", format = "simple"))

    return(res)
}