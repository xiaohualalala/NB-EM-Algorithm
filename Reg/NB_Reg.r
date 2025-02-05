#' Negative Binomial Regression using EM and IWLS
#'
#' Performs inference for negative binomial models using an Expectation-Maximization (EM) algorithm
#' combined with Iteratively Reweighted Least Squares (IWLS) for parameter estimation. The algorithm
#' iteratively estimates the regression coefficients.
#'
#' @param formula A formula object describing the model to be fitted. The formula should include the
#' response variable and explanatory variables.
#' @param data A data frame containing the variables used in the model. This data frame must contain
#' all variables specified in the formula.
#' @param iterMax Maximum number of iterations for the EM algorithm. Defaults to 15000. The algorithm will stop
#' if the maximum number of iterations is reached.
#' @param EPSILON Convergence threshold. If the sum of squared changes in the parameters between iterations
#' is less than this value, the algorithm stops. Defaults to 1e-7.
#'
#' @return A list containing several components:
#' \itemize{
#'   \item \code{coefficients}: The estimated coefficients of the regression model.
#'   \item \code{r_squared}: The R-squared value for the model, representing the proportion of variance explained by the model.
#'   \item \code{convergence}: A logical value indicating whether the algorithm converged (TRUE) or reached the maximum number of iterations (FALSE).
#'   \item \code{iter}: The number of iterations performed by the EM algorithm.
#'   \item \code{call}: The matched call to the function, showing the formula and data used.
#' }
#'
#' @details
#' The algorithm first initializes the parameters using basic moments (mean and variance) of the response variable.
#' Then, it iterates through the Expectation (E) step, where the necessary quantities are computed based on current
#' parameter estimates, followed by the Maximization (M) step, where the parameters are updated using Iteratively Reweighted
#' Least Squares (IWLS). The process continues until either convergence is reached (i.e., the parameter change is smaller than \code{EPSILON})
#' or the maximum number of iterations (\code{iterMax}) is reached.
#'
#' @export

# Calculate p-values
NB_summary <- function(fit, type = "Wald", verbose) {
    if (type == "Wald") {
        star1 <- as.character(stats::symnum(fit$wald$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        Wald_mat <- data.frame(cbind(
            fit$coefficients,
            fit$wald$SE,
            fit$wald$pval
        ), star1)
        if (verbose == TRUE) {
            cat("          NB Wald  \n")
            colnames(Wald_mat) <- c("Estimation", "SE", "pval", "")
            print(Wald_mat, digits = 3)
        }
        return(Wald_mat)
    } else if (type == "bWald") {
        star1 <- as.character(stats::symnum(fit$bwald$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        Wald_mat <- data.frame(cbind(
            fit$coefficients,
            fit$bwald$SE,
            fit$bwald$pval
        ), star1)
        if (verbose == TRUE) {
            cat("          NB bWald   \n")
            colnames(Wald_mat) <- c("Estimation", "SE", "pval", "")
            print(Wald_mat, digits = 3)
        }
        return(Wald_mat)
    } else {
        stop("Type must be one of 'Wald' or 'bWald'! ")
    }
}
# Compute R-squared
r_squared <- function(Y, X, esCoef) {
    mu <- exp(X %*% esCoef) # Compute predicted values
    ss_total <- sum((Y - mean(Y))^2) # Total sum of squares
    ss_residual <- sum((Y - mu)^2) # Residual sum of squares
    r_squared <- 1 - (ss_residual / ss_total) # Compute R-squared
    return(r_squared)
}

# Negative Binomial Regression function
NB_Reg_main <- function(Y, X, iterMax = 15000, EPSILON = 1e-7) {
    # Initialize regression parameters
    esCoef <- numeric(ncol(X))

    # Initial probability parameter estimates
    esAvg <- mean(Y)
    esVar <- var(Y)
    if (esAvg >= esVar) {
        stop("Invalid data: Variance should be greater than mean.")
    }
    eProb <- esAvg / esVar
    # E-step
    E_step <- function(Y, X, esCoef, eProb) {
        # E-step: Expectation calculation
        mu <- exp(X %*% esCoef) # Compute the expected value

        alpha <- -1 / log(eProb)
        gamma <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        temp <- -eProb * log(eProb) / (1 - eProb)
        lambda <- temp * mu
        eSize <- alpha * lambda
        Dvec <- eSize * (digamma(eSize + Y) - digamma(eSize))

        # Weight matrix and response vector
        Wmat <- diag(as.numeric(lambda)) # Weight matrix
        Rvec <- log(mu) - 1 + Dvec / lambda # Response vector

        return(list(Wmat = Wmat, Rvec = Rvec, Dvec = Dvec, gamma = gamma))
    }

    # M-step
    M_step <- function(Y, X, Wmat, Rvec, Dvec, gamma) {
        # M-step: Maximization to update parameters
        eProb <- sum(gamma * Dvec) / sum(Y + gamma * Dvec - Dvec) # Compute eProb
        if (is.na(eProb)) {
            warning("NB_EM_IWLS - NA appeared!")
            break
        }
        eProb <- pmin(pmax(eProb, 1e-10), 1 - 1e-10) # Ensure eProb is within valid bounds

        product <- t(X) %*% Wmat %*% X # Compute matrix product
        if (det(product) <= 0) {
            warning("NB_EM_IWLS - Matrix is not invertible!")
            break
        }
        esCoef <- solve(product) %*% t(X) %*% Wmat %*% Rvec # Solve for regression coefficients
        return(list(esCoef = esCoef, eProb = eProb)) # Return updated coefficients and eProb
    }
    
    # Start iterations
    for (iter in 1:iterMax) {
        prevParam <- c(esCoef = esCoef) # Save previous parameters

        # E-step: Expectation calculation
        e_results <- E_step(Y, X, esCoef, eProb)

        # M-step: Maximization to update parameters
        m_results <- M_step(Y, X,
            Wmat = e_results$Wmat,
            Rvec = e_results$Rvec,
            Dvec = e_results$Dvec,
            gamma = e_results$gamma
        )

        esCoef <- m_results$esCoef # Update regression coefficients
        eProb <- m_results$eProb # Update eProb

        # Update iterative parameters
        newParam <- c(esCoef = esCoef)

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break # Exit if converged
        }
    }
    wald <- function(X, esCoef, eProb) {
        mu <- exp(X %*% esCoef)
        Wmat <- diag(as.numeric(-eProb * log(eProb) / (1 - eProb) * mu))
        coef_var <- solve(t(X) %*% Wmat %*% X) # Compute the variance of the coefficients
        coef_se <- sqrt(diag(coef_var)) # Standard errors
        zval <- esCoef / coef_se
        twald <- esCoef^2 / diag(coef_var)
        pval <- stats::pchisq(twald, df = 1, lower.tail = FALSE)
        return(list(
            par = esCoef,
            SE = coef_se,
            zval = zval,
            twald = twald,
            pval = pval
        ))
    }
    wald <- wald(X, esCoef, eProb)
    # Create result list
    fit <- list()
    fit$coefficients <- esCoef
    fit$prob <- eProb
    fit$wald <- wald
    fit$convergence <- (iter < iterMax)
    fit$iter <- iter

    return(fit) # Return the fitted model
}

NB_Reg <- function(formula, data, iterMax = 15000, EPSILON = 1e-7, bWald_list = list(B = 100), type = "Wald", verbose = TRUE) {
    X <- model.matrix(formula, data) # Design matrix
    Y <- model.response(model.frame(formula, data)) # Response variable

    fit <- NB_Reg_main(Y = Y, X = X)
    esCoef <- fit$coefficients
    r_squared <- r_squared(Y, X, esCoef)
    fit$r_squared <- r_squared
    bwald <- function(Y, X, esCoef, bWald_list = list(B = 100)) {
        B <- bWald_list$B
        n_sample <- length(Y)
        par_b <- sapply(1:B, function(b) {
            resample <- sample(1:n_sample, n_sample, replace = TRUE)
            Y <- Y[resample]
            X <- X[resample, ]
            res_b <- NB_Reg_main(Y = Y, X = X)
            return(c(as.numeric(res_b$coefficients), res_b$wald$SE))
        })
        n_parms <- length(esCoef)
        par_SE <- t(par_b)[, -c(1:n_parms)]
        par_b <- t(par_b)[, 1:n_parms]
        varb <- stats::var(par_b)
        SE <- sqrt(diag(varb))
        twald <- esCoef^2 / diag(varb)
        pval <- stats::pchisq(twald, df = 1, lower.tail = FALSE)
        return(list(
            B = B,
            par_b = par_b,
            par_SE = par_SE,
            vcov = varb,
            SE = SE,
            twald = twald,
            pval = pval
        ))
    }
    if (type == "bWald") {
        fit$bwald <- bwald(Y, X, esCoef)
        fit$summary <- NB_summary(fit, type = "bWald", verbose)
    } else {
        fit$summary <- NB_summary(fit, type = "Wald", verbose)
    }
    fit$call <- match.call()
    class(fit) <- "NB_EM_IWLS" # Set class for the result
    return(fit) # Return the fitted model
}
