ZINB_summary <- function(fit, type = "Wald", verbose) {
    if (type == "Wald") {
        star1 <- as.character(stats::symnum(fit$wald$count$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        star2 <- as.character(stats::symnum(fit$wald$zero$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        Wald_mat_count <- data.frame(cbind(
            fit$coefficients_count,
            fit$wald$count$SE,
            fit$wald$count$pval
        ), star1)
        Wald_mat_zero <- data.frame(cbind(
            fit$coefficients_zero,
            fit$wald$zero$SE,
            fit$wald$zero$pval
        ), star2)
        if (verbose == TRUE) {
            cat("          NB Wald  \n")
            colnames(Wald_mat_count) <- c("Estimation", "SE", "pval", "")
            colnames(Wald_mat_zero) <- c("Estimation", "SE", "pval", "")
            cat("          Count Part  \n")
            print(Wald_mat_count, digits = 3)
            cat("          Zero Part  \n")
            print(Wald_mat_zero, digits = 3)
        }
        return(list(count = Wald_mat_count, zero = Wald_mat_zero))
    } else if (type == "bWald") {
        star1 <- as.character(stats::symnum(fit$bwald$count$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        star2 <- as.character(stats::symnum(fit$bwald$zero$pval,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", " ")
        ))
        Wald_mat_count <- data.frame(cbind(
            fit$coefficients_count,
            fit$bwald$count$SE,
            fit$bwald$count$pval
        ), star1)
        Wald_mat_zero <- data.frame(cbind(
            fit$coefficients_zero,
            fit$bwald$zero$SE,
            fit$bwald$zero$pval
        ), star2)
        if (verbose == TRUE) {
            cat("          NB bWald   \n")
            colnames(Wald_mat_count) <- c("Estimation", "SE", "pval", "")
            colnames(Wald_mat_zero) <- c("Estimation", "SE", "pval", "")
            cat("          Count Part  \n")
            print(Wald_mat_count, digits = 3)
            cat("          Zero Part  \n")
            print(Wald_mat_zero, digits = 3)
        }
        return(list(count = Wald_mat_count, zero = Wald_mat_zero))
    } else {
        stop("Type must be one of 'Wald' or 'bWald'! ")
    }
}

ZINB_Reg_main <- function(Y, X, Z, iterMax = 15000, EPSILON = 1e-5) {
    # Initialization of regression parameters
    esCoef_count <- numeric(ncol(X))
    esCoef_zero <- numeric(ncol(Z))

    # Initialization of zero and probability parameters
    eZero <- sum(Y == 0) / length(Y) # Initial zero inflation coefficient
    nonZeroY <- Y[Y > 0]
    esAvg <- mean(nonZeroY)
    esVar <- var(nonZeroY)
    eProb <- (1 - eZero) * esAvg / esVar

    E_step <- function(Y, X, Z, esCoef_count, esCoef_zero, eProb) {
        # count part
        mu <- exp(X %*% esCoef_count)
        alpha <- -1 / log(eProb)
        gamma <- 1 - 1 / (1 - eProb) - 1 / log(eProb)
        temp <- -eProb * log(eProb) / (1 - eProb)
        lambda <- temp * mu
        eSize <- alpha * lambda
        Dvec <- eSize * (digamma(eSize + Y) - digamma(eSize))

        # zero part
        pi <- 1 / (1 + exp(-Z %*% esCoef_zero))
        E1Vec <- ifelse(Y == 0, pi, 0)
        E2Vec <- (1 - pi) * dnbinom(Y, eSize, eProb)
        xi <- E1Vec / (E1Vec + E2Vec)

        # weight matrix and response vector
        W_count <- diag(as.numeric((1 - xi) * lambda))
        R_count <- log(mu) - 1 + Dvec / lambda
        W_zero <- diag(as.numeric(pi * (1 - pi)))
        R_zero <- log(pi / (1 - pi)) + xi / pi - (1 - xi) / (1 - pi)
        return(list(
            W_count = W_count, R_count = R_count,
            W_zero = W_zero, R_zero = R_zero,
            xi = xi, Dvec = Dvec, gamma = gamma
        ))
    }

    M_step <- function(Y, X, Z, W_count, R_count, W_zero, R_zero, xi, Dvec, gamma) {
        # update prob
        eProb <- sum((1 - xi) * gamma * Dvec) / sum((1 - xi) * (Y + gamma * Dvec - Dvec))
        if (is.na(eProb)) {
            warning("ZINB_EM_IWLS - NA is appeared!")
            eProb <- prevParam["eProb"]
            break
        }
        eProb <- pmin(pmax(eProb, 1e-10), 1 - 1e-10)

        # update regression parameter of count part
        product_count <- t(X) %*% W_count %*% X
        if (det(product_count) <= 0) {
            warning("ZINB_EM_IWLS - Invertible Matrix!")
            break
        }
        esCoef_count <- solve(product_count) %*% t(X) %*% W_count %*% R_count

        # update regression parameter of zero part
        product_zero <- t(Z) %*% W_zero %*% Z
        if (det(product_zero) <= 0) {
            warning("ZINB_EM_IWLS - Invertible Matrix!")
            break
        }
        esCoef_zero <- solve(product_zero) %*% t(Z) %*% W_zero %*% R_zero
        return(list(esCoef_count = esCoef_count, esCoef_zero = esCoef_zero, eProb = eProb))
    }

    # Iteration starts
    for (iter in 1:iterMax) {
        prevParam <- c(esCoef_count, esCoef_zero)

        # E-Step: Expectation calculation
        e_results <- E_step(Y, X, Z, esCoef_count, esCoef_zero, eProb)

        # M-Step: Maximization to update parameters
        m_results <- M_step(Y, X, Z,
            W_count = e_results$W_count, R_count = e_results$R_count,
            W_zero = e_results$W_zero, R_zero = e_results$R_zero,
            xi = e_results$xi, Dvec = e_results$Dvec, gamma = e_results$gamma
        )

        esCoef_count <- m_results$esCoef_count
        esCoef_zero <- m_results$esCoef_zero
        eProb <- m_results$eProb

        newParam <- c(esCoef_count = esCoef_count, esCoef_zero = esCoef_zero)

        # Check for convergence
        D <- sum((newParam - prevParam)^2)
        if (D < EPSILON) {
            break
        }
    }

    wald <- function(Y, X, Z, esCoef_count, esCoef_zero, eProb) {
        e_results <- E_step(Y, X, Z, esCoef_count, esCoef_zero, eProb)
        W_count <- e_results$W_count
        W_zero <- e_results$W_zero
        count_var <- solve(t(X) %*% W_count %*% X)
        zero_var <- solve(t(Z) %*% W_zero %*% Z)
        count_se <- sqrt(diag(count_var))
        zero_se <- sqrt(diag(zero_var))
        count_twald <- esCoef_count^2 / diag(count_var)
        zero_twald <- esCoef_zero^2 / diag(zero_var)
        count_pval <- stats::pchisq(count_twald, df = 1, lower.tail = FALSE)
        zero_pval <- stats::pchisq(zero_twald, df = 1, lower.tail = FALSE)
        return(list(
            count = list(par = esCoef_count, SE = count_se, twald = count_twald, pval = count_pval),
            zero = list(par = esCoef_zero, SE = zero_se, twald = zero_twald, pval = zero_pval)
        ))
    }

    wald <- wald(Y, X, Z, esCoef_count, esCoef_zero, eProb)
    fit <- list()
    fit$coefficients_count <- esCoef_count
    fit$coefficients_zero <- esCoef_zero
    fit$prob <- eProb
    fit$wald <- wald
    fit$convergence <- (iter < iterMax)
    fit$iter <- iter

    return(fit)
}

ZINB_Reg <- function(formula_count, formula_zero, data, iterMax = 15000, EPSILON = 1e-5, bWald_list = list(B = 100), type = "Wald", verbose = TRUE) {
    # Creating model matrix and response vector
    X <- model.matrix(formula_count, data) # Design matrix of count part
    Z <- model.matrix(formula_zero, data) # Design matrix of zero part
    Y <- model.response(model.frame(formula_count, data)) # Response variable
    fit <- ZINB_Reg_main(Y, X, Z)
    esCoef_count <- fit$coefficients_count
    esCoef_zero <- fit$coefficients_zero
    if (type == "bWald") {
        bwald <- function(Y, X, Z, esCoef_count, esCoef_zero, bWald_list = list(B = 100)) {
            B <- bWald_list$B
            n_sample <- length(Y)
            par_b <- sapply(1:B, function(b) {
                resample <- sample(1:n_sample, n_sample, replace = TRUE)
                Y <- Y[resample]
                X <- X[resample, ]
                Z <- Z[resample, ]
                res_b <- ZINB_Reg_main(Y = Y, X = X, Z = Z)
                return(c(
                    as.numeric(res_b$coefficients_count),
                    as.numeric(res_b$coefficients_zero),
                    res_b$wald$count$SE,
                    res_b$wald$zero$SE
                ))
            })
            n_parms_count <- length(esCoef_count)
            n_parms_zero <- length(esCoef_zero)
            n_parms <- n_parms_count + n_parms_zero
            par_SE <- t(par_b)[, -c(1:n_parms)]
            par_SE_count <- par_SE[, c(1:n_parms_count)]
            par_SE_zero <- par_SE[, -c(1:n_parms_count)]
            par_b <- t(par_b)[, c(1:n_parms)]
            par_b_count <- par_b[, c(1:n_parms_count)]
            par_b_zero <- par_b[, -c(1:n_parms_count)]
            varb_count <- stats::var(par_b_count)
            SE_count <- sqrt(diag(varb_count))
            varb_zero <- stats::var(par_b_zero)
            SE_zero <- sqrt(diag(varb_zero))
            twald_count <- esCoef_count^2 / diag(varb_count)
            twald_zero <- esCoef_zero^2 / diag(varb_zero)
            pval_count <- stats::pchisq(twald_count, df = 1, lower.tail = FALSE)
            pval_zero <- stats::pchisq(twald_zero, df = 1, lower.tail = FALSE)
            return(list(
                count = list(
                    B = B,
                    par_b = par_b_count,
                    par_SE = par_SE_count,
                    vcov = varb_count,
                    SE = SE_count,
                    twald = twald_count,
                    pval = pval_count
                ),
                zero = list(
                    B = B,
                    par_b = par_b_zero,
                    par_SE = par_SE_zero,
                    vcov = varb_zero,
                    SE = SE_zero,
                    twald = twald_zero,
                    pval = pval_zero
                )
            ))
        }
        fit$bwald <- bwald(Y, X, Z, esCoef_count, esCoef_zero)
        fit$summary <- ZINB_summary(fit, type = "bWald", verbose)
    } else {
        fit$summary <- ZINB_summary(fit, type = "Wald", verbose)
    }
    fit$call <- match.call()
    class(fit) <- "ZINB_EM_IWLS"
    return(fit)
}

predict <- function(fit, formula_count, formula_zero, data) {
    Z <- model.matrix(formula_zero, data)
    X <- model.matrix(formula_count, data)
    Y <- model.response(model.frame(formula_count, data))

    pi <- 1 / (1 + exp(-Z %*% fit$coefficients_zero))
    mu <- exp(X %*% fit$coefficients_count)
    prob <- fit$prob
    size <- prob / (1 - prob) * mu
    prob_nb <- (1 - pi) * dnbinom(Y, size = size, prob = prob)
    prob_zero <- pi * (Y == 0)
    tau <- prob_zero / (prob_zero + prob_nb)
    pred_source <- ifelse(tau > 0.5, "zero_part", "count_part")
    table <- table(pred_source, data$source)
    accuracy <- sum(pred_source == data$source) / length(pred_source)
    return(list(pred_source = pred_source, table = table, accuracy = accuracy))
}
