if (!require("numDeriv")) install.packages("numDeriv", dependencies = TRUE)
library(numDeriv)

# Function to estimate Zero-Inflated Negative Binomial parameters using Newton-Raphson algorithm
zinb_newton <- function(X, starting_values = NULL, iterMax = 500, EPSILON = 1e-5) {
    start_time <- Sys.time()  # Start time
    nbT <- length(X)
    converged <- FALSE
    iterations <- 0 # Sum of inner iterations

    # Initial parameter estimates
    if (is.null(starting_values)) {
        eZero <- sum(X == 0) / length(X)  # Initial estimate for zero inflation
        esAvg <- mean(X)  # Mean of non-zero observations
        esVar <- var(X)  # Variance of non-zero observations
        eSize <- esAvg * esAvg / ((1 - eZero) * (esVar - esAvg) - eZero * esAvg * esAvg)
        eProb <- esAvg / (esVar - eZero / (1 - eZero) * esAvg * esAvg)
    } else {
        eProb <- starting_values[1]
        eSize <- starting_values[2]
        eZero <- starting_values[3]
    }

    params <- c(eProb = eProb, eSize = eSize, eZero = eZero)

    # EM Algorithm
    for (iter in 1:iterMax) {
        # E-Step
        E1Vec <- numeric(nbT)
        E2Vec <- numeric(nbT)
        
        for (i in 1:nbT) {
            E1Vec[i] <- ifelse(X[i] == 0, eZero, 0)
            E2Vec[i] <- (1 - eZero) * dnbinom(X[i], eSize, eProb)
        }
        
        E1Vec <- E1Vec / (E1Vec + E2Vec)
        E2Vec <- E2Vec / (E1Vec + E2Vec)

        # M-Step: Update zero inflation parameter
        eZero <- sum(E1Vec) / sum(E1Vec + E2Vec)
        
        params_nb <- c(eProb, eSize)
        converged_inner <- FALSE
        
        # Inner iteration using Newton-Raphson
        for (iter_inner in 1:iterMax) {
            log_likelihood <- function(params_nb) {
                sum(E2Vec * dnbinom(X, prob = params_nb[1], size = params_nb[2], log = TRUE))
            }
            
            grad <- grad(log_likelihood, params_nb)
            hess <- hessian(log_likelihood, params_nb)
            
            if (any(is.na(hess)) || any(is.infinite(hess))) {
                warning("Hessian matrix contains NA or infinite values. Iteration stopped.")
                break
            }
            
            params_nb_new <- params_nb - solve(hess) %*% grad
            
            if (params_nb_new[1] <= 0 || params_nb_new[1] >= 1 || params_nb_new[2] <= 0) {
                warning("Parameters out of bounds. Iteration stopped.")
                break
            }
            
            if (sum((params_nb_new - params_nb)^2) < EPSILON) {
                converged_inner <- TRUE
                break
            }
            
            params_nb <- params_nb_new
        }
        iterations <- iterations + iter_inner 

        if (converged_inner) {
            eProb <- params_nb[1]
            eSize <- params_nb[2]
        } else {
            warning("M-Step not converged. Iteration stopped.")
            break
        }
        
        params_new <- c(eProb = eProb, eSize = eSize, eZero = eZero)
        
        if (sum((params_new - params)^2) < EPSILON) {
            converged <- TRUE
            break
        }
        
        params <- params_new
    }
    
    end_time <- Sys.time()  # End time
    duration <- difftime(end_time, start_time, units = "secs")  # Compute duration

    return(list(params = params, iterations = iterations + iter , converged = converged, duration = duration))
}

# Function to estimate Zero-Inflated Negative Binomial parameters using an Alternative EM Algorithm
zinb_alem <- function(X, starting_values = NULL, iterMax = 500, EPSILON = 1e-5) {
    start_time <- Sys.time()  # Start time
    nbT <- length(X)
    converged <- FALSE
    
    # Initial parameter estimates
    if (is.null(starting_values)) {
        eZero <- sum(X == 0) / length(X)  # Initial estimate for zero inflation
        esAvg <- mean(X)  # Mean of non-zero observations
        esVar <- var(X)  # Variance of non-zero observations
        eSize <- esAvg * esAvg / ((1 - eZero) * (esVar - esAvg) - eZero * esAvg * esAvg)
        eProb <- esAvg / (esVar - eZero / (1 - eZero) * esAvg * esAvg)
    } else {
        eProb <- starting_values[1]
        eSize <- starting_values[2]
        eZero <- starting_values[3]
    }
    
    params <- c(eProb = eProb, eSize = eSize, eZero = eZero)
    
    # EM Algorithm
    for (iter in 1:iterMax) {
        # E-Step
        beta <- as.vector(1 - 1 / (1 - eProb) - 1 / log(eProb))
        E1Vec <- numeric(nbT)
        E2Vec <- numeric(nbT)
        Dvec <- numeric(nbT)
        
        E1Vec <- ifelse(X == 0, eZero, 0)
        E2Vec <- (1 - eZero) * dnbinom(X, eSize, eProb)
        Dvec <- eSize * (digamma(eSize + X) - digamma(eSize))
        
        E1Vec <- E1Vec / (E1Vec + E2Vec)
        E2Vec <- E2Vec / (E1Vec + E2Vec)
        
        # M-Step
        eZero <- sum(E1Vec) / sum(E1Vec + E2Vec)
        eProb <- sum(E2Vec * (beta * Dvec)) / sum(E2Vec * (X + beta * Dvec - Dvec))
        eSize <- -sum(E2Vec * Dvec) / sum(E2Vec) / log(eProb)
        
        if (eProb >= 1 || eProb <= 0 || eSize <= 0 || eZero >= 1 || eZero <= 0) {
            warning("Parameters out of bounds. Iteration stopped.")
            break
        }
        
        params_new <- c(eProb = eProb, eSize = eSize, eZero = eZero)
        
        if (sum((params_new - params)^2) < EPSILON) {
            converged <- TRUE
            break
        }
        
        params <- params_new
    }
    
    end_time <- Sys.time()  # End time
    duration <- difftime(end_time, start_time, units = "secs")  # Compute duration
    
    return(list(params = params, iterations = iter, converged = converged, duration = duration))
}

compare_zinb <- function(zero_values, prob_values, size_values, starting_values = NULL, num_samples = 200, sample_size = 2000) {
    df <- data.frame()

    for (zero in zero_values) {
        for (prob in prob_values) {
            for (size in size_values) {
                samples <- SamplesZINB(prob, size, zero, num_samples, sample_size)$samples

                new_size <- numeric(num_samples)
                new_prob <- numeric(num_samples)
                new_zero <- numeric(num_samples)
                new_iterations <- numeric(num_samples)
                new_duration <- numeric(num_samples)
                new_converged <- numeric(num_samples)

                alem_size <- numeric(num_samples)
                alem_prob <- numeric(num_samples)
                alem_zero <- numeric(num_samples)
                alem_iterations <- numeric(num_samples)
                alem_duration <- numeric(num_samples)
                alem_converged <- numeric(num_samples)

                for (i in 1:num_samples) {
                    data <- samples[[i]]
                    starting_values <- c(runif(1, 0.01, 0.99), runif(1, 0.1, 50), runif(1, 0.01, 0.99))
                    new_res <- zinb_newton(data, starting_values, iterMax = 500)
                    alem_res <- zinb_alem(data, starting_values, iterMax = 500)

                    new_prob[i] <- new_res$params[1]
                    new_size[i] <- new_res$params[2]
                    new_zero[i] <- new_res$params[3]

                    new_iterations[i] <- new_res$iterations
                    new_duration[i] <- new_res$duration
                    new_converged[i] <- new_res$converged

                    alem_prob[i] <- alem_res$params[1]
                    alem_size[i] <- alem_res$params[2]
                    alem_zero[i] <- alem_res$params[3]

                    alem_iterations[i] <- alem_res$iterations
                    alem_duration[i] <- alem_res$duration
                    alem_converged[i] <- alem_res$converged
                }
                valid_new <- which(new_converged==1)
                valid_alem <- which(alem_converged==1)

                temp_new <- data.frame(
                    "Zero" = zero,
                    "Prob" = prob,
                    "Size" = size,
                    "Method" = "Newton",
                    "Probability" = sum(new_converged)/length(new_converged),
                    "Zero Error mean" = mean(new_zero[valid_new]),
                    "Zero Error dev" = sd(new_zero[valid_new]),
                    "Prob Error mean" = mean(new_prob[valid_new]),
                    "Prob Error dev" = sd(new_prob[valid_new]),
                    "Size Error mean" = mean(new_size[valid_new]),
                    "Size Error mean" = sd(new_size[valid_new]),
                    "CPU time mean" = mean(new_duration[valid_new]),
                    "CPU time dev" = sd(new_duration[valid_new]),
                    "Iterations mean" = mean(new_iterations[valid_new]),
                    "Iterations dev" = sd(new_iterations[valid_new])
                )
                temp_alem <- data.frame(
                    "Zero" = zero,
                    "Prob" = prob,
                    "Size" = size,
                    "Method" = "ALEM",
                    "Probability" = sum(alem_converged)/length(alem_converged),
                    "Zero Error mean" = mean(alem_zero[valid_alem]),
                    "Zero Error dev" = sd(alem_zero[valid_alem]),
                    "Prob Error mean" = mean(alem_prob[valid_alem]),
                    "Prob Error dev" = sd(alem_prob[valid_alem]),
                    "Size Error mean" = mean(alem_size[valid_alem]),
                    "Size Error mean" = sd(alem_size[valid_alem]),
                    "CPU time mean" = mean(alem_duration[valid_alem]),
                    "CPU time dev" = sd(alem_duration[valid_alem]),
                    "Iterations mean" = mean(alem_iterations[valid_alem]),
                    "Iterations dev" = sd(alem_iterations[valid_alem])
                )
                temp <- rbind(temp_new,temp_alem)
                df <- rbind(df, temp)
            }
        }
    }
    return(df)
}

library(extraDistr)
X <- rzinb(n = 2000, size = 10, prob = 0.7, pi = 0.2)
res_new <- zinb_newton(X)
res_alem <- zinb_alem(X)

zero_values <- c(0.2, 0.4, 0.6)
prob_values <- c(0.2, 0.4, 0.6)
size_values <- c(6, 8, 10)

res <- compare_zinb(zero_values, prob_values, size_values)

table <- xtable(res)

set.seed(123)

res1 <- compare_zinb(0.2, 0.2, 2)
res2 <- compare_zinb(0.4, 0.4, 4)
res3 <- compare_zinb(0.6, 0.6, 6)
res <- rbind(res1,res2,res3)
table <- xtable(res)
