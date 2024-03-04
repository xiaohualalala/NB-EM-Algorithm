library(Rcpp)
library(ggplot2)

sourceCpp("Generate_Samples.cpp")
sourceCpp("SimulateNB.cpp")

set.seed(123)


# NB
prob_values <- seq(0.2, 0.8, 0.2)
size_values <- seq(1, 9, 2)
res_NB <- SimulateNBs(prob_values, size_values, 20, 2000)

ggplot(res_NB, aes(x = as.factor(size), y = as.factor(prob))) +
    geom_tile(aes(fill = avg_t)) +
    scale_fill_gradient(low = "#ffeda0", high = "#f03b20") +
    labs(x = "Size", y = "Probability", fill = "Iterations") +
    theme_minimal() +
    ggtitle("Mean of Iterations Across Different Parameters")

## MixNB
source('SimulateMixNB.r')
weights <- c(0.5, 0.5)
means <- c(5, 50)
vars <- c(10, 90)
prob <- means / vars
sizes <- means^2 / (vars - means)

res <- SimulateMixNB(probs, sizes, weights, 1, 2000)






NB <- SamplesNB(0.9, 5, 10, 1000)
ZINB <- SamplesZINB(0.5, 1.0, 0.2, 10, 1000)
probs <- c(0.5, 0.3)
sizes <- c(4, 8)
weights <- c(0.6, 0.4)
zero <- c(0.2, 0.4)
MixNB <- SamplesMixNB(probs, sizes, weights, 1, 2000)
MixZINB <- SamplesMixZINB(probs, sizes, zero, weights, 1, 2000)

# MLE









# ZINB
sourceCpp("/Users/guosa/Desktop/日常/科研/负二项分布/code/cpp_code/Analyse/SimulateZINB.cpp")

# MixNB










MixNB <- SamplesMixNB(probs, sizes, weights, 1, 2000)
observations <- MixNB$samples
true_labels <- MixNB$labels
K <- length(weights)

predicted_labels <- list()
accuracy <- numeric(length(observations))
for (i in 1:length(observations)) {
    observation <- unlist(observations[i])
    true_label <- unlist(true_labels[i])
    predicted_label <- numeric(length(observation))

    res <- mixNB_MLE(observation, K, verbose = FALSE)

    for (j in 1:length(observation)) {
        probability <- sapply(1:K, function(k) {
            dnbinom(observation[j], size = res$eSize[k], prob = res$eProb[k])
        })
        predicted_label[j] <- which.max(probability)
    }

    accuracy[i] <- sum(true_label == predicted_label) / length(true_label)
    predicted_labels[[i]] <- predicted_label

    range <- seq(min(observation), max(observation))

    density <- mixNB_density(range, res$eSize, res$eProb, res$eWeight)

    data <- data.frame(Value = observation)
    density_data <- data.frame(x = range, density = density)

    ggplot(data, aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = diff(range(observation))/40, fill = "skyblue", color = "white") +
    geom_line(data = density_data, aes(x = x, y = density), color = "red", size = 1) +
    theme_minimal() +
    labs(title = "Observations with Estimated MixNB Density", x = "Value", y = "Density")
}

accuracy

# 定义混合负二项分布的密度函数
mixNB_density <- function(x, sizes, probs, weights) {
    K <- length(weights)
    density <- numeric(length(x))

    for (i in 1:length(x)) {
        for (j in 1:K) {
            density[i] <- density[i] + weights[j] * dnbinom(x[i], size = sizes[j], prob = probs[j])
        }
    }
    return(density)
}





library(caret)
library(ggplot2)
confusionMatrix <- confusionMatrix(factor(eslabel), factor(label))







source("/Users/guosa/Desktop/日常/科研/负二项分布/code/NB-EM-ALGORITHM/MLE/MixNB_MLE.r")

# MixZINB
source("/Users/guosa/Desktop/日常/科研/负二项分布/code/cpp_code/MixZINB/MixZINB_MLE.r", chdir = TRUE)
res <- mixZINB_MLE(unlist(MixZINB), 2)

# Reg
# NB
source("/Users/guosa/Desktop/日常/科研/负二项分布/code/cpp_code/NB/NB_Reg.r")
n <- 4000
X1 <- rnorm(n)
X2 <- rnorm(n)

mu <- exp(1 + 0.2 * X1 - 0.4 * X2)
Y <- rnbinom(n, size = 10, mu = mu)
data <- data.frame(Y, X1, X2)
NB_EM_IWLS(Y ~ X1 + X2, data)$coefficients

prob <- exp(1 + 0.2 * X1 - 0.4 * X2) / (1 + exp(1 + 0.2 * X1 - 0.4 * X2))
Y <- rnbinom(n, size = 10, prob = prob)
data <- data.frame(Y, X1, X2)
NB_EM_IWLS(Y ~ X1 + X2, data, type = "probability", link = "logit")$coefficients

# ZINB
source("/Users/guosa/Desktop/日常/科研/负二项分布/code/cpp_code/ZINB/ZINB_Reg.r")
mu <- exp(1 + 0.3 * X1 - 0.4 * X2)
size <- 5
pi <- 0.2
prob <- size * (1 - pi) / (mu + size * (1 - pi))
Y <- rzinb(n, size = size, prob = prob, pi = pi)
data <- data.frame(Y, X1, X2)
ZINB_EM_IWLS(Y ~ X1 + X2, data)$coefficients

prob <- exp(1 + 0.2 * X1 - 0.4 * X2) / (1 + exp(1 + 0.2 * X1 - 0.4 * X2))
Y <- rzinb(n, size = 10, prob = prob, pi = 0.4)
data <- data.frame(Y, X1, X2)
ZINB_EM_IWLS(Y ~ X1 + X2, data, type = "probability", link = "logit")$coefficients
