# NB Reg
setwd("/Users/guosa/Desktop/日常/科研/负二项分布/code/NB-EM-Algorithm/Simulation")
source('../Reg/NB_Reg.r', chdir = TRUE)
library(MASS)
library(xtable)

# 数据生成函数
generate_nb_data <- function(n, beta, X_intercept = TRUE) {
  k <- length(beta) - X_intercept  # 协变量个数
  X <- matrix(rnorm(n * k), nrow = n)  # 随机生成协变量矩阵
  if (X_intercept) {
    X <- cbind(1, X)  # 添加截距项
  }
  log_mu <- X %*% beta  # 计算线性预测器
  mu <- exp(log_mu)  # 使用 logit 函数计算概率 p
  theta <- 2  # 设置过度离散参数，theta > 0
  y <- rnbinom(n, size = theta, mu = mu)  # 生成负二项分布的响应变量
  data <- data.frame(y = y, X = X[, -1, drop = FALSE])  # 将数据集转换为 data frame
  return(data)
}

set.seed(42)
data <- generate_nb_data(n = 2000, beta = c(4, -0.5, 0.2))
fit <- NB_Reg(y ~ ., data = data)
fit2 <- NB_Reg(y ~ ., data = data, type = "bWald")
fit

fit3 <- glm.nb(y ~ ., data = data)
summary(fit3)

formula <- y ~ .
X <- model.matrix(formula, data) # Design matrix
Y <- model.response(model.frame(formula, data)) # Response variable
r_squared(Y,X,fit2$coefficients)
r_squared(Y,X,fit$coefficients)

data$X_random <- rnorm(n = 2000, 10, 1)  # 添加一个随机列
fit4 <- glm.nb(y ~ X.2 + X_random, data = data)
summary(fit4)

fit5 <- NB_Reg(y ~ X.2 + X_random, data = data, type = "bWald")
fit5

compare_nb <- function(beta, n = 2000) {
  set.seed(42)
  data <- generate_nb_data(n, beta)
  # 第一种算法（bWald 方法）
  fit1 <- NB_Reg(y ~ ., data = data, type = "bWald")
  # 第二种算法（glm.nb 方法）
  fit2 <- glm.nb(y ~ ., data = data)
  
  bias1 <- mean(fit1$coefficients - beta)
  bias1_sd <- sd(fit1$coefficients - beta)
  bias2 <- mean(fit2$coefficients - beta)
  bias2_sd <- sd(fit2$coefficients - beta)

  sd1 <- mean(fit1$bwald$SE)  # bWald 的标准误
  sd2 <- mean(summary(fit2)$coefficients[, 2])  # glm.nb 的标准误
  
  rmse1 <- sqrt(mean((fit1$coefficients - beta)^2))
  rmse2 <- sqrt(mean((fit2$coefficients - beta)^2))
  
  results <- list(
    Bias = c(Alg_bWald = bias1, Alg_glmnb = bias2),
    Bias_SD = c(Alg_bWald = bias1_sd, Alg_glmnb = bias2_sd),
    SD = c(Alg_bWald = sd1, Alg_glmnb = sd2),
    RMSE = c(Alg_bWald = rmse1, Alg_glmnb = rmse2)
  )
  return(results)
}

beta = c(-2, 0.2, -0.3, 1)
res <- compare_nb(beta)
xtable(as.data.frame(res), digits = 3)


