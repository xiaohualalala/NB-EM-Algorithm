library(Rcpp)
library(ggplot2)
require(pscl)
library(xtable)
library(gridExtra)

setwd("/Users/guosa/Desktop/日常/科研/负二项分布/code/NB-EM-Algorithm/Simulation")

sourceCpp("../Simulation/Samples/Generate_Samples.cpp")
sourceCpp("../Simulation/NB/SimulateNB.cpp")
sourceCpp("../Simulation/MixNB/SimulateMixNB.cpp")
sourceCpp("../Simulation/ZINB/SimulateZINB.cpp")
sourceCpp("../MLE/C/ZINB_MLE.cpp")
sourceCpp("../MLE/C/NB_MLE.cpp")

# 设置随机数种子
set.seed(123)

# NB
# 模拟参数
prob_values <- seq(0.2, 0.8, 0.2)
size_values <- seq(1, 9, 2)
res_NB <- SimulateNBs(prob_values, size_values, 200, 2000)

res_NB_latex <- xtable(res_NB)
# 迭代次数与参数的热力图
p1 <- 
ggplot(res_NB, aes(x = as.factor(size), y = as.factor(prob))) +
    geom_tile(aes(fill = avg_t)) +
    scale_fill_gradient(low = "#ffeda0", high = "#f03b20") +
    labs(x = "size", y = "prob", fill = "Iterations") +
    theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/hotplot.pdf", p1)

# 不同初始点对收敛次数的影响
source('../MLE/R/NB_MLE.r', chdir = TRUE)

# 生成模拟数据
set.seed(123)
X <- rnbinom(2000, size = 5, prob = 0.5)

# random
random_conver <- logical(200)
random_time <- numeric(200)
random_iter <- numeric(200)

for (i in 1:length(conver)) {
    init_values <- c(runif(1, 0.1, 50), runif(1, 0.01, 0.99))
    result <- NB_MLE(X, init_values = init_values, verbose = FALSE)
    random_conver[i] <- result$converged
    random_time[i] <- result$duration
    random_iter[i] <- result$iter
}

random_rate <- sum(random_conver) / length(random_conver)
random_cpu <- mean(random_time)
sd(random_time)
mean(random_iter)
sd(random_iter)

# extreme 

extreme_conver <- logical(4)
extreme_time <- numeric(4)
extreme_iter <- numeric(4)

init_values <- list(c(0.1,0.001),c(0.1,0.999),c(1000,0.001),c(1000,0.999))

for (i in 1:length(init_values)) {
    result <- NB_MLE(X, init_values = unlist(init_values[i]), verbose = FALSE)
    extreme_conver[i] <- result$converged
    extreme_time[i] <- result$duration
    extreme_iter[i] <- result$iter
}
extreme_rate <- sum(extreme_conver) / length(extreme_conver)
extreme_cpu <- mean(extreme_time)
sd(extreme_time)
mean(extreme_iter)
sd(extreme_iter)

# moment 
result <- NB_MLE(X, verbose = FALSE)
sum(result$converged)
result$duration
result$iter

library(plotly)

# 定义初始值的范围
prob_values <- seq(0.1, 0.9, length.out = 40)  # 取10个概率值
size_values <- seq(1, 20, length.out = 40)       # 取10个size值

init_values <- c(runif(1, 0.1, 0.9), runif(1, 1, 20))
# 记录结果
results <- expand.grid(Prob = prob_values, Size = size_values)
results$Iterations <- NA

for (i in 1:nrow(results)) {
    init_values <- c(results$Size[i], results$Prob[i])
    result <- NB_MLE(X, init_values = init_values, verbose = FALSE)
    results$Iterations[i] <- result$iter
}

# 确保数据格式正确
Iterations <- matrix(results$Iterations, nrow = length(prob_values), ncol = length(size_values))

# 创建三维曲面图
p2 <- 
plot_ly(
    x = ~prob_values,
    y = ~size_values,
    z = ~Iterations,
    type = 'surface'
) %>%
layout(
    scene = list(
        xaxis = list(title = 'Prob'),
        yaxis = list(title = 'Size'),
        zaxis = list(title = 'Iterations')
    )
)
library(htmlwidgets)
saveWidget(p2, "/Users/guosa/Desktop/毕业论文/figures/初始值与迭代次数.html")

## MixNB
source("../Simulation/MixNB/VisualMixNB.r")

# 双样本
# 模拟参数
weights <- c(0.4, 0.6)
means <- c(100, 200)
vars <- c(250, 250)
probs <- means / vars
sizes <- means^2 / (vars - means)

set.seed(123)
probs <- c(0.25, 0.5)
sizes <- c(20, 150)

# 生成数据并进行分组
res2 <- SimulateMixNB(probs, sizes, weights, 2000)

# 对比图
p3 <- plotHistDensity(res$observations, res$estimation, res$accuracy)
ggsave("/Users/guosa/Desktop/毕业论文/figures/mix2.pdf", p3)

# 运行100次，观察准确率与迭代次数
data <- SimulateMixNBs(probs, sizes, weights, 100, 2000)

df_accuracy <- data.frame(Value = data$accuracy, Type = "Accuracy")
df_iterations <- data.frame(Value = data$iterations, Type = "Iterations")
mean(df_iterations$Value)
sd(df_iterations$Value)

K <- length(weights)
df_size <- data$eSize
size_mean <- sapply(1:K, function(i) mean(sapply(df_size, `[`, i)))
size_sd <- sapply(1:K, function(i) sd(sapply(df_size, `[`, i)))

df_prob <- data$eProb
prob_mean <- sapply(1:K, function(i) mean(sapply(df_prob, `[`, i)))
prob_sd <- sapply(1:K, function(i) sd(sapply(df_prob, `[`, i)))

df_weight <- data$eWeight
weight_mean <- sapply(1:K, function(i) mean(sapply(df_weight, `[`, i)))
weight_sd <- sapply(1:K, function(i) sd(sapply(df_weight, `[`, i)))

# accuracy箱线图
p5 <- 
ggplot(df_accuracy, aes(x = Type, y = Value)) +
    geom_boxplot() +
    theme_minimal() +
    ylim(0.965, 0.975) + 
    annotate("text", x = Inf, y = Inf, label = paste("Accuracy:", round(mean(df_accuracy$Value),4),"\u00B1", round(sd(df_accuracy$Value),4)), hjust = 1.1, vjust = 2, size = 8, color = "blue") + 
    labs(x = "", y = "Accuracy")

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix4.pdf", p5)

# iterations箱线图
p6 <- 
ggplot(df_iterations, aes(x = Type, y = Value)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "", y = "Iterations")

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix5.pdf", p6)

# 模拟参数
weights <- c(0.2, 0.5, 0.3)
means <- c(100, 200, 300)
vars <- c(1000, 1000, 1000)
probs <- means / vars
sizes <- means^2 / (vars - means)

probs <- c(0.1, 0.2, 0.3)
sizes <- c(10, 50, 150)

# 生成数据并进行分组
res3 <- SimulateMixNB(probs, sizes, weights, 2000)

# 对比图
p4 <- plotHistDensity(res3$observations, res3$estimation, res3$accuracy)

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix3.pdf", p4)

# 运行100次，观察准确率与迭代次数
data <- SimulateMixNBs(probs, sizes, weights, 100, 2000)

df_accuracy <- data.frame(Value = data$accuracy, Type = "Accuracy")
mean(df_accuracy$Value)
sd(df_accuracy$Value)
df_iterations <- data.frame(Value = data$iterations, Type = "Iterations")
mean(df_iterations$Value)
sd(df_iterations$Value)

K <- length(weights)
df_size <- data$eSize
size_mean <- sapply(1:K, function(i) mean(sapply(df_size, `[`, i)))
size_sd <- sapply(1:K, function(i) sd(sapply(df_size, `[`, i)))

df_prob <- data$eProb
prob_mean <- sapply(1:K, function(i) mean(sapply(df_prob, `[`, i)))
prob_sd <- sapply(1:K, function(i) sd(sapply(df_prob, `[`, i)))

df_weight <- data$eWeight
weight_mean <- sapply(1:K, function(i) mean(sapply(df_weight, `[`, i)))
weight_sd <- sapply(1:K, function(i) sd(sapply(df_weight, `[`, i)))

p5 <- grid.arrange(p3, p4, nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/mix6.pdf", p5)

# confusion matrix
conf_matrix <- confusionMatrix(as.factor(res$predicted_labels), as.factor(res$true_labels))
conf_matrix

# roc
library(pROC)
library(caret)

roc_curve2 <- roc(res2$true_labels, res2$prediction, plot = TRUE, print.auc = TRUE)
roc_curve3 <- roc(res3$true_labels, res3$prediction, plot = TRUE, print.auc = TRUE)
pdf("/Users/guosa/Desktop/毕业论文/figures/ROC_Curve.pdf")
plot(roc_curve2, col = "coral", main = "ROC Curves", print.auc = TRUE, auc.polygon = FALSE)
plot(roc_curve3, col = "navy", add = TRUE, print.auc = TRUE, auc.polygon = FALSE, print.auc.y = 0.4)
legend("bottomright", legend = c("K=2", "K=3"), col = c("coral", "navy"), lwd = 2)
dev.off()

## ZINB
# 模拟参数
prob_values <- seq(0.2, 0.6, 0.2)
size_values <- seq(6, 10, 2)
zero_values <- seq(0.2, 0.6, 0.2)

res_ZINB <- SimulateZINBs(prob_values, size_values, zero_values, 20, 2000)
res_ZINB_latex <- xtable(res_ZINB)

# 分类效果模拟
set.seed(123)
prob <- 0.1
size <- 2
zero <- 0.2

ZINB <- SamplesZINB(prob, size, zero, 1, 2000)
X <- unlist(ZINB$samples)
true_labels <- unlist(ZINB$labels)

res_ZINB <- ZINB_MLE(X)
eZero <- res_ZINB$Parameter_Estimates$eZero
eSize <- res_ZINB$Parameter_Estimates$eSize
eProb <- res_ZINB$Parameter_Estimates$eProb
prob_nb <- (1 - eZero) * dnbinom(X, size = eSize, prob = eProb)
prob_zero <- eZero * (X == 0)
tau <- prob_zero / (prob_zero + prob_nb)
predict_labels <- ifelse(tau > 0.5, 0, 1)

eZero / ((1 - eZero) * dnbinom(0, size = eSize, prob = eProb) + eZero)
sum(X==0) dnbinom(0, size = eSize, prob = eProb)*2000

accuracy <- mean(true_labels == predict_labels)

NB_part <- X[predict_labels == 1]
Zero_part <- X[predict_labels == 0]
Source <- ifelse(tau > 0.5, "Negative Binomial", "Zero-Inflated")

# 计算估计的负二项分布的概率密度
x_values <- seq(0, max(NB_part), by = 1)
nb_density <- dnbinom(x_values, size = eSize, prob = eProb)

# 创建数据框用于绘制概率密度曲线
density_data <- data.frame(x_values, nb_density)

data <- data.frame(Value = X)

# 绘制直方图
p7 <- 
    ggplot(data, aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = diff(range(X)) / 70, fill = "skyblue", color = "white") +
    labs(x = "Value",
    y = "Density") +
    theme_minimal() 

ggsave("/Users/guosa/Desktop/毕业论文/figures/zinb-hist.pdf", p7)

  +
  annotate("text", x = Inf, y = Inf, label = paste("Accuracy:", accuracy), hjust = 1.1, vjust = 2, size = 5, color = "blue")


# 添加概率密度曲线
p8 <- 
p7 + geom_line(data = density_data, aes(x = x_values, y = nb_density), color = "red", size = 1)

ggsave("/Users/guosa/Desktop/毕业论文/figures/zinb1.pdf", p8)

# Rate 
set.seed(123)

probs <- seq(0.2, 0.8, 0.2)
sizes <- seq(1, 10, 1)
rates <- SimulateRates(probs,sizes,20,2000)

p9 <- 
ggplot(rates, aes(x = size, y = avg_rate, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_rate - sd_rate, ymax = avg_rate + sd_rate), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Size", 
       y = "Convergence Rate", 
       color = "Prob", 
       title = "Convergence Rate vs Size") +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate1.pdf", p9)

probs <- seq(0.1, 0.9, 0.1)
sizes <- seq(2, 8, 2)
rates <- SimulateRates(probs,sizes,20,2000)
 
p10 <- 
ggplot(rates, aes(x = prob, y = avg_rate, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_rate - sd_rate, ymax = avg_rate + sd_rate), width = 0.01, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Prob", y = "Convergence Rate", color = "Size", title = "Convergence Rate vs Prob") +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate2.pdf", p10)

p11 <- grid.arrange(p9, p10, nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/rates.pdf", p11)

prob <- 0.6
size <- 1
sample_sizes <- seq(500,100000,500)
rates <- data.frame(n = sample_sizes, rate = numeric(200))
for (i in 1:200){
    rates$rate[i] <- SimulateRate(prob,size,1,sample_sizes[i])$conv_rates
}

p12 <- 
ggplot(rates, aes(x = n, y = rate)) +
  geom_point(color='red', size=0.8) +
  labs(x = "n", y = "Convergence Rate") +
  theme_minimal()+
  ylim(0, 1) +
  annotate("text", x = max(rates$n), y = 1, label = "Size=1, Prob=0.6", hjust=1, vjust=1, color = "blue", size = 7)

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate3.pdf", p12)

sourceCpp("../Simulation/NB/SimulateNB.cpp")

set.seed(123)

probs <- seq(0.2, 0.8, 0.2)
sizes <- seq(1, 10, 1)
errors <- SimulateErrors(probs,sizes,20,2000)

p13 <- 
ggplot(errors, aes(x = size, y = avg_error_2, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.2) + 
  labs(x = "Size", 
       y = bquote("Standard Error (" ~ phi[1] ~ ")"), 
       color = "Prob", 
       title = expression("Standard Error (" ~ phi[1] ~ ") vs Size" )) +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error1.pdf", p13)

p14 <- 
ggplot(errors, aes(x = size, y = avg_error_2, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.5) + 
  labs(x = "Size", 
       y = bquote("Standard Error (" ~ phi[2] ~ ")"), 
       color = "Prob",
       title = expression("Standard Error (" ~ phi[2] ~ ") vs Size" )) +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error2.pdf", p14)

probs <- seq(0.1, 0.9, 0.1)
sizes <- seq(2, 8, 2)
errors <- SimulateErrors(probs,sizes,20,2000)

p15 <- 
ggplot(errors, aes(x = prob, y = avg_error_2, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.01,  linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.2) + 
  labs(x = "Prob", 
       y = bquote("Standard Error (" ~ phi[1] ~ ")"), 
       color = "Size",
       title = expression("Standard Error (" ~ phi[1] ~ ") vs Prob" )) +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error3.pdf", p15)

p16 <- 
ggplot(errors, aes(x = prob, y = avg_error_2, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.01, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.5) + 
  labs(x = "Prob", 
       y = bquote("Standard Error (" ~ phi[2] ~ ")"), 
       color = "Size",
       title = expression("Standard Error (" ~ phi[2] ~ ") vs Prob" )) +
  theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error4.pdf", p16)

p17 <- grid.arrange(p13,p15,nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/error_2.pdf", p17)

p18 <- grid.arrange(p14,p16,nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/error_2.pdf", p18)







# EM与BFGS比较

library(bbmle)

# 生成一些负二项分布的样本数据
set.seed(123)
data <- rnbinom(100, size = 5, prob = 0.1)

# 定义负二项分布的负对数似然函数
nll <- function(size, prob) {
  -sum(dnbinom(data, size = size, prob = prob, log = TRUE))
}

# 使用mle2函数进行参数估计
fit <- mle2(nll, start = list(size = 1, prob = 0.1))
summary(fit)$Coefficients


compare <- function(prob, size, zero, sample_size) {
    sample <- unlist(SamplesZINB(prob, size, zero, 1, sample_size))

    res_BFGS <- zeroinfl(count ~ 1, data = data.frame(count = sample), dist = "negbin")
    size_BFGS <- res_BFGS$theta
    mu_BFGS <- exp(res_BFGS$coefficients$count[["(Intercept)"]])
    zero_BFGS <- exp(res_BFGS$coefficients$zero[["(Intercept)"]]) / (1 + exp(res_BFGS$coefficients$zero[["(Intercept)"]]))
    prob_BFGS <- size_BFGS / (size_BFGS + mu_BFGS)

    res_EM <- ZINB_MLE(sample)

    df <- data.frame(
        prob = prob,
        EM_prob = res_EM$Parameter_Estimates$eProb,
        BFGS_prob = prob_BFGS,
        size = size,
        EM_size = res_EM$Parameter_Estimates$eSize,
        BFGS_size = size_BFGS,
        zero = zero,
        EM_zero = res_EM$Parameter_Estimates$eZero,
        BFGS_zero = zero_BFGS
    )

    return(df)
}

probs <- seq(0.2, 0.4, 0.2)
sizes <- seq(6, 10, 2)
zeros <- seq(0.1, 0.2, 0.1)
parameters <- expand.grid(prob = probs, size = sizes, zero = zeros)
EM_BFGS = data.frame()
for (i in 1:nrow(parameters)) {
    df <- compare(parameters$prob[i], parameters$size[i], parameters$zero[i], sample_size = 1000)
    EM_BFGS <- rbind(EM_BFGS, df)
}
EM_BFGS





# gene data
library(readxl)
gene_data <- read_excel("../Simulation/ZINB/gene.xlsx")
head(gene_data)
# 行为基因，列为样本

# 计算每个基因在所有样本中的平均表达量
average_expression_before <- rowMeans(gene_data[, -1])

# lop1p=log(1+x)
ggplot(data.frame(Expression = log1p(average_expression_before)), aes(x = Expression)) +
    geom_histogram(bins = 80, fill = "skyblue", color = "black") +
    theme_minimal() +
    ggtitle("Histogram of Average Log-transformed Expression") +
    xlab("Log-transformed Average Expression") +
    ylab("Frequency")

# 直方图中接近零的表达量区域高频:存在大量低表达或未表达的基因
# 直方图的尾部延伸较远，且尾部的条形较低但持续存在:有一些基因具有较高的表达水平
# 数据分布
# 峰值靠近左侧:大多数基因表达量较低

# 过滤掉低表达基因
keep <- rowSums(gene_data[, -1] >= 5) >= 2 # 至少在两个样本中有5个以上的计数
filtered_data <- gene_data[keep, ]

# 输出过滤后的结果
cat("Number of genes before filtering:", nrow(gene_data), "\n")
cat("Number of genes after filtering:", nrow(filtered_data), "\n")

average_expression_after <- rowMeans(filtered_data[, -1])

# 绘制过滤后直方图
ggplot(data.frame(Expression = log1p(average_expression_after)), aes(x = Expression)) +
    geom_histogram(bins = 80, fill = "skyblue", color = "black") +
    theme_minimal() +
    ggtitle("Histogram of Average Log-transformed Expression After Filtering") +
    xlab("Log-transformed Average Expression") +
    ylab("Frequency")


# 标准化
# 若已知样本信息和实验信息
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.16")
BiocManager::install("DESeq2")

# CPM标准化
# 计算每个样本的总计数
total_counts <- colSums(filtered_data[, -1])

# 用每个样本的总计数除以每个表达值，然后乘以10^6进行标准化
cpm_data <- t(t(filtered_data[, -1]) / total_counts * 10^6)
head(cpm_data)

# 进行对数转换以使数据更加正态分布
log_cpm_data <- log1p(cpm_data)

# 基因数量很大，随机选择一些基因
set.seed(123)
sample <- log_cpm_data[sample(1:nrow(log_cpm_data), 10), ]

# 主成分分析
pca_res <- prcomp(t(log_cpm_data), scale. = TRUE)
plot(pca_res$x[, 1:2], asp = 1, xlab = "PC1", ylab = "PC2", main = "PCA of High-variance Genes")

# 假设gene_data是你的数据框架，且第一列为GeneID
results_list <- list() # 用于存储每个基因的拟合结果
gene_idx <- 2
for (gene_idx in 1:100) {
    # 提取该基因的表达数据
    X <- as.numeric(filtered_data[gene_idx, -1]) # 排除GeneID列

    # 调用ZINB_MLE函数拟合模型
    result <- ZINB_MLE(X)$Parameter_Estimates

    # 将结果存储在列表中
    results_list[[gene_idx]] <- result
}

# 检查前几个基因的结果
results_list[1] # 假设你想查看第一个基因的结果





# ZINB回归模型对比
zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")
zinb <- within(zinb, {
    nofish <- factor(nofish)
    livebait <- factor(livebait)
    camper <- factor(camper)
})

summary(zinb)

m1 <- zeroinfl(count ~ 1,
    data = zinb, dist = "negbin"
)
m1$coefficients
summary(m1)

source("../Reg/ZINB_Reg.r")

ZINB_EM_IWLS(count ~ child, data = zinb)$coefficients

m <-
    NB <- SamplesNB(0.9, 5, 10, 1000)
ZINB <- SamplesZINB(0.5, 1.0, 0.2, 10, 1000)
probs <- c(0.5, 0.3)
sizes <- c(4, 8)
weights <- c(0.6, 0.4)
zero <- c(0.2, 0.4)
MixNB <- SamplesMixNB(probs, sizes, weights, 1, 2000)
MixZINB <- SamplesMixZINB(probs, sizes, zero, weights, 1, 2000)






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
