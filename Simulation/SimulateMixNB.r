# 定义混合负二项分布的密度函数
mixNB_density <- function(x, sizes, probs, weights) {
  # 初始化密度向量
  density <- numeric(length(x))
  # 计算混合密度
  for (i in 1:length(x)) {
    for (j in 1:length(weights)) {
      density[i] <- density[i] + weights[j] * dnbinom(x[i], size = sizes[j], prob = probs[j])
    }
  }
  return(density)
}

# 分析混合负二项分布数据的函数
SimulateMixNB <- function(probs, sizes, weights, num_samples, sample_size) {
  # 生成混合负二项分布样本
  MixNB <- SamplesMixNB(probs, sizes, weights, num_samples, sample_size)
  observations <- MixNB$samples
  true_labels <- MixNB$labels
  K <- length(weights)
  
  predicted_labels <- list()
  accuracy <- numeric(length(observations))
  
  # 对每组样本进行分析
  for (i in 1:length(observations)) {
    observation <- unlist(observations[i])
    true_label <- unlist(true_labels[i])
    
    # 使用MLE估计参数
    res <- mixNB_MLE(observation, K, verbose = FALSE)
    
    # 预测标签
    predicted_label <- numeric(length(observation))
    for (j in 1:length(observation)) {
      probability <- sapply(1:K, function(k) {
        dnbinom(observation[j], size = res$eSize[k], prob = res$eProb[k])
      })
      predicted_label[j] <- which.max(probability)
    }
    
    # 计算准确率
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
  
  return(list(accuracy = accuracy, true_labels = true_labels, predicted_labels = predicted_labels))
}