DEtype <- function(results, threshold){
  # Invalid input judge
  if(class(results) != "data.frame")
    stop("Invalid input of wrong data type of results")
  if(ncol(results) != 22)
    stop("Invalid input of wrong column number of results")
  if(colnames(results)[21] != "pvalue.adj.FDR" | colnames(results)[16] != "pvalue_LR2" | colnames(results)[17] != "pvalue_LR3")
    stop("Invalid input of wrong column name of results")
  if(class(threshold) != "numeric")
    stop("Invalid input of wrong data type of threshold")
  if(threshold <= 0 | threshold > 0.1)
    stop("Invalid input of wrong range of threshold")

  # Classify the types of DE genes
  results <- cbind(results, NA, NA)
  colnames(results)[c(ncol(results)-1, ncol(results))] <- c("Type", "State")
  for(i in 1:nrow(results)){
    if(results[i,"pvalue.adj.FDR"] < threshold)
    {
      if(results[i,"pvalue_LR2"] < threshold & results[i,"pvalue_LR3"] < threshold){
        results[i,"Type"] <- "DEg"
        if(results[i,"mu_1"] * (1 - results[i,"theta_1"]) >= results[i,"mu_2"] * (1 - results[i,"theta_2"]))
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"pvalue_LR2"] < threshold){
        results[i,"Type"] <- "DEs"
        if(results[i,"theta_1"] <= results[i,"theta_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"pvalue_LR3"] < threshold){
        results[i,"Type"] <- "DEa"
        if(results[i,"mu_1"] >= results[i,"mu_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else{
        results[i,"Type"] <- "DEg"
        if(results[i,"mu_1"] * (1 - results[i,"theta_1"]) >= results[i,"mu_2"] * (1 - results[i,"theta_2"]))
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
    }
    else
      next;
  }
  results
}