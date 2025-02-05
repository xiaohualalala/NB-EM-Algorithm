# 安装并加载必要的包
if (!require("numDeriv")) install.packages("numDeriv", dependencies = TRUE)
library(numDeriv)

# 定义Newton-Raphson算法估计负二项分布参数的函数
newton_raphson <- function(data, starting_values = NULL, tol = 1e-5, max_iter = 50) {
    start_time <- Sys.time() # 开始时间

    # 定义负二项分布的对数似然函数
    log_likelihood <- function(params) {
        size <- params[1]
        prob <- params[2]
        if (size <= 0 || prob <= 0 || prob >= 1) {
            return(-Inf)
        }
        sum(dnbinom(data, size = size, prob = prob, log = TRUE))
    }

    if (is.null(starting_values)) {
        avg <- mean(data)
        var <- var(data)
        starting_values[1] <- avg * avg / (var - avg)
        starting_values[2] <- avg / var
    }

    params <- starting_values
    converged <- FALSE
    for (i in 1:max_iter) {
        grad <- grad(log_likelihood, params)
        hess <- hessian(log_likelihood, params)

        # 解决Hessian矩阵的数值问题
        if (any(is.na(hess)) || any(is.infinite(hess))) {
            warning("Hessian matrix contains NA or infinite values. Iteration stopped.")
            break
        }

        params_new <- params - solve(hess) %*% grad

        # 检查新参数是否在合理范围内
        if (params_new[1] <= 0 || params_new[2] <= 0 || params_new[2] >= 1) {
            warning("Parameters out of bounds. Iteration stopped.")
            break
        }

        if (sum(abs(params_new - params)) < tol) {
            params <- params_new
            converged <- TRUE
            break
        }

        params <- params_new
    }

    end_time <- Sys.time() # 结束时间
    duration <- difftime(end_time, start_time, units = "secs") # 计算时间

    return(list(params = params, iterations = i, converged = converged, duration = duration))
}

