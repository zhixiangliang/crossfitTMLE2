#'  Estimate Average Treatment Effect (ATE) using from TMLE estimator using cross-fit algorithm (generalization 1)
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param family.y it is the family for outcome model. It can `binomial() (default)` or `"gaussian"`
#' @param learners similar as\code{SL.library()} in `SuperLearner` package.
#' @param control similar as \code{cvControl()} in `SuperLearner` package.
#' @param num_cf number of repetition done. The default is 5.
#' @param n_split number of splits used, default `n_split = 3`
#' @param rand_split logical value; if be FALSE `(default)`, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.
#' @param gbounds value between (0,1) for truncation of predicted probabilities of exposure. The defaults are 5/sqrt(n)log(n) and 1-5/sqrt(n)log(n). See \code{tmle::tmle()} for more information.
#' @param Qbounds used to keep predicted probabilities of outcomes values bounded away from (0,1). The defaults are 5e-04 and 1-5e-04.
#' @param seed numeric value to reproduce the splits distribution
#' @param conf.level confidence limit for confidence interval, `default = 0.95`.
#' @param stat name of the summary statistics, "median" (default) or "mean" to calculate overall ATE from repetitions specific estimates.
#' @return It return a list of two elements. The first element `ATE` is a tibble of the estimates.
#'
#' @import dplyr tibble tidyr purrr furrr
#'
#' @export
#'
#'
#' @examples
#'
#' # See the README file for details
#'
#' sum(1:5)
#'
#'
crossfitTMLE <- function(data,
                         exposure,
                         outcome,
                         covarsT,
                         covarsO,
                         family.y="binomial",
                         learners,
                         control,
                         n_split,
                         num_cf,
                         rand_split=FALSE,
                         gbounds = NULL,
                         Qbounds = 5e-04,
                         seed=146,
                         conf.level=0.95,
                         stat = "median"){

  runs <- list()
  # placeholder_output <- generate_placeholder_output(learners, n_split)
  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)
  n = nrow(data)
  if(is.null(gbounds)) gbounds = 5/sqrt(n)/log(n)
  #################### step 2 and 3 ########################

  for(cf in 1:num_cf){
    seed1 = cf_seed[cf]
    fit_result = try({
        tmle_single(data,
                       exposure,
                       outcome,
                       covarsT,
                       covarsO,
                       family.y,
                       learners,
                       control,
                       n_split,
                       rand_split,
                       gbounds,
                       Qbounds,
                       seed=seed1)}, silent = TRUE)

    if (inherits(fit_result, "try-error")) {
      fit_sngle <-  data.frame(or=NA, log_or_se=NA)
    } else {
      fit_sngle <- fit_result
    }

    runs[[cf]] <- fit_sngle
  }

  res = dplyr::bind_rows(runs)

  # 计算对数OR
  res$log_or <- log(res$or)

  if(stat == "mean"){
    # 在对数尺度上计算均值
    mean_log_or <- mean(res$log_or, na.rm = TRUE)
    # 调整方差：平均的渐近方差 + 交叉拟合间的方差
    # 总方差 = 平均的log_or_se^2 + (log_or - mean_log_or)^2 的均值
    mean_asymptotic_variance <- mean(res$log_or_se^2, na.rm = TRUE)
    between_variance <- mean((res$log_or - mean_log_or)^2, na.rm = TRUE)
    total_variance <- mean_asymptotic_variance + between_variance
    results <- c(mean_log_or, total_variance)
  }
  if(stat == "median"){
    # 在对数尺度上计算中位数
    median_log_or <- median(res$log_or, na.rm = TRUE)
    # 对于中位数，我们如何调整方差？这里我们使用中位数绝对偏差（MAD）吗？
    # 但是原代码中使用了中位数，然后调整方差的方式与均值类似，但使用中位数代替均值。
    # 注意：原代码中对于中位数的调整是：var0 = var + (rd - medians[1])^2，然后取中位数。
    # 这里我们类似地做：
    median_asymptotic_variance <- median(res$log_or_se^2, na.rm = TRUE)
    # 计算每个log_or与中位数的偏差的平方，然后取中位数
    median_between_variance <- median((res$log_or - median_log_or)^2, na.rm = TRUE)
    total_variance <- median_asymptotic_variance + median_between_variance
    results <- c(median_log_or, total_variance)
  }

  # 点估计：取指数
  point_estimate <- exp(results[1])
  # 标准误：sqrt(total_variance)
  se <- sqrt(results[2])

  # 置信区间：在对数尺度上计算，然后取指数
  t.value = qt((1-conf.level)/2, nrow(data), lower.tail = F)
  l_ci_log = results[1] - t.value * se
  u_ci_log = results[1] + t.value * se

  # 转换回OR尺度
  l_ci = exp(l_ci_log)
  u_ci = exp(u_ci_log)

  res1 = tibble(OR=point_estimate, se = se, lower.ci = l_ci, upper.ci = u_ci)

  return(res1)
}
