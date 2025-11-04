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
      fit_sngle <-  data.frame(rd=NA, var = NA , S1 = NA, S0 = NA)
    } else {
      fit_sngle <- fit_result
    }

    runs[[cf]] <- fit_sngle
  }

  res = dplyr::bind_rows(runs)

  if(stat == "mean"){
    medians <- apply(res, 2, mean, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    results <- apply(res, 2, mean, na.rm = TRUE)
  }
  if(stat == "median"){
    medians <- apply(res, 2, median, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    results <- apply(res, 2, median, na.rm = TRUE)
  }
  t.value = qt((1-conf.level)/2, nrow(data), lower.tail = F)

  l_ci = results[1] - t.value*sqrt(results[3])
  u_ci = results[1] + t.value*sqrt(results[3])

  res1 = tibble(ATE=results[1], se = sqrt(results[3]), lower.ci = l_ci, upper.ci = u_ci, S1 = results[3] , S0 = results[4])

  return(res1)
}
