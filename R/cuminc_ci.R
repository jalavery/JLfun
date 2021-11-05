#' Compute confidence intervals for cumulative incidence
#'
#' @param ci_fit The object returned by
#' @param times A vector of times at which to estimate the confidence interval
#' @param conf_level The confidence interval; default is 0.95.
#'
#' @return
#' @export
#'
#' @examples See Tringale-Second Cancer project
cuminc_ci <- function(ci_fit, times, conf_level = 0.95) {
  # get estimates at each time point
  est <- cmprsk::timepoints(ci_fit, times)$est %>%
    as_tibble() %>%
    rownames_to_column("event") %>%
    filter(event == 1) %>%
    pivot_longer(cols = -event,
                 names_to = "time",
                 values_to = "est")

  # get variance at each time point
  var <- cmprsk::timepoints(ci_fit, times)$var %>%
    as_tibble() %>%
    rownames_to_column("event") %>%
    filter(event == 1) %>%
    pivot_longer(cols = -event,
                 names_to = "time",
                 values_to = "var")

  # get z based on % CI
  z <- qnorm(1-(1-conf_level)/2)

  # table of estimate + CI
  full_join(est, var,
            by = c("event", "time")) %>%
    mutate(lower = round(est ^ exp(-z*sqrt(var)/(est*log(est))), 3),
           upper = round(est ^ exp(z*sqrt(var)/(est*log(est))), 3),
           est_ci = paste0(round(est, 3), " (", lower, ", ", upper, ")"))
}
