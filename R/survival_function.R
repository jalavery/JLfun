#' Survival Analyes
#'
#' A function for running common pieces of a survival analyses: (1) the KM
#' model (2) the median survival (3) the KM curve
#'
#' @param df The data frame containing time, status and covariate variables
#' @param time1 If the analysis involves left truncation, this is the time
#' from the origin to entry into the risk set
#' @param time The time from the origin to the event or end of follow-up
#' @param covar The covariate of interest
#' @param time_units An optional character string indicating the units of time.
#' The default is "Months". This is used to properly label the units for median
#' survival and the x-axis on the survival plot.
#' @param time_origin An optional character string indicating the time origin.
#' The default is "Diagnosis". This is used to label the x-axis on the survival
#' plot.
#'
#' @return
#' @export
#'
#' @examples
survival_function <- function(df, time1 = NULL, time, status,
                              covar,
                              time_units = "Months",
                              time_origin = "Diagnosis",
                              ...){

  # if no covariate, just include 1
  if (missing(covar)){
    covar <- rep(1, nrow(df))
  }

  # km
  # if left truncation
  if (!is.null(time1)){
    km <- survminer::surv_fit(Surv(
      time = get(time1),
      time2 = get(time),
      get(status)
    ) ~ covar,
    data = df
    )
  } else {
    km <- survminer::surv_fit(Surv(
      get(time),
      get(status)
    ) ~ get(covar),
    data = df
    )
  }

  # table of median survival
  tbl_median_surv <- gtsummary::tbl_survfit(km,
                                 probs = 0.5,
                                 label_header =
                                   paste0("**Median survival (95% CI; ", time_units, ")**"))

  # plot
  plot_survival <- survminer::ggsurvplot(
    fit = km,
    data = df,
    legend.title = "",
    legend.labs = stringr::word(names(km$strata), 2, sep = "="),
    risk.table = TRUE,
    tables.y.text = FALSE,
    break.x.by = 12,
    ylab = "Survival Probability",
    xlab = paste0("Time (", time_units, ") From ", time_origin)
  )

  # return data frames, tables and survival analyses
  return(list("km" = km,
              "tbl_median_surv" = tbl_median_surv,
              "plot_survival" = plot_survival))
}
