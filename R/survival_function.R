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
#' survival_function(df = gtsummary::trial,
#' time = ttdeath, status = death)
survival_function <- function(df, time1 = NULL, time, status,
                              covar = NULL,
                              time_units = "Months",
                              time_origin = "Diagnosis"){

  #### km
  # if left truncation, no covariate
  if (!is.null(time1) & is.null(covar)){
    message("Left truncation, no covariate")
    km <- survminer::surv_fit(survival::Surv(
      time = get(time1),
      time2 = get(time),
      get(status)
    ) ~ NULL,
    data = df
    )
    # left truncation, with covariate
  } else if (!is.null(time1) & !is.null(covar)) {
    message("Left truncation, with covariate")
    km <- survminer::surv_fit(survival::Surv(
      time = get(time1),
      time2 = get(time),
      get(status)
    ) ~ get(covar),
    data = df
    )
    # no left truncation, no covariate
  } else if (is.null(time1) & is.null(covar)) {
    message("No left truncation, no covariate")
    km <- survminer::surv_fit(survival::Surv(
      get(time),
      get(status)
    ) ~ NULL,
    data = df
    )
    # no left truncation, with covariate
  } else if (is.null(time1) & !is.null(covar)) {
    message("No left truncation, with covariate")
    # without left truncation
    km <- survminer::surv_fit(survival::Surv(
      get(time),
      get(status)
    ) ~ get(covar),
    data = df
    )
  }

  #### table of median survival
  # table and legend labels
  # if strata
  if (!is.null(covar)){
    # tbl_covar_lab = {{covar}} # deprecated and no longer accepts a single string
    tbl_covar_lab = eval(parse(text = paste0("`get(covar)` ~ \'", covar, "\'")))
    labs = stringr::word(names(km$strata), 2, sep = "=")
  } else { # if no strata
    tbl_covar_lab = NULL
    labs = NULL
  }

  tbl_median_surv <- gtsummary::tbl_survfit(km,
                                 probs = 0.5,
                                 label = tbl_covar_lab,
                                 label_header =
                                   paste0("**Median survival (95% CI; ", time_units, ")**"))

  #### plot
  # if months vs years
  if (stringr::str_to_upper(time_units) == "MONTHS"){
    break_x = 12
  } else if (stringr::str_to_upper(time_units) == "YEARS"){
    break_x = 2
  }

  plot_survival <- survminer::ggsurvplot(
    fit = km,
    data = df,
    legend.title = "",
    legend.labs = labs,
    risk.table = TRUE,
    tables.y.text = FALSE,
    break.x.by = break_x,
    ylab = "Survival Probability",
    xlab = paste0("Time (", time_units, ") From ", time_origin)
  )

  # return data frames, tables and survival analyses
  return(list("km" = km,
              "tbl_median_surv" = tbl_median_surv,
              "plot_survival" = plot_survival))
}
