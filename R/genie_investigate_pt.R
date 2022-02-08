#' Investigate GENIE BPC Patient
#'
#' @param pt record_id to investigate
#'
#' @return
#'
#' @export
#'
#' @examples
genie_investigate_pt <- function(pt){
  df_pt <- pt_derived %>%
    filter(record_id == pt) %>%
    select(cohort, record_id, redacted)

  df_dx <- ca_dx_derived %>%
    filter(record_id == pt) %>%
    select(cohort, record_id, ca_seq, redcap_ca_index, ca_cadx_int,
           stage_dx, dmets_stage_i_iii)

  df_drugs <- ca_drugs_derived %>%
    filter(record_id == pt) %>%
    select(record_id, ca_seq, regimen_number,
           dob_reg_start_int, pfs_i_g_status, tt_pfs_i_g_days)

  return(list("pt_derived" = df_pt,
              "ca_dx_derived" = df_dx,
              "ca_drugs_derived" = df_drugs))
}
