#' Simulate sample metadata for two organisms
#'
#' @importFrom magrittr %>%
#' @param org Character: organism label ("EC" or "PP").
#' @param times Character vector of time points.
#' @param groups Character vector of group labels.
#' @param batches Character vector of batch labels.
#' @param reps Integer: number of replicates per condition.
#' @return A data.frame of sample metadata.
#' @importFrom dplyr arrange mutate
#' @importFrom tidyr expand_grid
#' @export
simulate_samples <- function(org,
                             times   = c("0h", "2h", "4h", "8h"),
                             groups  = if (org == "EC") c("less", "equal", "more", "all")
                             else c("less", "equal", "more", "none"),
                             batches = paste0("B", 1:2),
                             reps    = 4) {
  tidyr::expand_grid(
    organism = org,
    time     = times,
    group    = groups,
    batch    = batches,
    rep      = seq_len(reps),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(batch, time, group, rep) %>%
    dplyr::mutate(
      sample_id = paste(organism, group, time, batch, rep, sep = "_"),
      ratio0    = group
    )
}
