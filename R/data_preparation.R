#' Dataframe with no lagged column
#'
#' This function allows to turn data in the format with lagged values for a
#' chosen column (i.e. there are two columns with the same quantity, but one
#' column is lagged in time) into the format with just one column
#'
#' @param df Dataframe with data with a column with lagged values
#' @param col Column with quantity not lagged
#' @param col_lagged Column with the same quantity as \code{col}, but the values
#' are lagged in time
#' @param timestamp_col Column with timestamps (e.g. years)
#' @param entity_col Column with entities (e.g. countries)
#' @param timestep Difference between timestamps (e.g. 10)
#'
#' @return
#' A dataframe with two columns merged, i.e. just one column with the desired
#' quantity is left.
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
join_lagged_col <- function(df, col, col_lagged, timestamp_col,
                            entity_col, timestep = NULL) {
  non_lagged_df <- df %>%
    dplyr::select(c({{ timestamp_col }}, {{ entity_col }}, {{ col }}))
  lagged_df <- df %>%
    dplyr::select(c({{ timestamp_col }}, {{ entity_col }}, {{ col_lagged }})) %>%
    dplyr::mutate("{{timestamp_col}}" := {{ timestamp_col }} - timestep)

  lagged_df %>% dplyr::full_join(non_lagged_df,
                          by = dplyr::join_by(
                            {{ timestamp_col }} == {{ timestamp_col }},
                            {{ entity_col }} == {{ entity_col }}
                          )) %>%
    dplyr::arrange({{ entity_col }}, {{ timestamp_col }}) %>%
    dplyr::mutate("{{col}}" := dplyr::coalesce({{ col }}, {{ col_lagged }})) %>%
    dplyr::select(!{{ col_lagged }}) %>%
    dplyr::left_join(dplyr::select(df, !{{ col }} & !{{ col_lagged }}),
              by = dplyr::join_by(
                {{ timestamp_col }} == {{ timestamp_col }},
                {{ entity_col }} == {{ entity_col }}
              ))
}
