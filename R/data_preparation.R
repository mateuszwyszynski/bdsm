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
#' @examples
#' df <- data.frame(
#'   year = c(2000, 2001, 2002, 2003, 2004),
#'   country = c("A", "A", "B", "B", "C"),
#'   gdp = c(1, 2, 3, 4, 5),
#'   gdp_lagged = c(NA, 1, 2, 3, 4)
#' )
#'
#' join_lagged_col(df, gdp, gdp_lagged, year, country, 1)
#'
#' @importFrom rlang :=
#'
#' @export
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

#' Perform feature standarization
#'
#' @description
#' This function performs
#' \href{https://en.wikipedia.org/wiki/Feature_scaling}{feature standarization}
#' (also known as z-score normalization), i.e. the features are centered around
#' the mean and scaled with standard deviation.
#'
#' @param df Dataframe with data that should be prepared for LIML estimation
#' @param timestamp_col Column with timestamps (e.g. years)
#' @param entity_col Column with entities (e.g. countries)
#' @param time_effects Whether to introduce time fixed effects
#' (by cross-sectional demeaning)
#' @param scale Whether to divide by the standard deviation \code{TRUE} or not
#' \code{FALSE}. Default is \code{TRUE}.
#'
#' @return A dataframe with standardized features
#'
#' @examples
#' df <- data.frame(
#'   year = c(2000, 2001, 2002, 2003, 2004),
#'   country = c("A", "A", "B", "B", "C"),
#'   gdp = c(1, 2, 3, 4, 5),
#'   ish = c(2, 3, 4, 5, 6),
#'   sed = c(3, 4, 5, 6, 7)
#' )
#'
#' feature_standardization(df, year, country)
#'
#' @export
feature_standardization <- function(df, timestamp_col, entity_col,
                                    time_effects = FALSE, scale = TRUE) {
  if (!time_effects) {
    df %>%
      dplyr::mutate(dplyr::across(!({{ timestamp_col }}:{{ entity_col }}),
                                  function(x) c(scale(x, scale = scale))))
  } else {
    df %>% dplyr::group_by({{ timestamp_col }}) %>%
      dplyr::reframe("{{entity_col}}" := {{ entity_col }},
                     dplyr::across(!{{ entity_col }},
                                   function(x) c(scale(x, scale = scale)))) %>%
      dplyr::arrange({{ entity_col }}) %>% dplyr::ungroup()
  }
}

#' Perform standardization of variables and prepears fixed effects estiamtion
#'
#' @description
#' This function performs
#' \href{https://en.wikipedia.org/wiki/Feature_scaling}{feature standarization}
#' (also known as z-score normalization), i.e. the features are centered around
#' the mean and scaled with standard deviation. Additionally, it allows introduction
#' of cross sectional and time fixed effects through demeaning.
#'
#' @param df Dataframe with data that should be prepared for LIML estimation
#' @param timestamp_col Column with timestamps (e.g. years)
#' @param entity_col Column with entities (e.g. countries)
#' @param standardize Whether to standardize the data (by mean subtraction)
#' @param scale Whether to divide by the standard deviation \code{TRUE} or not
#' \code{FALSE} during standardization. Default is \code{TRUE}
#' @param time_effects Whether to introduce time fixed effects
#' (by cross-sectional demeaning)
#' @param entity_effects Whether to introduce time cross-section effects
#' (by time demeaning)
#'
#'
#' @return A dataframe with standardized variables or/and prepared for fixed effects
#' estimation
#'
#' @examples
#' df <- data.frame(
#'   year = c(2000, 2001, 2002, 2003, 2004),
#'   country = c("A", "A", "B", "B", "C"),
#'   gdp = c(1, 2, 3, 4, 5),
#'   ish = c(2, 3, 4, 5, 6),
#'   sed = c(3, 4, 5, 6, 7)
#' )
#'
#' data_prep(df, year, country, entity_effects = TRUE)
#'
#' @export
data_prep <- function(df, timestamp_col, entity_col,standardize = TRUE,
                      entity_effects = FALSE, time_effects = FALSE,
                      scale = TRUE) {
  if (standardize == TRUE) {
    df %>%
      dplyr::mutate(dplyr::across(!({{ timestamp_col }}:{{ entity_col }}),
                                  function(x) c(scale(x, scale = scale))))
  }
  if (time_effects == TRUE) {
    df %>% dplyr::group_by({{ timestamp_col }}) %>%
      dplyr::reframe("{{entity_col}}" := {{ entity_col }},
                     dplyr::across(!{{ entity_col }},
                                   function(x) c(scale(x, scale = FALSE)))) %>%
      dplyr::arrange({{ entity_col }}) %>% dplyr::ungroup()
  }
  if (entity_effects == TRUE) {
    df %>% dplyr::group_by({{ entity_col }}) %>%
      dplyr::reframe("{{timestamp_col}}" := {{ timestamp_col }},
                     dplyr::across(!{{ timestamp_col }},
                                   function(x) c(scale(x, scale = FALSE)))) %>%
      dplyr::arrange({{ timestamp_col }}) %>% dplyr::ungroup()
  }

}
