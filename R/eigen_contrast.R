#' Eigen Contrasts
#'
#' @param object list object with trait module information
#' @param contr_object data frame of class `conditionContrasts`
#'
#' @return data frame
#' @export
#' @importFrom dplyr bind_rows mutate
#' @importFrom purrr map transpose
#' @importFrom tidyr pivot_longer separate_wider_delim
#' @importFrom tibble rownames_to_column
#'
eigen_contrast <- function(object, contr_object) {
  object <- dplyr::mutate(
    dplyr::bind_rows(
      purrr::map(purrr::transpose(object)$eigen, eigen_df)),
    trait = factor(.data$trait, unique(.data$trait)),
    p.value = 10 ^ -match(.data$trait, rev(levels(.data$trait))))
  
  class(object) <- c("conditionContrasts", class(object))
  attr(object, "conditions") <- attr(contr_object, "conditions")
  attr(object, "termname") <- attr(contr_object, "termname")
  object
}
eigen_df <- function(object) {
  tidyr::separate_wider_delim(
    tidyr::pivot_longer(
      tibble::rownames_to_column(object, "dataset_strain"),
      -dataset_strain,
      names_to = "trait", values_to = "value"),
    dataset_strain, delim = "_",
    names = c("dataset", "strain","sex"))
}
