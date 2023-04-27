#' Extract module information from kME and colors
#'
#' @param kME data frame of kME correlations of traits with eigentraits
#' @param colors color for modules represented by eigentraits
#'
#' @return data frame with trait, module, kME
#' @export
#' @importFrom dplyr arrange count desc left_join mutate rename tibble
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_remove
#' @importFrom tibble rownames_to_column
#'
module_factors <- function(kME, colors) {
  kME <- 
    dplyr::rename(
      tibble::rownames_to_column(kME),
      trait = rowname)
  
  # Left join with `trait` and `module` to extract kME value.
  object <-
    dplyr::left_join(
      
      # Create tibble with `trait` names (rownames of kME) and `module` colors.
      dplyr::tibble(
        trait = kME$trait,
        module = colors),
      
      # Mutate `module` by removing leading "kME".
      dplyr::mutate(
        # Pivot kME to long form 
        tidyr::pivot_longer(
          kME,
          -trait,
          names_to = "module", values_to = "kME"),
        module = stringr::str_remove(.data$module, "^kME")),
      
      by = c("trait", "module"))
  
  # Arrange `module` colors in descending order by count.
  modcolors <- dplyr::arrange(
    dplyr::count(
      object,
      .data$module),
    dplyr::desc(n))$module
  
  object$module <- factor(object$module, modcolors)
  
  object
}