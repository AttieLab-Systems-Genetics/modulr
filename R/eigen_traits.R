#' Compare Eigen Traits with Original Traits
#'
#' @param object module list object
#' @param sexname name of sex
#' @param modulename 
#' @param contr_object 
#'
#' @return data frame
#' @export
#' @importFrom dplyr filter left_join mutate rename select
#'
eigen_traits <- function(object,
                          sexname = sexnames,
                          modulename,
                          contr_object,
                          eigen_object = eigen_contrast(object, contr_object)) {
  
  eigen_object$p.value <- 1e-10
  
  object <- dplyr::filter(
    object[[sexname]]$modules,
    .data$module == modulename)
  
  object <- 
    dplyr::mutate(
      # Join contrast traits with `1 - |kME|` as p.value.
      dplyr::left_join(
        dplyr::mutate(
          dplyr::select(
            # Filter contrast object to include module traits.
            dplyr::filter(
              contr_object,
              .data$sex == sexname,
              .data$trait %in% object$trait),
            -p.value),
          trait = as.character(trait)),
        dplyr::mutate(
          dplyr::rename(
            dplyr::select(object, -module),
            p.value = "kME"),
          p.value = 1 - abs(.data$p.value))),
      trait = reorder(.data$trait, .data$p.value))
  
  object <- 
    # Reorder traits to be in increasing `p.value` order.
    dplyr::mutate(
      dplyr::bind_rows(
      object,
      # Make sure `p.value` for module is smaller than all others.
      dplyr::mutate(
        # Filter eigen object to `sex` and `module`.
        dplyr::filter(
          eigen_object,
          .data$sex == sexname,
          .data$trait == modulename),
        p.value = min(object$p.value) / 10)),
      trait = reorder(as.character(trait), p.value))
  
  class(object) <- c("conditionContrasts", class(object))
  attr(object, "conditions") <- attr(contr_object, "conditions")
  attr(object, "termname") <- attr(contr_object, "termname")
  attr(object, "datatype") <- "eigen"
  object
}