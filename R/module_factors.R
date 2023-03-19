module_factors <- function(object, module) {
  modcolors <- dplyr::arrange(
    dplyr::count(
      object,
      .data[[module]]),
    dplyr::desc(n))[[module]]
  
  object[[module]] <- factor(object[[module]], modcolors)
  object
}