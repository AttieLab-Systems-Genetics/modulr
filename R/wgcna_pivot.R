#' Pivot WGCNA object
#'
#' @param object of class `wgcnaModules`
#'
#' @return pivoted object in wide format for WGCNA routines
#' @export
#' @importFrom dplyr arrange distinct
#' @importFrom tidyr all_of unite
wgcna_pivot <- function(object) {
  # Pivot object
  IDcols <- c("dataset", "strain", "sex", "condition", "animal")
  m <- match(IDcols, names(object), nomatch = 0)
  IDcols <- IDcols[m > 0]
  
  IDobj <- wgcna_ID(object)

  object <- 
    as.data.frame(
      tidyr::pivot_wider(
        tidyr::unite(
          object,
          ID, tidyr::all_of(IDcols)),
        names_from = "trait", values_from = "value"))
  
  rownames(object) <- object$ID
  list(
    matrix = as.matrix(object[,-1]),
    ID = IDobj)
}
