#' ID data frame from 
#'
#' @param object data frame containing ID information
#' @param condition_under condition has embedded "_" if `TRUE`
#'
#' @return data frame with column `ID` and optional column `animal`
#' @export
#' @importFrom dplyr arrange distinct mutate
#' @importFrom tidyr all_of separate_wider_delim unite
#' @importFrom stringr str_replace
#' @importFrom rlang .data
#'
wgcna_ID <- function(object, condition_under = TRUE) {
  if(inherits(object, "wgcnaModules"))
    wgnca_ID_object(object)
  else
    wgcna_ID_ME(object, condition_under)
}

# ID data frame from eigentraits rownames
wgcna_ID_ME <- function(MEs, condition_under = TRUE) {
  IDobj <- dplyr::tibble(ID = rownames(MEs))
  if(condition_under) {
    IDobj <- dplyr::mutate(
      IDobj,
      ID = stringr::str_replace(.data$ID, "_([A-Z]+)$", ":\\1"))
  }
  # Change "." after strain to "_"
  IDobj <- dplyr::mutate(
    IDobj,
    ID = stringr::str_replace(.data$ID, "\\.", "_"))
  
  IDcols <- c("strain", "animal", "sex", "condition")
  
  IDobj <- 
    tidyr::unite(
      dplyr::mutate(
        tidyr::separate_wider_delim(
          IDobj,
          ID,
          delim = "_",
          names = IDcols),
        strain = ifelse(.data$strain == "X129", "129", strain)),
      ID,
      strain, sex, condition)
  
  # Change back
  if(condition_under) {
    IDobj <- dplyr::mutate(
      IDobj,
      ID = stringr::str_replace(.data$ID, ":", "_"))
  }
  
  IDobj
}

# ID data frame from module object rownames
wgcna_ID_object <- function(object) {
  # ID columns other than animal
  IDcols <- c("dataset", "strain", "sex", "condition")
  m <- match(IDcols, names(object), nomatch = 0)
  IDcols <- IDcols[m > 0]
  
  # Unite IDcols into `ID`
  IDobj <-
    tidyr::unite(
      object,
      ID, tidyr::all_of(IDcols))
  
  # Add column for `animal` if present.
  if("animal" %in% names(object)) {
    IDobj <- 
      dplyr::arrange(
        dplyr::distinct(
          IDobj,
          ID, animal),
        ID, animal)
  } else {
    IDobj <- 
      dplyr::arrange(
        dplyr::distinct(
          IDobj,
          ID),
        ID)
  }
  
  IDobj
}
  