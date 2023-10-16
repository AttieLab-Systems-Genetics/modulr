#' Parameters of WGCNA object
#'
#' @param params list of updated parameters (defaults if `NULL`)
#'
#' @return list with parameters
#' @export
wgcna_params <- function(params = NULL) {
  defaults <- list(
    signType = "unsigned",
    power = 6, 
    minSize = 4,
    method = "average",
    cutHeight = 0.995,
    split = 2,
    thresholdKME = 0.365,
    thresholdMEDiss = 0.25,
    verbose = 0
  )
  if(is.null(params))
    return(defaults)
  
  for(i in names(defaults)) {
    if(is.null(params[[i]]))
      params[[i]] <- defaults[[i]]
  }
  if(!(params$signType %in% c("unsigned","signed")))
    params$signType <- "unsigned"
  
  params
}
