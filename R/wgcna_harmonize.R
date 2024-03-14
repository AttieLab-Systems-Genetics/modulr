#' Harmonize WGCNA module
#'
#' This loads `geneTree`, `kMEs`, `merge`, `params` and `dynamicColors_2`
#' from an `Rdata` object and creates and `RDS` object for later reuse.
#' 
#' @param mod name for module object
#' @param moddir director
#' @param modRdata file name ending with `.Rdata`
#' @param params list of parameters
#' @param annot annotation file (ignored if `NULL`)
#' @param harmonizedir name of directory to save `RDS` object in its `moddir`
#' @param force force creation if `TRUE`
#'
#' @return invisible file name for created object
#' @export
wgcna_harmonize <- function(mod = "WGCNAmodule",
                            moddir = ".",
                            modRdata = NULL,
                            params = NULL,
                            annot = NULL,
                            harmonizeddir = ".",
                            force = FALSE) {
  if(is.null(modRdata))
    return(NULL)
  
  # WGCNA Parameters
  if(is.null(params))
    params <- list(
      power = 12, 
      signType = "unsigned",
      minSize = 4)
  paramcode <- params
  if(is.null(paramcode$signType))
    paramcode$signType <- "U"
  else
    paramcode$signType <- ifelse(paramcode$signType == "unsigned", "U", "S")
  
  modname <- paste0(mod, "Module")
  if(file.exists(
    filename <- file.path(harmonizeddir, mod,
                          paste0(modname, "_", paste(paramcode, collapse = ""), ".rds"))) &
    !force) {
    assign(modname, readRDS(filename))
  } else {
    assign(
      modname,
      modulr::load_wgcna(
        moddir,
        modRdata,
        params = params,
        annot = annot))
    if(!dir.exists(file.path(harmonizeddir, mod)))
      dir.create(file.path(harmonizeddir, mod))
    saveRDS(get(modname), filename)
  }
}
