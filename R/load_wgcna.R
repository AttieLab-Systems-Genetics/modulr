#' Load WGCNA object
#'
#' @param moddir directory name containing module object
#' @param modobj name of module object, ending in `.Rdata` or `.RData`
#' @param params non-default parameters for WGCNA (see `wgcna_params`)
#'
#' @return list object that has all WGCNA components
#' @export
#'
load_wgcna <- function(moddir, modobj = "WGCNA_objects_ms10.Rdata",
                       params = list(signType = "unsigned",
                                     power = 12, 
                                     minSize = 4)) {
  # local() keeps loaded objects local.
  out <- local({
    load(file.path(moddir, modobj))
    # operations to create object, which is returned.
    
    list(
      ID = wgcna_ID(merge$newMEs),
      dendro = merge$dendro,
      eigen = merge$newMEs,
      modules = module_factors(kMEs, merge$colors),
      params = wgcna_params(params))
  })
  
  class(out) <- c("wgcnaModules", class(out))
  attr(out, "params") <- out$params
  out
}
