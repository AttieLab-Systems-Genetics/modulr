#' Create WGCNA Modules
#'
#' Data in `object` are assumed to be in long format with columns
#' a subset of "dataset", "strain", "sex", "animal", "condition"
#' 
#' @param object harmonized data frame from routine `foundr`
#' @param params list of parameters for WGCNA routines
#'
#' @return object of class wgcnaModules
#' 
#' @export
#' @importFrom WGCNA adjacency labels2colors TOMsimilarity
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stringr str_remove
#'
wgcnaModules <- function(object, params = NULL) {

  # Pivot object to have traits in columns and ID in rownames.
  object <- wgcna_pivot(object)
  
  ID <- object$ID
  object <- object$matrix
  
  # Check and assign parameters
  params <- wgcna_params(params)
  
  # Topological overlap (TOM)
  dissTOM <- wgcna_dist(object, params)

  # hierarchical cluster
  dendro <- fastcluster::hclust(
    as.dist(dissTOM), 
    method = params$method)
  
  # Create modules
  mods <- dynamicTreeCut::cutreeDynamic(
    dendro = dendro,
    cutHeight = params$cutHeight,
    distM = dissTOM,
    deepSplit = params$split,
    pamRespectsDendro = TRUE,
    minClusterSize = params$minSize,
    verbose = params$verbose)
  
  # Colors for modules
  colors <- WGCNA::labels2colors(mods)
  
  # Merge close modules based on 
  merge <- WGCNA::mergeCloseModules(
    object, 
    colors,
    unassdColor="grey",
    useAbs = FALSE,
    relabel = TRUE,
    getNewMEs = TRUE,
    getNewUnassdME = TRUE,
    iterate = TRUE,
    cutHeight = params$thresholdMEDiss,
    verbose = params$verbose)
  
  colors <- merge$colors
  eigen <- merge$newMEs

  # determine `kME` and add color by `module`.
  kME <- module_factors(
    WGCNA::signedKME(object, eigen),
    colors)

  names(eigen) <- stringr::str_remove(names(eigen), "^ME")
  eigen <- eigen[levels(kME$module)]
  
  out <- list(
    ID = ID,
    dendro = dendro,
    eigen = eigen,
    modules = kME)
  
  class(out) <- c("wgcnaModules", class(out))
  attr(out, "params") <- params
  out
}



#' Create List of WGCNA Modules
#'
#' @param traitData data frame from `foundr::traitData()`
#' @param traitSignal data frame from `foundr::partition()`
#' @param params list of parameters for WGCNA routines
#'
#' @return object of class `listof_wgcnaModules`
#' @export
#' @importFrom foundr join_signal
#' @importFrom dplyr rename select
#' 
#' @rdname wgcnaModules
#'
listof_wgcnaModules <- function(traitData, traitSignal, params = NULL) {
  out <- list(
    value      = wgcnaModules(traitData, params),
    cellmean   = wgcnaModules(
      dplyr::select(
        dplyr::rename(traitSignal, value = "cellmean"),
        -signal),
      params),
    signal     = wgcnaModules(
      dplyr::select(
        dplyr::rename(traitSignal, value = "signal"),
        -cellmean),
      params),
    rest       = wgcnaModules(
        foundr::join_signal(
          traitData,
          traitSignal,
          "rest"),
        params),
    noise      = wgcnaModules(
      foundr::join_signal(
        traitData,
        traitSignal,
        "noise"),
      params))
  class(out) <- c("listof_wgcnaModules", class(out))
  attr(out, "params") <- attr(out[[1]], "params")
  out
}
#' @export
#' @rdname wgcnaModules
#' @importFrom WGCNA plotDendroAndColors
#' 
plot_wgcnaModules <- function(x,
                              main = "Gene dendrogram and module colors",
                              ...) {
  WGCNA::plotDendroAndColors(
    x$dendro,
    x$modules$module,
    "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = FALSE, guideHang = 0.05,
    main = main)
}
#' @export
#' @rdname wgcnaModules
#' @method plot wgcnaModules
plot.wgcnaModules <- function(x, ...)
  plot_wgcnaModules(x, ...)
