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
  geneTree <- fastcluster::hclust(
    as.dist(dissTOM), 
    method = params$method)
  
  # Create modules
  mods <- dynamicTreeCut::cutreeDynamic(
    dendro = geneTree,
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

  # determine kMEs
  kME <- WGCNA::signedKME(object, eigen)
  kME <-
    # Left join with trait names and colors to extract kME value.
    module_factors(
      dplyr::left_join(
        dplyr::tibble(
          trait = row.names(kME),
          module = colors),
        # Pivot kME to long form 
        dplyr::mutate(
          tidyr::pivot_longer(
            dplyr::mutate(
              kME,
              trait = row.names(kME)),
            -trait,
            names_to = "module", values_to = "kME"),
          module = stringr::str_remove(module, "^kME")),
        by = c("trait", "module")),
      "module")

  names(eigen) <- stringr::str_remove(names(eigen), "^ME")
  eigen <- eigen[levels(kME$module)]
  
  out <- list(
    ID = ID,
    geneTree = geneTree,
    eigen = eigen,
    modules = kME)
  
  class(out) <- c("wgcnaModules", class(out))
  out
}
#' Create List of WGCNA Modules
#'
#' @param traitData data frame from `foundr::traitData()`
#' @param traitSignal data frame from `foundr::partition()`
#'
#' @return object of class `listof_wgcnaModules`
#' @export
#' @importFrom foundr join_signal
#' @importFrom dplyr rename select
#' 
#' @rdname wgcnaModules
#'
listof_wgcnaModules <- function(traitData, traitSignal) {
  out <- list(
    individual = wgcnaModules(traitData),
    cellmean   = wgcnaModules(
      dplyr::select(
        dplyr::rename(traitSignal, value = "cellmean"),
        -signal)),
    signal     = wgcnaModules(
      dplyr::select(
        dplyr::rename(traitSignal, value = "signal"),
        -cellmean)),
    rest       = wgcnaModules(
        foundr::join_signal(
          traitData,
          traitSignal,
          "rest")),
    noise      = wgcnaModules(
      foundr::join_signal(
        traitData,
        traitSignal,
        "noise")))
  class(out) <- c("listof_wgcnaModules", class(out))
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
    x$geneTree,
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
