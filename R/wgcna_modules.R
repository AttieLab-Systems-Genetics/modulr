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
#'
wgcnaModules <- function(object, params = NULL) {

  # Pivot object to have traits in columns and ID in rownames.
  object <- wgcna_pivot(object)
  
  # Check and assign parameters
  params <- wgcna_params(params)
  
  # Topological overlap (TOM)
  dissTOM <- 
    1 - WGCNA::TOMsimilarity(
      WGCNA::adjacency(
        object, 
        type = params$signType,
        power = params$power), 
      TOMType = params$signType,
      verbose = params$verbose)
  
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
    left_join(
      tibble(trait = row.names(kME),
             module = colors),
      # Pivot kME to long form 
      kME %>%
        mutate(trait = row.names(kME)) %>%
        pivot_longer(-trait, names_to = "module", values_to = "kME") %>%
        mutate(module = str_remove(module, "^kME")),
      by = c("trait", "module"))

  names(eigen) <- stringr::str_remove(names(eigen), "^ME")
  
  out <- list(
    dissTOM = dissTOM,
    geneTree = geneTree,
    eigen = eigen,
    modules = kME)
  
  class(out) <- c("wgcnaModules", class(out))
  out
}

wgcna_pivot <- function(object) {
  # Pivot object
  IDcols <- c("dataset", "strain", "sex", "animal", "condition")
  IDcols <- IDcols[IDcols %in% names(object)]
  
  object <- 
    as.data.frame(
      tidyr::pivot_wider(
        tidyr::unite(
          object,
          ID, tidyr::all_of(IDcols)),
        names_from = "trait", values_from = "value"))
  
  rownames(object) <- object$animalID
  as.matrix(object[,-1])
}

wgcna_params <- function(params = NULL) {
  defaults <- list(
    signType = "unsigned",
    power = 6, 
    minSize = 20,
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
  if(!(params$type %in% c("unsigned","signed")))
    params$type <- "unsigned"

  params
}

#' @param object,x object of class `wgcnaModules`
#' @param main title for plot
#' @param ... additional parameters

#' @export
#' @rdname wgcnaModules
summary_wgcnaModules <- function(object, ...) {
  dplyr::arrange(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          object$modules,
          module),
        count = dplyr::n(),
        maxkME = signif(max(kME), 4),
        minkME = signif(min(kME), 4))),
    dplyr::desc(count))
}
#' @export
#' @rdname wgcnaModules
#' @method summary wgcnaModules
summary.wgcnaModules <- function(object, ...)
  summary_wgcnaModules(object, ...)

#' @export
#' @rdname wgcnaModules
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
