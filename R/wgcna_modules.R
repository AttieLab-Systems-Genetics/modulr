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
          module = str_remove(module, "^kME")),
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
#' @importFrom foundr calc_ind_signal
#' 
#' @rdname wgcnaModules
#'
listof_wgcnaModules <- function(traitData, traitSignal) {
  out <- list(
    individual = wgcnaModules(traitData),
    cellmean   = wgcnaModules(traitSignal %>%
                                rename(value = "cellmean") %>%
                                select(-signal)),
    signal     = wgcnaModules(traitSignal %>%
                                rename(value = "signal") %>%
                                select(-cellmean)),
    ind_signal = wgcnaModules(foundr::calc_ind_signal(traitData, traitSignal) %>%
                                select(-signal)))
  class(out) <- c("listof_wgcnaModules", class(out))
  out  
}
#' Summary of List of WGCNA Modules
#'
#' @param object object of class `listof_wgcnaModules`
#' @param ... additional parameters
#'
#' @return data frame
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom purrr map set_names
#' 
#' @rdname wgcnaModules
#' @method summary listof_wgcnaModules
#'
summary.listof_wgcnaModules <- function(object, ...) {
  dplyr::bind_rows(
    purrr::set_names(
      purrr::map(
        names(object),
        function(x) summary(object[[x]])),
      names(object)),
    .id = "response")
}

wgcna_dist <- function(object, params) {
  
  # Pivot object to have traits in columns and ID in rownames.
  if(!is.matrix(object))
    object <- wgcna_pivot(object)$matrix
  
  # Check and assign parameters
  params <- wgcna_params(params)
  
  # Topological overlap (TOM)
  1 - WGCNA::TOMsimilarity(
    WGCNA::adjacency(
      object, 
      type = params$signType,
      power = params$power), 
    TOMType = params$signType,
    verbose = params$verbose)
}

wgcna_pivot <- function(object) {
  # Pivot object
  IDcols <- c("dataset", "strain", "sex", "condition")
  m <- match(IDcols, names(object), nomatch = 0)
  IDcols <- IDcols[m > 0]
  
  if("animal" %in% names(object)) {
    IDobj <- 
      dplyr::arrange(
        dplyr::distinct(
          tidyr::unite(
            object,
            ID, tidyr::all_of(IDcols)),
          ID, animal),
        ID, animal)
    IDcols <- c(IDcols, "animal")
  } else {
    IDobj <- 
      dplyr::arrange(
        dplyr::distinct(
          tidyr::unite(
            object,
            ID, tidyr::all_of(IDcols)),
          ID),
        ID)
  }
  
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
  if(!(params$signType %in% c("unsigned","signed")))
    params$signType <- "unsigned"

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
          dplyr::mutate(
            object$modules,
            module = as.character(module)),
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
