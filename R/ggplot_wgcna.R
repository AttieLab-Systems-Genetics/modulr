#' GGplot of WGCNA Modules
#'
#' @param object object of class `wgcnaModules`
#' @param main title of plot
#' @param ... additional parameters
#'
#' @return ggplot2 object
#' @importFrom ggdendro ggdendrogram
#' @importFrom ggplot2 autoplot ylim
#' 
#' @rdname wgcnaModules
#'
ggplot_wgcnaModules <- function(object,
                              main = "Gene dendrogram and module colors",
                              ...) {
  ggplot_DendroAndColors(
    object$geneTree,
    object$modules$module,
    "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = FALSE, guideHang = 0.05,
    main = main)
}
#' Autoplot of wgcnaModules
#'
#' @param object 
#' @param ... 
#'
#' @return ggplot2 object
#' @rdname wgcnaModules
#' @method autoplot wgcnaModules
#'
autoplot.wgcnaModules <- function(object, ...) {
  gplot_wgcnaModules(object, ...)
}

ggplot_DendroAndColors <- function (dendro, colors, groupLabels = NULL, rowText = NULL, 
    rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
    textPositions = NULL, setLayout = TRUE, autoColorHeight = TRUE, 
    colorHeight = 0.2, colorHeightBase = 0.2, colorHeightMax = 0.6, 
    rowWidths = NULL, dendroLabels = NULL, addGuide = FALSE, 
    guideAll = FALSE, guideCount = 50, guideHang = 0.2, addTextGuide = FALSE, 
    cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8, 
    marAll = c(1, 5, 3, 1), saveMar = TRUE, abHeight = NULL, 
    abCol = "red", ...) 
{
    if (!is.null(dim(colors))) {
        nRows = dim(colors)[2]
    }
    else nRows = as.numeric(length(colors) > 0)
    
    if (!is.null(rowText)) 
        nRows = nRows + if (is.null(textPositions)) 
            nRows
        else length(textPositions)
    
    if (autoColorHeight) 
        colorHeight = colorHeightBase + (colorHeightMax - colorHeightBase) * 
            (1 - exp(-(nRows - 1)/6))
    
    if (setLayout) 
        layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, 
            colorHeight))
    
    p1 <- ggdendro::ggdendrogram(dendro, labels = dendroLabels, size = cex.dendroLabels, 
        ...) +
      ggplot2::ylim(c(0.1,1)) # This causes message; how to pass as arg?

    if (!is.null(abHeight)) 
      p1 <- p1 + ggplot2::geom_hline(yintercept = abHeight, col = abCol)

    plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels, 
        rowText = rowText, rowTextAlignment = rowTextAlignment, 
        rowTextIgnore = rowTextIgnore, textPositions = textPositions, 
        cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
    if (saveMar) 
        par(mar = oldMar)
}
