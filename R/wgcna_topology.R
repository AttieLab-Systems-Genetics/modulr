#' WGCNA topology
#'
#' @param object data frame
#' @param powers vector of beta values
#' @param verbose level of verbose messages
#' @param ... ignored
#'
#' @return data frame of class `wgcna_topology`
#' @export
#' @importFrom WGCNA pickSoftThreshold
#'
wgcna_topology <- function(object,
                           powers = c(c(1:10), seq(from = 12, to=20, by=2)),
                           verbose = 0, ...) {
  
  # Pivot object to have traits in columns and ID in rownames.
  object <- wgcna_pivot(object, ...)
  
  ID <- object$ID
  object <- object$matrix

  # Call the network topology analysis function
  out <- WGCNA::pickSoftThreshold(object, powerVector = powers,
                                  verbose = verbose)$fitIndices
  class(out) <- c("wgcna_topology", class(out))
  out
}
#' GGplot of WGCNA Topology
#'
#' @param object object of class  
#' @param cutoff explained variation cutoff
#' @param ... 
#'
#' @return ggplot object
#' @export
#' @rdname wgcna_topology
#' @importFrom ggplot2 aes geom_hline geom_text ggplot ggtitle xlab ylab
#' @importFrom cowplot plot_grid
#'
ggplot_wgcna_topology <- function(object, cutoff = 0.9, ...) {
  # Scale-free topology fit index as a function of the soft-thresholding power
  p1 <- ggplot2::ggplot(object) +
    ggplot2::aes(.data$Power, -sign(.data$slope) * .data$SFT.R.sq,
                 label = .data$Power) +
    ggplot2::xlab("Soft Threshold (power)") +
    ggplot2::ylab("Scale Free Topology Model Fit,signed R^2") +
    ggplot2::ggtitle(paste("Scale independence")) +
    ggplot2::geom_text(col = "red") +
    # this line corresponds to R^2 cutoff.
    ggplot2::geom_hline(yintercept = cutoff, col = "red")
  
  # Mean connectivity as a function of the soft-thresholding power
  p2 <- ggplot2::ggplot(object) +
    ggplot2::aes(.data$Power, .data$mean.k., label = .data$Power) +
    ggplot2::xlab("Soft Threshold (power)") +
    ggplot2::ylab("Mean Connectivity") +
    ggplot2::ggtitle(paste("Mean connectivity")) +
    ggplot2::geom_text(col = "red")
  cowplot::plot_grid(plotlist = list(p1, p2), nrow = 1)
}