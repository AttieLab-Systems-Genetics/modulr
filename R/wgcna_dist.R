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
