#' This function takes in a seurat object and returns a named list with 2 objects, 1st the metadata + the coordinates of the 2D embeding, the later with names coord_x and coord_y, and second the expression matrix selected.
#'
#' @param se_obj Object of class Seurat from which we want to extract the information.
#' @param assay Object of class Character indicating from which assay to extract expression data.
#' @param slot Object of class Character indicating from which slot to extract expression data, by default data.
#' @param reduction Object of class Character indicating from which dimensionality reduction we want to extract the coordinates.
#' @return This function returns a named list, the first position contains the joint metadata + 2D embeding and the second contains the expression data matrix. 2D dimensions are labelles as coord_x and coord_y.
#' @export
#' @examples
#'

prepare_se_obj <- function(se_obj,
                           assay,
                           slot = "data",
                           reduction = "umap") {
  
  suppressMessages(library(dplyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tibble))
  suppressMessages(library(SummarizedExperiment))
  
  # Input check
  if (! is(se_obj, "Seurat")) stop("Object passed is not a Seurat object;")
  if (! assay %in% Seurat::Assays(se_obj)) stop("assay not in the Seurat object's available assays.")
  if (! slot %in% c("counts", "data", "scale.data")) warning("slot not in the Seurat object's assays.")
  if (! reduction %in% names(se_obj@reductions)) stop("reduction not in the Seurat object's available reductions.")
  
  # Extract 2D coordinates
  embed_df <- se_obj@reductions[[reduction]]@cell.embeddings %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode")
  
  # Change 2D coordinate names to coord_x and coord_y
  colnames(embed_df) <- c("barcode", "coord_x", "coord_y")
  
  # Join metadata with coordinates
  metadata_df <- se_obj@meta.data %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(embed_df, by = "barcode")
  
  # Extract expression data
  assay_data <- Seurat::GetAssayData(se_obj,
                                     slot = slot,
                                     assay = assay)
  
  return(list(metadata = metadata_df,
              expr_data = assay_data))
  }
