#!/usr/bin/env R



annotate_read_fractions <- function(seurat_object) {
  loop <- c("percentage_mitochondrial_reads" = "^MT-",
    "percentage_ribosomal_reads" = "^RPL")
  
  for(slot in names(loop)) {
    regex <- loop[[slot]]
    print(slot)
    print(regex)
    print("")
    
    
    target_features <- seurat_object |>
      row.names() |>
      data.frame(hugo_symbol = _) |>
      dplyr::mutate(status = grepl(regex, hugo_symbol)) |>
      (function(.) {
        assertthat::assert_that(sum(.$status) > 1) # must contain at least one gene
        return(.)
      })()
    
    reads_total <- Matrix::colSums(Seurat::GetAssayData(object = seurat_object, slot="counts"))
    reads_target <- Matrix::colSums(Seurat::GetAssayData(object = seurat_object, slot="counts")[target_features$status,])
    
    seurat_object[[slot]] <- reads_target / reads_total * 100.0
  }
  
  
  # colSums crashes with only 1 column, so do MALAT1 separately
  reads_total <- Matrix::colSums(Seurat::GetAssayData(object = seurat_object, slot="counts"))
  reads_target <- Seurat::GetAssayData(object = seurat_object, slot="counts")[c("MALAT1"),]
  seurat_object[["percentage_MALAT1_reads"]] <- reads_target / reads_total * 100.0
  
  
  return(seurat_object)
  
}




