
find_doublets <- function(object, seed = 123456) {
  
  
  
  dbl <- Seurat::as.SingleCellExperiment(object)
  
  
  
  top.mam <- scran::getTopHVGs(dbl, prop=0.1)
  
  dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                                 
                                                 subset.row=top.mam,
                                                 
                                                 d=ncol(reducedDim(dbl)))
  
  
  
  dbl$DoubletScore <- dbl.dens
  
  
  
  stopifnot(colnames(dbl) == colnames(object))
  
  object$DoubletScore <- dbl$DoubletScore
  
  
  
  rm(dbl, dbl.dens, top.mam)
  
  gc()
  
  object$DoubletScoreLog <- log1p(object$DoubletScore)
  
  return(object)
  
}



####loading datasets####



#####1#####

# Loading and pre-processing dataset 1

object_1 = Read10X(data.dir = "/outs/filtered_feature_bc_matrix/")



doubletcheck = CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

doubletcheck = NormalizeData(object = doubletcheck, normalization.method = "LogNormalize", scale.factor = 1e4)

doubletcheck = FindVariableFeatures(object = doubletcheck, selection.method = "vst", nfeatures = 2000)

doubletcheck = ScaleData(doubletcheck, features = rownames(doubletcheck))

doubletcheck = RunPCA(doubletcheck)

doubletcheck = FindNeighbors(doubletcheck, dims = 1:30)

doubletcheck = FindClusters(doubletcheck, resolution = 1)

doubletcheck = RunUMAP(doubletcheck, dims = 1:30)

doubletcheck = find_doublets(object = doubletcheck)



a = FeaturePlot(doubletcheck, "DoubletScore")

b = ggplot(data.frame("doublets" = doubletcheck$DoubletScore), aes(doublets)) +
  
  geom_density() +
  
  theme_bw()

a+b+plot_layout(ncol = 2)



doublets_1 = names(doubletcheck$DoubletScoreLog[doubletcheck$DoubletScoreLog > 2])



object_1  = object_1[,!colnames(object_1) %in% doublets_1]

object_1 = CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object = grep(pattern = "^MT-", x=rownames(x=object_1), value=T)

percent.mito_object = Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))

object_1[["percent.mito"]] = percent.mito_object

VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5)

object_1 = subset(x = object_1, subset = nFeature_RNA > 1000 & nCount_RNA > 200 & percent.mito <0.1)

object_1 = NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)

object_1 = FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)

object_1[["state"]] = "P1"


