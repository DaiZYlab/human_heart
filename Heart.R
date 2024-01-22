library(Seurat)
library(patchwork)

# fetal heart
GSM4008686_Fetal_Heart1_dge.data=read.table("~/R/Heart/GSM4008686_Fetal-Heart1_dge.txt",sep="\t",header=TRUE,row.names=1)
GSM4008687_Fetal_Heart2_dge.data=read.table("~/R/Heart/GSM4008687_Fetal-Heart2_dge.txt",sep="\t",header=TRUE,row.names=1)
Fetal_Heart1 <- CreateSeuratObject(counts = GSM4008686_Fetal_Heart1_dge.data, project = "Fetal_Heart1", min.cells = 5, min.features = 100)
Fetal_Heart2 <- CreateSeuratObject(counts = GSM4008687_Fetal_Heart2_dge.data, project = "Fetal_Heart2", min.cells = 5, min.features = 100)

Fetal_Heart1[["percent.mt"]] <- PercentageFeatureSet(object = Fetal_Heart1, pattern = "^MT-")
VlnPlot(Fetal_Heart1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Fetal_Heart2[["percent.mt"]] <- PercentageFeatureSet(object = Fetal_Heart2, pattern = "^MT-")
VlnPlot(Fetal_Heart2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Fetal_Heart1  <- subset(Fetal_Heart1, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 10)
Fetal_Heart2  <- subset(Fetal_Heart2, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 10)
Fetal_Heart1 <- NormalizeData(object =Fetal_Heart1, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal_Heart2 <- NormalizeData(object =Fetal_Heart2, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal_Heart1 <- FindVariableFeatures(object = Fetal_Heart1, selection.method = "vst", nfeatures = 2000)
Fetal_Heart2 <- FindVariableFeatures(object = Fetal_Heart2, selection.method = "vst", nfeatures = 2000)
saveRDS(Fetal_Heart1, file = "~/R/Heart/Fetal_Heart1.rds")
saveRDS(Fetal_Heart2, file = "~/R/Heart/Fetal_Heart2.rds")

Fetal_Heart.anchors <- FindIntegrationAnchors(object.list = list(Fetal_Heart1,Fetal_Heart2), dims = 1:30)
Fetal_Heart <- IntegrateData(anchorset = Fetal_Heart.anchors, dims = 1:20)
DefaultAssay(object = Fetal_Heart) <- "integrated"
Fetal_Heart <- ScaleData(object = Fetal_Heart, verbose = FALSE)
Fetal_Heart<- RunPCA(object = Fetal_Heart, npcs = 30, verbose = FALSE)
ElbowPlot(Fetal_Heart)
Fetal_Heart <- RunUMAP(Fetal_Heart, reduction = "pca", dims = 1:14)
Fetal_Heart <- RunTSNE(Fetal_Heart, dims = 1:14, method = "FIt-SNE")
Fetal_Heart <- FindNeighbors(Fetal_Heart, reduction = "pca", dims = 1:14)
Fetal_Heart <- FindClusters(Fetal_Heart, resolution = 0.5)
DimPlot(Fetal_Heart, reduction = "umap", label = TRUE)
DimPlot(Fetal_Heart, reduction = "tsne", label = TRUE)
DimPlot(Fetal_Heart, reduction = "umap", split.by = "orig.ident", label = TRUE)
DimPlot(Fetal_Heart, reduction = "tsne", split.by = "orig.ident", label = TRUE)
DefaultAssay(Fetal_Heart) <- "RNA"
Fetal_Heart.markers <- FindAllMarkers(Fetal_Heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Fetal_Heart.markers, file = "~/R/Heart/Fetal_Heart.csv")
saveRDS(Fetal_Heart, file = "~/R/Heart/Fetal_Heart.rds")
markers.to.plot <- c("CDH5", "ACTA2","CD68", "NPPA", "NPPB", "PROX1","TNNT2", "MYH7", "COLA1", "VIM", "SOX4", "PDGFRB", "CD3G", "CD19", "MYH11", "FN1", "TOP2A", "KLRC1", "IGHM", "CD14", "CD3E", "NCAM1", "CSPG4", "RGS5", "DES", "FABP4", "TMEM100", "CD74")
DotPlot(Fetal_Heart, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

Fetal_Heart <- RenameIdents(Fetal_Heart, '0' = "Fib", '1' = "Cardiomyocyte", '2' = "Cardiomyocyte",  '3' = "EC", '4' = "Fib", '5' = "EC",  '6' = "SMC", '7' = "Cardiomyocyte", '8' = "Fib",  '9' =  "Mac", '10' = "Cardiomyocyte", '11' = "EC", '12' = "Cardiomyocyte")
Fetal_Heart <- RenameIdents(Fetal_Heart, 'Cardiomyocyte' = "CM")

VlnPlot(Fetal_Heart, features = c("EGLN1", "EGLN2", "EGLN3", "HIF1A", "EPAS1", "VHL"), ncol = 2, pt.size = 0.1)
markers.to.plot <- c("EGLN1", "EGLN2", "EGLN3", "HIF1A", "EPAS1", "VHL")

# Adult
GSM4850581_Heart_Counts.data=read.table("~/R/Heart/GSM4850581_Heart_Counts.csv",sep=",",header=TRUE, row.names=1)
Adult_Heart <- CreateSeuratObject(counts = GSM4850581_Heart_Counts.data, project = "Adult", min.cells = 5, min.features = 100)

Adult_Heart[["percent.mt"]] <- PercentageFeatureSet(object = Adult_Heart, pattern = "^MT-")
VlnPlot(Adult_Heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Adult_Heart  <- subset(Adult_Heart, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 10)
Adult_Heart <- NormalizeData(object =Adult_Heart, normalization.method = "LogNormalize", scale.factor = 10000)
Adult_Heart <- NormalizeData(object =Adult_Heart, normalization.method = "LogNormalize", scale.factor = 10000)
Adult_Heart <- FindVariableFeatures(object = Adult_Heart, selection.method = "vst", nfeatures = 2000)
saveRDS(Adult_Heart, file = "~/R/Heart/Adult_Heart.rds")

Adult_Heart <- ScaleData(object = Adult_Heart, verbose = FALSE)
Adult_Heart<- RunPCA(object = Adult_Heart, npcs = 30, verbose = FALSE)
ElbowPlot(Adult_Heart)
Adult_Heart <- RunUMAP(Adult_Heart, reduction = "pca", dims = 1:7)
Adult_Heart <- RunTSNE(Adult_Heart, dims = 1:7, method = "FIt-SNE")
Adult_Heart <- FindNeighbors(Adult_Heart, reduction = "pca", dims = 1:7)
Adult_Heart <- FindClusters(Adult_Heart, resolution = 0.5)
DimPlot(Adult_Heart, reduction = "umap", label = TRUE)
DimPlot(Adult_Heart, reduction = "tsne", label = TRUE)

Adult_Heart.markers <- FindAllMarkers(Adult_Heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Adult_Heart.markers, file = "~/R/Heart/Adult_Heart.csv")
saveRDS(Adult_Heart, file = "~/R/Heart/Fetal_Heart.rds")
markers.to.plot <- c("CDH5", "ACTA2","CD68", "NPPA", "NPPB", "PROX1","TNNT2", "MYH7", "COLA1", "VIM", "SOX4", "PDGFRB", "CD3G", "CD19", "MYH11", "FN1", "TOP2A", "KLRC1", "IGHM", "CD14", "CD3E", "NCAM1", "CSPG4", "RGS5", "DES", "FABP4", "TMEM100", "CD74")
DotPlot(Adult_Heart, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
RidgePlot(Adult_Heart, features = c("EGLN1", "EGLN2", "EGLN3", "HIF1A", "EPAS1", "VHL"), ncol = 2)


Adult_Heart <- RenameIdents(Adult_Heart, '0' = "Fib", '1' = "Fib", '2' = "Fib",  '3' = "EC", '4' = "EC", '5' = "Mac",  '6' = "SMC", '7' = "Fib", '8' = "T cell",  '9' =  "Fib")
VlnPlot(Adult_Heart, features = c("EGLN1", "EGLN2", "EGLN3", "HIF1A", "EPAS1", "VHL"), ncol = 2, pt.size = 0.1)
markers.to.plot <- c("EGLN1", "EGLN2", "EGLN3", "HIF1A", "EPAS1", "VHL")

FeaturePlot(Adult_Heart, features = c("EPAS1"))
