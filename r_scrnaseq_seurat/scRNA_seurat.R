library(tidyverse)
library(ggplot2)
library(patchwork)
library(Seurat)
library(dplyr)
library(clustree)

# input data

pbmc.data <- Read10X(data.dir = "scRNA/filtered_feature_bc_matrix/")
class(pbmc.data)
View(pbmc.data)
# initial filter, genes have to express in more than 3 cells, and 
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`,
                               project = "pbmc3k")

pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

Assays(pbmc)
#RNA is the default already, don't need to change it yet

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
View(pbmc[[]])

#do some plot for mt genes
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# set some filters
VlnPlot(pbmc, features = c("nFeature_RNA"))+
        geom_line(yintercept = 1250)
        
VlnPlot(pbmc, features = c("nCount_RNA"))+
        geom_line(yintercept = 1250)
                
VlnPlot(pbmc, features = c("percent.mt"))+
        geom_line(yintercept = 14)

# look at the ribsome gene expression
pbmc[["percent.rib"]] <- PercentageFeatureSet(pbmc,
                                              pattern =  "^RPS|^RPL")
VlnPlot(pbmc, features = c("percent.rib"))
# pbmc has less then 50% of ribsome gene

# QC plot of nCount_RNA by nFeature_RNA colored by percent.mt(use ggplot)

ggplot(pbmc[[]], aes(x = nCount_RNA, y=nFeature_RNA, color = percent.mt))+
    geom_point()
# # Suggestion: if just looking at this plot, can cut off at nFeature_RNA = 6000 # But this is an iterative process 

# Cutoffs:

# # 1000 for nFeature_RNA

# # 12.5 for percent.mt

# # Leave nCount_RNA for now

# # Leave percent.rib for now (but some people say that if it's above 50% in PBMCs then they should be filtered out, but be careful)



plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# normalize
pbmc <- SCTransform(pbmc,
                    assay = 'RNA',
                    seed.use = 1448145,
                    verbose = TRUE)

# This procedures a new assay which has all the sc transform stuff, lookc as ducumentation of SCTransform for defaults, particulary.
# New.assay.name  ="SCT", variable.features.n =3000, vars.to.regress = NULL(can regress out precent.mt if wanted)
# but do not try to regress out BIG batch effects as this is only meant for samll-y effecting variables
## Remember that SCTransform performs these three functions:

# NormalizeData

# FindVariableFeatures

# ScaleData

# warning: dont worry about this one,it seems to work fine
DefaultAssay(pbmc)
# SCTranform automatically changes default assay to SCT

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(pbmc, reduction = "pca", group.by = "percent.mt")+
    NoLegend()

print(pbmc[["pca"]], dims = 1:5, nFeatures = 5)
# print the PC loadings: genes that most contribute to PC

VisDimloadings(pbmc, dims = 1:3, reduction = "pca")
# plots above print statement essentially

DimHeatmap(pbmc,dims = 1, cells = 2000, balanced = TRUE )

# Dimheatmap can give an indication of which PCs we should included in the dim red
# JackStraw cannot bu run on SCTramsform normalized data
# pbmc <- Jackstraw(pbmc, num.replicate =100)


ElbowPlot(pbmc, ndims = 30)

# using SCTransform can include a few more PCs than otherwise, so error side of increasing
# in this section, we will choose between 20-30

# Clustering the cellls

pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)
#k.param play around with this , depending on size of clusters
# annoy = appromimate nearest neighbour is another P package

pbmc <- FindClusters(pbmc, resolution = c(0.5, 0.8,1,1.3,2))
# From vignette, we find that setting this parameter [resolution] between 0.4-1.2 returns good results for single cell
# dataset of around about 3k cells
# affected by the numbers of cells you used
# 
View(pbmc[[]])

Idents(pbmc) <- "SCT_snn_res.0.8"
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set "lable = TRUE " or use the lableClusters fucntion to help label
DimPlot(pbmc, reduction = "umap")

#highlight the cells of your interest
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nCount_RNA > 11000))

DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = percent.mt>12.5))

DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nFeature_RNA<1000))

#the thresholds look quite reasonable at this stage

ggplot(pbmc[[]], aes(x = nFeature_RNA))+
    geom_histogram(bins = 100)

pbmc <- subset(pbmc, 
               subset = nFeature_RNA >1000 &
                   nFeature_RNA <5500 &
                   percent.mt < 12.5)

#==========
# re-run the normalization
pbmc <- SCTransform(pbmc,
                    assay = 'RNA',
                    seed.use = 1448145,
                    verbose = TRUE)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)


pbmc <- FindClusters(pbmc, resolution = c(0.5, 0.8,1,1.3,2))

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

FeaturePlot(pbmc, features  = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A","LYZ", "PPBP", "CD8A"))

# post-filtering steps
#find out how many clusters we should use

clustree(pbmc)
# let's choose a resolution of 0.8
Idents(pbmc) <- "SCT_snn_res.0.8"



#run the protein 
pbmc
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay = "ADT")

DefaultAssay(pbmc) <- "ADT"
rownames(pbmc)
FeaturePlot(pbmc, features = c("CD3-TotalSeqB", "CD11b-TotalSeqB", "CD8a-TotalSeqB", "CD16-TotalSeqB", "PD-1-TotalSeqB", "TIGIT-TotalSeqB", "CD14-TotalSeqB", "CD4-TotalSeqB"))

RidgePlot(pbmc, features = c("CD3-TotalSeqB", "CD11c-TotalSeqB", "CD8a-TotalSeqB", "CD16-TotalSeqB"), ncol = 2)


# do pca on protein
IgG <- c("IgG1-control-TotalSeqB", "IgG2a-control-TotalSeqB", "IgG2b-control-TotalSeqB")

pbmc <- RunPCA(pbmc, 
               features = 
                   rownames(pbmc)[!rownames(pbmc) %in% IgG],
               npcs = 10,
               reduction.name = "pca_adt", 
               reduction.key = "pca_adt_", 
               verbose = FALSE)

Reductions(pbmc) # numbers of reduction that you have made

DimPlot(pbmc, reduction = "pca_adt")

pbmc <- FindNeighbors(pbmc, features = rownames(pbmc)[!rownames(pbmc) %in% IgG], dims = NULL)

pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.2, 0.3, 0.4), 
                     graph.name = "ADT_snn")

DefaultAssay(pbmc) <- "ADT"
Idents(pbmc) <- "ADT_snn_res_0.1"

clustering.table <- table(Idents(pbmc), pbmc@meta.data$SCT_snn_res.0.8)

umap_rnaClusters <- DimPlot(pbmc, reduction = "umap_adt", group.by = "SCT_snn_res.0.8") + NoLegend()
umap_rnaClusters




