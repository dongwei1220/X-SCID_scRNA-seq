###################################################################################################
# Integrating X-SCID patient and healthy donor scRNA-seq data                                     #
###################################################################################################
rm(list=ls())

# load required packages
library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)

# https://satijalab.org/seurat/v3.0/future_vignette.html
library(future)
# it is tricky to set options
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 30 * 1024^3)

packageVersion("Seurat")
#[1] ‘3.2.2’
packageVersion("future")
#[1] ‘1.19.1’

# set working directory
setwd("/data/dongwei/projects/X-SCID_scRNA-seq/")

Patient_list <- readRDS("Patient_BMMC.rds")

MantonBM_list <- readRDS("MantonBMMC_healthy_donor.rds")

Healthy_list <- readRDS("Healthy_PBMC.rds")

data_list <- c(Patient_list, MantonBM_list, Healthy_list)
data_list
###################################################################################################
# integration with SCTransform
data_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000, verbose = FALSE)

data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = data_features, verbose = FALSE)

# integrate data with selected reference
data_anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", 
                                       anchor.features = data_features, reference = c(1, 3:5, 8), verbose = FALSE)

data_integrated <- IntegrateData(anchorset = data_anchors, normalization.method = "SCT", verbose = FALSE)

# pca dimensional reduction
data_integrated <- RunPCA(object = data_integrated, verbose = FALSE)

#pdf(file='Integrated_DimHeatmap.pdf',height=20,width=12)
#DimHeatmap(data_integrated, dims = 1:15, cells = 500, balanced = TRUE)
#dev.off()
pdf(file='Integrated_ElbowPlot.pdf',height=8,width=11)
ElbowPlot(data_integrated, ndims=50)
dev.off()
pdf(file='Integrated_PCA_DimPlot.pdf',height=8,width=10)
DimPlot(data_integrated, reduction = "pca", group.by="sample")
dev.off()

# Cluster the cells
# The clusters are saved in the object@ident slot
data_integrated <- FindNeighbors(data_integrated, reduction = "pca", dims = 1:50)

# with increased resolution values leading to a greater number of clusters
data_integrated <- FindClusters(data_integrated, resolution = 0.4)
data_integrated <- FindClusters(data_integrated, resolution = 0.5)
data_integrated <- FindClusters(data_integrated, resolution = 0.6)
data_integrated <- FindClusters(data_integrated, resolution = 0.7)
data_integrated <- FindClusters(data_integrated, resolution = 0.8)

# Run non-linear dimensional reduction
#The goal of these algorithms is to learn the underlying manifold of the data in order to place 
# similar cells together in low-dimensional space.
data_integrated <- RunUMAP(object = data_integrated, dims = 1:50, verbose = FALSE)

# chose resolution=0.6 -> 30 clusters
Idents(object = data_integrated) <- 'integrated_snn_res.0.6'
levels(data_integrated) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29")

###################################################################################################
# data visualization
DefaultAssay(data_integrated) <- "integrated"

# group plot
p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "sample") & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))
p2 <- DimPlot(data_integrated, reduction = "umap", label = T) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))
pdf(file='Integrated_DimPlot.umap_group.all.pdf',height=12,width=26)
p1 + p2
dev.off()

# group plot by cluster
cluster_colors <- c("#5A2955", "#D51F26", "#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#DCABCF", "#FDDD03", "#F69421", "#9C82BA", "#AA875A", "#3EB3A7", "#245359", "#711E21", "#B36B45", "#DB7D8D", "#BF5D58", "#216069", "#9A7456", "#8ED2D0", "#3D3D3D", "#B4B883", "#C39C6F", "#9B8EC4", "#6264A0", "#BFCADF", "#8AC972", "#D51F26", "#272E6A", "#208A42")
pdf(file='Integrated_DimPlot.umap_group.cluster.pdf',height=12,width=13)
p <- DimPlot(data_integrated, reduction = "umap", label=TRUE, label.size = 8) + scale_color_manual(values=cluster_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

# group plot by sample
sample_colors <- rev(c("#a6cee3","#1f78b4","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#b15928"))
pdf(file='Integrated_DimPlot.umap_group.sample.pdf',height=12,width=14)
DimPlot(data_integrated, reduction = "umap", group.by = "sample", label=FALSE) + xlim(-14,12) + ylim(-11,18) + scale_color_manual(values=sample_colors)
dev.off()

# split plot by sample
pdf(file='Integrated_DimPlot.umap_split.sample.pdf',height=12,width=98)
DimPlot(data_integrated, reduction = "umap", split.by = "sample", label=FALSE) + xlim(-14,12) + ylim(-11,18) + NoLegend()
dev.off()

# split plot by cluster
# integrated_snn_res.0.6 is set for Idents
pdf(file='Integrated_DimPlot.umap_split.cluster.pdf',height=6,width=98)
DimPlot(data_integrated, reduction = "umap", split.by = "integrated_snn_res.0.6", label=FALSE) + xlim(-14,12) + ylim(-11,18) + NoLegend() + scale_color_manual(values=cluster_colors)
dev.off()

########################################################################################################################
DefaultAssay(data_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
data_integrated <- NormalizeData(data_integrated, verbose = FALSE)

saveRDS(data_integrated, file="All_patient_and_healthy_donors_BMMC_PBMC.rds")
