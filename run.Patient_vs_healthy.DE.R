###################################################################################################
# DE analysis for patient(child and mother)-vs-healthy sample                                     #
###################################################################################################

# load required packages
library(Seurat)
library(patchwork)
library(ggplot2)

# https://satijalab.org/seurat/v3.0/future_vignette.html
library(future)
# it is tricky to set options
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 30 * 1024^3)

# set working directory
setwd("/data/dongwei/projects/X-SCID_scRNA-seq/")

#load 10X RDS
data_integrated <- readRDS("All_patient_and_healthy_donors_BMMC_PBMC.rds")

# get patient and healthy pbmc data
Idents(data_integrated) <- "sample"
patient_healthy_obj <- subset(data_integrated,idents=c("Patient","Healthy_PBMC"))

# add patient cell source meta data
patient_cells_source_df <- read.csv("Patient_cell_source.meta_data.csv",header=T,row.names=1)

healthy_cells <- grep("rep",colnames(patient_healthy_obj),invert=T,value=T)
healthy_cells_source_df <- data.frame(row.names=healthy_cells,Source=rep("Healthy",10964))

patient_healthy_cells_source_df <- rbind(patient_cells_source_df,healthy_cells_source_df)

# add metadata for patient and healthy source
patient_healthy_obj <- AddMetaData(object = patient_healthy_obj, metadata = patient_healthy_cells_source_df)

cluster_colors <- c("#5A2955", "#D51F26", "#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#DCABCF", "#FDDD03", "#F69421", "#9C82BA", "#AA875A", "#3EB3A7", "#245359", "#711E21", "#B36B45", "#DB7D8D", "#BF5D58", "#216069", "#9A7456", "#8ED2D0", "#3D3D3D", "#B4B883", "#C39C6F", "#9B8EC4", "#6264A0", "#BFCADF", "#8AC972", "#D51F26", "#272E6A", "#208A42")
Idents(patient_healthy_obj) <- "integrated_snn_res.0.6"
levels(patient_healthy_obj) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29") 

pdf("DimPlot_patient_healthy_split_by_source.pdf",height=7,width=25)
DimPlot(patient_healthy_obj,split.by="Source",reduction ="umap",label=F) + scale_color_manual(values=cluster_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

source_colors <- c("#2b83ba","#d7191c","#fdae61","#7570b3")
pdf("DimPlot_patient_healthy_group_by_source.pdf",height=7,width=13)
DimPlot(patient_healthy_obj,split.by="sample",group.by="Source",reduction ="umap",label=F) + scale_color_manual(values=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

saveRDS(patient_healthy_obj,"integrated_patient_healthy_obj.rds")

#######################################################################################################################
#Patient vs. Healthy whole T cells DE
Idents(patient_healthy_obj) <- "integrated_snn_res.0.6"

# T cells cluster(2,3,4,5,14,15,19)
patient_obj_cluster <- subset(patient_healthy_obj, idents=c(2,3,4,5,14,15,19))

pdf("DimPlot_Patient_Healthy_Tcells_group_by_source.pdf",height=8,width=9)
DimPlot(patient_obj_cluster,group.by="Source",reduction ="umap",label=F) + scale_color_manual(values=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

cluster_colors <- c("#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#B36B45", "#DB7D8D", "#8ED2D0", "#208A42")
pdf("DimPlot_Patient_Healthy_Tcells_group_by_cluster.pdf",height=8,width=8.5)
DimPlot(patient_obj_cluster,reduction ="umap",label=F) + scale_color_manual(values=cluster_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

DefaultAssay(patient_obj_cluster) <- "RNA"
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother","Healthy"))

# DE analysis for whole T cells(FindMarkers)
# patient-vs-healthy
dir_Tcells_DE_Patient_Healthy <- paste0("DE_Patient_Healthy_T_cells")
dir.create(dir_Tcells_DE_Patient_Healthy)

FindDE_cluster <- FindMarkers(patient_obj_cluster, ident.1=c("Child","Mother"), ident.2="Healthy", logfc.threshold = 0.25, min.pct = 0.1)
FindDE_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_cluster))

# filPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}

# filtering pval adjust significant DEs(p_val_adj <= 0.05)
FindDE_cluster_res <- FindDE_cluster[FindDE_cluster$p_val_adj <= 0.05,]

output=paste0(dir_Tcells_DE_Patient_Healthy, "/", "DE_Patient_Healthy.Tcells.pValAdj_0.05.csv")
write.csv(FindDE_cluster_res, file=output, quote=F)

# DE visualization
# FeaturePlot
for(i in 1:10){
        pdf(paste0(dir_Tcells_DE_Patient_Healthy, "/FeaturePlot_DE_Patient_Healthy_",rownames(FindDE_cluster_res)[i],".pdf"), width=14, height=7)
        print(FeaturePlot(patient_obj_cluster, features=rownames(FindDE_cluster_res)[i], order=T, split.by="sample", cols=c("grey","red"), max.cutoff="q95", pt.size=0.5, combine=F) + xlim(-14,12) + ylim(-11,18))
        dev.off()
}

# Patient-vs-Healthy average expression scatter plot
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"sample"]
patient_healthy_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}

colnames(patient_healthy_overall_avg_expr) <- c("Patient", "Healthy")

library(ggplot2)
p1 <- ggplot(patient_healthy_overall_avg_expr, aes(Patient, Healthy)) + geom_point(size=2) + ggtitle( paste0("Patient vs. Healthy average expression (T Cells)" )) + theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = patient_healthy_overall_avg_expr[rownames(subset(FindDE_cluster_res,(avg_logFC > 0))),], col="red") + geom_point(data = patient_healthy_overall_avg_expr[rownames(subset(FindDE_cluster_res,(avg_logFC <= 0))),], col="green") + geom_abline(slope=1,intercept=0.25,col="gray",lty=2) + geom_abline(slope=1,intercept=-0.25,col="gray",lty=2)
genes_to_label <- c("CCL4","IFITM1","IFI27","CCL5","CD74","CST7","RGS1","CXCR4","SRGN","LGALS1","FGFBP2","GZMB","NKG7","PRF1","HBA1","HLA-DPB1","LTB","FOSB","FOS","IL7R","SELL","TCF7","CD27","LEF1","MYC","CCR7","KLRB1","TPT1","LYZ","EGR1","FOXP1")
p1 <- LabelPoints(plot = p1, points = genes_to_label, size=5, repel = TRUE)
pdf(file=paste0(dir_Tcells_DE_Patient_Healthy, "/ScatterPlot_Patient_vs_Healthy_AverageExpression.pdf"), width=8, height=8)
print(p1)
dev.off()

Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_healthy_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}

colnames(patient_healthy_overall_avg_expr) <- c("Child", "Mother", "Healthy")

#patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[rownames(FindDE_cluster_res[abs(FindDE_cluster_res$avg_logFC)>1,]),]
select_gene <- c("CCL4","FGFBP2","GZMH","NKG7","GZMB","CCL5","PRF1","CCL4L2","CXCR4","GNLY","LY6E","FCGR3A","CTSW","KLRD1","HLA-DPB1","HLA-DRB1","HLA-DPA1","HLA-DQA2","IFITM1","IFI27","CST7","RGS1","LGALS1","SRGN","TNFAIP3","TSC22D3","DDIT4","IFI6","ZFP36","ZEB2","EFHD2","IFI44L","ISG15","IFITM2","XAF1","LTB","FOSB","FOS","KLRB1","IL7R","TCF7","CD27","SELL","CCR7","LYZ","LEF1","MYC","EGR1","FOXP1","ATM","TLE4","SATB1","NR4A1","BCL2","JUN")
patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[select_gene,]

library(pheatmap)
pdf(paste0(dir_Tcells_DE_Patient_Healthy, "/Heatmap_Patient_vs_Healthy_DE_AverageExpression.pdf"), width=6, height=10)
pheatmap(patient_healthy_DE_avg_expr,scale="row",show_rownames=T,cluster_row=F,border=F,gaps_row=c(18,35),cellwidth=80,cellhight=50,angle_col=0,fontsize=12)
dev.off()

select_gene2 <- c("NKG7","GZMA","GZMB","GNLY","PRF1","IFNG","LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","ICOS","TNFRSF9","TNFRSF14","TNFRSF4","TNFRSF18","CD27","CD28","CD40LG","SELL","CCR7","TCF7","LEF1","FOSB","FOS","NR4A1","BACH2","ID3","FOXO1","BCL6","GATA3","HIF1A","EOMES","RUNX3","TBX21","PRDM1","ID2","ZEB2","STAT4","ZNF683","TOX","TOX2","TOX4")
patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[select_gene2,]

library(pheatmap)
pdf(paste0(dir_Tcells_DE_Patient_Healthy, "/Heatmap_Patient_vs_Healthy_DE_AverageExpression2.pdf"), width=6, height=10)
pheatmap(patient_healthy_DE_avg_expr,scale="row",show_rownames=T,cluster_row=F,border=F,gaps_row=c(6,11,19),cellwidth=80,cellhight=50,angle_col=0,fontsize=12)
dev.off()

# IL2RG related gene
library(ggplot2)
library(viridis)
feature_input_list = c("IL2", "IL2RA","IL2RB","IL2RG","IL4","IL4R","IL7","IL7R","IL9R","IL15","IL15RA","IL21R","JAK1","JAK2","JAK3","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6")
Idents(patient_obj_cluster) <- "Source"
pdf(file=paste0(dir_Tcells_DE_Patient_Healthy, "/DotPlots_group_patient_source_IL2RG_related_gene.pdf"), width=5, height=6)
print(DotPlot(patient_obj_cluster, assay="SCT", group.by="Source", idents=c("Child","Mother","Healthy"), features = feature_input_list, cols=c("green","red"), scale.min=0, scale.max=40) + coord_flip() + scale_color_viridis())
dev.off()

# exhausted related gene
feature_input_list2 = c("PDCD1","HAVCR2","LAG3","CTLA4","TIGIT","TOX","TOX2","TOX4")
Idents(patient_obj_cluster) <- "Source"
pdf(file=paste0(dir_Tcells_DE_Patient_Healthy, "/DotPlots_group_patient_source_exhausted_related_gene.pdf"), width=6, height=4)
print(DotPlot(patient_obj_cluster, assay="SCT", group.by="Source", idents=c("Child","Mother","Healthy"), features = feature_input_list2, cols=c("green","red"), scale.min=0, scale.max=20) + scale_color_viridis(direction=1) + coord_flip())
dev.off()

pdf(file=paste0(dir_Tcells_DE_Patient_Healthy, "/VlnPlots_Patient_Healthy_exhausted_related_gene.pdf"), width=8, height=6)
print(VlnPlot(patient_obj_cluster, assay="SCT", pt.size=0, group.by="Source", ncol=4, idents=c("Child","Mother","Healthy"), features = feature_input_list2, cols=source_colors[-3]))
dev.off()

saveRDS(patient_obj_cluster,paste0(dir_Tcells_DE_Patient_Healthy, "/Patient_Healthy_Tcells_obj.rds"))

#######################################################################################################################
### Patient vs. Healthy NK cells DE
Idents(patient_healthy_obj) <- "integrated_snn_res.0.6"

# NK cells cluster(0)
patient_obj_cluster <- subset(patient_healthy_obj, idents=0)

pdf("DimPlot_Patient_Healthy_NKcells_group_by_source.pdf",height=8,width=9)
DimPlot(patient_obj_cluster,group.by="Source",reduction ="umap",label=F) + scale_color_manual(values=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

DefaultAssay(patient_obj_cluster) <- "RNA"
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother","Healthy"))
patient_obj_cluster

# DE analysis for whole NK cells(FindMarkers)
# patient-vs-healthy
dir_NKcells_DE_Patient_Healthy <- paste0("DE_Patient_Healthy_NK_cells")
dir.create(dir_NKcells_DE_Patient_Healthy)

FindDE_cluster <- FindMarkers(patient_obj_cluster, ident.1=c("Child","Mother"), ident.2="Healthy", logfc.threshold = 0.25, min.pct = 0.1)
FindDE_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_cluster))

# filPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}

# filtering pval adjust significant DEs
FindDE_cluster_res <- FindDE_cluster[FindDE_cluster$p_val_adj <= 0.05,]

output=paste0(dir_NKcells_DE_Patient_Healthy, "/", "DE_Patient_Healthy.NKcells.pValAdj_0.05.csv")
write.csv(FindDE_cluster_res, file=output, quote=F)

# DE visualization
# FeaturePlot
for(i in 1:10){
        pdf(paste0(dir_NKcells_DE_Patient_Healthy, "/FeaturePlot_DE_Patient_Healthy_",rownames(FindDE_cluster_res)[i],".pdf"), width=14, height=7)
        print(FeaturePlot(patient_obj_cluster, features=rownames(FindDE_cluster_res)[i], order=T, split.by="sample", cols=c("grey","red"), max.cutoff="q95", pt.size=0.5) + xlim(-14,12) + ylim(-11,18))
        dev.off()
}

# Patient-vs-Healthy average expression scater plot
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"sample"]
patient_healthy_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}

colnames(patient_healthy_overall_avg_expr) <- c("Patient", "Healthy")

library(ggplot2)
p1 <- ggplot(patient_healthy_overall_avg_expr, aes(Patient, Healthy)) + geom_point(size=2) + ggtitle( paste0("Patient vs. Healthy average expression (NK Cells)" )) + theme_bw(base_size=16) + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = patient_healthy_overall_avg_expr[rownames(subset(FindDE_cluster_res,(avg_logFC > 0))),], col="red") + geom_point(data = patient_healthy_overall_avg_expr[rownames(subset(FindDE_cluster_res,(avg_logFC <= 0))),], col="green") + geom_abline(slope=1,intercept=0.25,col="gray",lty=2) + geom_abline(slope=1,intercept=-0.25,col="gray",lty=2)  
genes_to_label <- c("IFI27","IFITM1","CXCR4","IFI6","IL32","LAIR2","LGALS1","TNFAIP3","LAG3","DDIT4","GZMH","HLA-DPB1","ISG15","CD8B","CD3D","TRAC","HBB","FOS","FOSB","KLRB1","KLRG1","LYZ","IL2RB","IL7R","JUN","LTB","TRDC","FCER1G","GZMK","CXXCR5","CEBPD")
p1 <- LabelPoints(plot = p1, points = genes_to_label, size=5, repel = TRUE)
pdf(file=paste0(dir_NKcells_DE_Patient_Healthy, "/ScatterPlot_Patient_vs_Healthy_AverageExpression.pdf"), width=8, height=8)
print(p1)
dev.off()

Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_healthy_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(patient_healthy_overall_avg_expr), perl=T)
if(length(tmp) > 0){patient_healthy_overall_avg_expr <- patient_healthy_overall_avg_expr[-tmp, ]}

colnames(patient_healthy_overall_avg_expr) <- c("Child", "Mother", "Healthy")

#patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[rownames(FindDE_cluster_res[abs(FindDE_cluster_res$avg_logFC)>1,]),]
select_gene <- c("GZMH","CXCR4","CCL4L2","HLA-DRB1","HLA-DPB1","HLA-DQA2","HLA-DPA1","TRAC","IL32","CCL4","CCL5","TRBC2","CD74","IRF9","LY6E","IFI27","IFITM1","IFI6","ISG15","TNFAIP3","LGALS1","TSC22D3","DDIT4","CST7","IFI44L","LITAF","IFIT3","FOS","KLRB1","FCER1G","FOSB","LYZ","S100A9","TRDC","CMC1","CD7","KLRG1","LTB","KLRC1","IL2RB","JUN","KLRF1","IGFBP7","CEBPD","CXXC5","GZMK","XCL1")
patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[select_gene,]

library(pheatmap)
pdf(paste0(dir_NKcells_DE_Patient_Healthy, "/Heatmap_Patient_vs_Healthy_DE_AverageExpression.pdf"), width=6, height=10)
pheatmap(patient_healthy_DE_avg_expr,scale="row",show_rownames=T,cluster_row=F,gaps_row=c(15,27),border=F,cellwidth=80,cellhight=50,angle_col=0,fontsize=12)
dev.off()

select_gene2 <- c("NKG7","GZMA","GZMB","GNLY","PRF1","IFNG","LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","ICOS","TNFRSF9","TNFRSF14","TNFRSF4","TNFRSF18","CD27","CD28","CD40LG","SELL","CCR7","TCF7","LEF1","FOSB","FOS","NR4A1","BACH2","ID3","FOXO1","BCL6","GATA3","HIF1A","EOMES","RUNX3","TBX21","PRDM1","ID2","ZEB2","STAT4","ZNF683","TOX","TOX2","TOX4")
patient_healthy_DE_avg_expr <- patient_healthy_overall_avg_expr[select_gene2,]

library(pheatmap)
pdf(paste0(dir_NKcells_DE_Patient_Healthy, "/Heatmap_Patient_vs_Healthy_DE_AverageExpression2.pdf"), width=6, height=10)
pheatmap(patient_healthy_DE_avg_expr,scale="row",show_rownames=T,cluster_row=F,border=F,gaps_row=c(6,11,19),cellwidth=80,cellhight=50,angle_col=0,fontsize=12)
dev.off()

# IL2RG related gene
library(ggplot2)
library(viridis)
feature_input_list = c("IL2", "IL2RA","IL2RB","IL2RG","IL4","IL4R","IL7","IL7R","IL9R","IL15","IL15RA","IL21R","JAK1","JAK2","JAK3","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6")
Idents(patient_obj_cluster) <- "Source"
pdf(file=paste0(dir_NKcells_DE_Patient_Healthy, "/DotPlots_group_patient_source_IL2RG_related_gene.pdf"), width=5, height=6)
print(DotPlot(patient_obj_cluster, assay="SCT", group.by="Source", idents=c("Child","Mother","Healthy"), features = feature_input_list, cols=c("green","red"), scale.max=40) + coord_flip() + scale_color_viridis())
dev.off()

# exhausted related gene
feature_input_list2 = c("PDCD1","HAVCR2","LAG3","CTLA4","TIGIT","TOX","TOX2","TOX4")
Idents(patient_obj_cluster) <- "Source"
pdf(file=paste0(dir_NKcells_DE_Patient_Healthy, "/DotPlots_group_patient_source_exhausted_related_gene.pdf"), width=6, height=4)
print(DotPlot(patient_obj_cluster, assay="SCT", group.by="Source", idents=c("Child","Mother","Healthy"), features = feature_input_list2, cols=c("green","red"), scale.max=30) + scale_color_viridis(direction=1) + coord_flip())
dev.off()

pdf(file=paste0(dir_NKcells_DE_Patient_Healthy, "/VlnPlots_Patient_Healthy_exhausted_related_gene.pdf"), width=8, height=6)
print(VlnPlot(patient_obj_cluster, assay="SCT", pt.size=0, group.by="Source", ncol=3, idents=c("Child","Mother","Healthy"), features = feature_input_list2, cols=source_colors[-3]))
dev.off()

saveRDS(patient_obj_cluster,paste0(dir_NKcells_DE_Patient_Healthy, "/Patient_Healthy_NKcells_obj.rds"))
