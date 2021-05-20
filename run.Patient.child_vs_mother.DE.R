###################################################################################################
# DE analysis for patient child-vs-mother sample                                                  #
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

#load 10X RDS file
data_integrated <- readRDS("All_patient_and_healthy_donors_BMMC_PBMC.rds")

# get patient data
patient_obj <- subset(data_integrated,subset=sample=="Patient")
patient_obj

# add patient cell source meta data
patient_cells_source_df <- read.csv("Patient_cell_source.meta_data.csv",header=T,row.names=1)
patient_obj <- AddMetaData(object = patient_obj,metadata = patient_cells_source_df)
table(patient_obj$Source)

patient_stat <- table(Idents(patient_obj),patient_obj$Source)

library(reshape2)
library(ggplot2)
patient_stat <- melt(patient_stat)
colnames(patient_stat) <- c("Cluster","Sample","Value")
patient_stat$Cluster <- factor(patient_stat$Cluster,levels=rev(c(15,0,14,3,4,29,5,19,2,28,26,23,27,10,8,1,11,22,18,16,17,25,20,9,6,13,7,24,21,12)))

source_colors <- c("#2b83ba","#d7191c","#fdae61")
pdf("BarPlot_patient_source_stat.pdf",height=9,width=6)
ggplot(patient_stat,aes(Cluster,Value,fill=Sample)) + geom_bar(stat="identity",position= "fill") + theme_bw() + coord_flip() + scale_fill_manual(values=source_colors)
dev.off()

pdf("DimPlot_patient_group_by_source.pdf",height=8,width=9)
DimPlot(patient_obj,group.by="Source",reduction ="umap",label=F,cols=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

pdf("DimPlot_patient_split_by_source.pdf",height=8,width=23)
DimPlot(patient_obj,group.by="Source",split.by="Source",reduction ="umap",label=F,cols=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

cluster_colors <- c("#5A2955", "#D51F26", "#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#DCABCF", "#FDDD03", "#F69421", "#9C82BA", "#AA875A", "#3EB3A7", "#245359", "#711E21", "#B36B45", "#DB7D8D", "#BF5D58", "#216069", "#9A7456", "#8ED2D0", "#3D3D3D", "#B4B883", "#C39C6F", "#9B8EC4", "#6264A0", "#BFCADF", "#8AC972", "#D51F26", "#272E6A", "#208A42")
pdf("DimPlot_patient_group_by_cluster.pdf",height=8,width=9.5)
DimPlot(patient_obj,reduction ="umap",label=F) + scale_color_manual(values=cluster_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

saveRDS(patient_obj,"integrated_patient_obj.rds")

#######################################################################################################################
### patient Child vs. Mother whole T cells DE
Idents(patient_obj) <- "integrated_snn_res.0.6"

# T cells cluster(2,3,4,5,14,15,19)
patient_obj_cluster <- subset(patient_obj, idents=c(2,3,4,5,14,15,19))

DefaultAssay(patient_obj_cluster) <- "RNA"
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother"))

pdf("DimPlot_Patient_Tcells_group_by_source.pdf",height=8,width=9)
DimPlot(patient_obj_cluster,group.by="Source",reduction ="umap",label=F,cols=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

levels(patient_obj_cluster) <- c("2","3","4","5","14","15","19")
cluster_colors <- c("#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#B36B45", "#DB7D8D", "#8ED2D0", "#208A42")
pdf("DimPlot_Patient_Tcells_group_by_cluster.pdf",height=8,width=8.5)
DimPlot(patient_obj_cluster,reduction ="umap",label=F) + scale_color_manual(values=cluster_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

# DE analysis for whole T cells(FindMarkers)
# child-vs-mother
dir_Tcells_DE_Child_Mother <- paste0("DE_Patient_Child_Mother_T_cells")
dir.create(dir_Tcells_DE_Child_Mother)

FindDE_cluster <- FindMarkers(patient_obj_cluster, ident.1="Child", ident.2="Mother", test.use = "wilcox", assay="RNA", slot="data", logfc.threshold = 0.25, min.pct = 0.1)
FindDE_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_cluster))

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}

output=paste0(dir_Tcells_DE_Child_Mother, "/", "DE_Child_Mother.Tcells.csv")
write.csv(FindDE_cluster, file=output, quote=F)

# DE visualization
pdf(file=paste0(dir_Tcells_DE_Child_Mother, "/VlnPlots_DE_Child_Mother.pdf"), width=16, height=8)
print(VlnPlot(patient_obj_cluster, assay="SCT", pt.size=0, group.by="Source", ncol=5, idents=c("Child","Mother"), features = rownames(FindDE_cluster), cols=source_colors[1:2]))
dev.off()

# Child-vs-Mother average expression scatter plot
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
child_mother_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(child_mother_overall_avg_expr), perl=T)
if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(child_mother_overall_avg_expr), perl=T)
if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}

colnames(child_mother_overall_avg_expr) <- c("Child", "Mother")

library(ggplot2)
p1 <- ggplot(child_mother_overall_avg_expr, aes(Child, Mother)) + geom_point(size=2) + ggtitle( paste0("Child vs. Mother average expression (T Cells)" )) + theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC > 0))),], col="red") + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC <= 0))),], col="green") + geom_abline(slope=1,intercept=0.25,col="gray",lty=2) + geom_abline(slope=1,intercept=-0.25,col="gray",lty=2)
genes_to_label <- rownames(FindDE_cluster)
p1 <- LabelPoints(plot = p1, points = genes_to_label, size= 5, repel = TRUE)
pdf(file=paste0(dir_Tcells_DE_Child_Mother, "/ScatterPlot_Child_vs_Mother_AverageExpression.pdf"), width=9, height=9)
print(p1)
dev.off()

###################################################################################################
### patient Child vs. Mother NK cells DE
Idents(patient_obj) <- "integrated_snn_res.0.6"

# NK cells cluster(0)
patient_obj_cluster <- subset(patient_obj, idents=0)

DefaultAssay(patient_obj_cluster) <- "RNA"
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother"))

pdf("DimPlot_Patient_NKcells_group_by_source.pdf",height=8,width=9)
DimPlot(patient_obj_cluster,group.by="Source",reduction ="umap",label=F, cols=source_colors) + xlim(-14,12) + ylim(-11,18)
dev.off()

# DE analysis for NK cells(FindMarkers)
# child-vs-mother
dir_NKcells_DE_Child_Mother <- paste0("DE_Patient_Child_Mother_NK_cells")
dir.create(dir_NKcells_DE_Child_Mother)

FindDE_cluster <- FindMarkers(patient_obj_cluster, ident.1="Child", ident.2="Mother", test.use = "wilcox", assay="RNA", slot="data", logfc.threshold = 0.25, min.pct = 0.1)
FindDE_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_cluster))

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}

output=paste0(dir_NKcells_DE_Child_Mother, "/", "DE_Child_Mother.NKcells.csv")
write.csv(FindDE_cluster, file=output, quote=F)

# DE visualization
pdf(file=paste0(dir_NKcells_DE_Child_Mother, "/VlnPlots_DE_Child_Mother.pdf"), width=16, height=8)
print(VlnPlot(patient_obj_cluster, assay="SCT", pt.size=0, group.by="Source", ncol=5, idents=c("Child","Mother"), features = rownames(FindDE_cluster), cols=source_colors[1:2]))
dev.off()

# Child-vs-Mother average expression scatter plot
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
child_mother_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[L|S]",rownames(child_mother_overall_avg_expr), perl=T)
if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}
tmp <- grep("^MT-",rownames(child_mother_overall_avg_expr), perl=T)
if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}

colnames(child_mother_overall_avg_expr) <- c("Child", "Mother")

library(ggplot2)
p1 <- ggplot(child_mother_overall_avg_expr, aes(Child, Mother)) + geom_point(size=2) + ggtitle( paste0("Child vs. Mother average expression (NK Cells)" )) + theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC > 0))),], col="red") + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC <= 0))),], col="green") + geom_abline(slope=1,intercept=0.25,col="gray",lty=2) + geom_abline(slope=1,intercept=-0.25,col="gray",lty=2)
genes_to_label <- rownames(FindDE_cluster)
p1 <- LabelPoints(plot = p1, points = genes_to_label, size=5, repel = TRUE)
pdf(file=paste0(dir_NKcells_DE_Child_Mother, "/ScatterPlot_Child_vs_Mother_AverageExpression.pdf"), width=9, height=9)
print(p1)
dev.off()

#############################################################################################################
# patient Child vs. Mother for selected cluster DE
Idents(patient_obj) <- "integrated_snn_res.0.6"

# T and NK cell cluster(5,2,4,3,15,19,14,0)
patient_obj_cluster <- subset(patient_obj, idents=c(0,2,3,4,5,14,15,19))
patient_obj_cluster

DefaultAssay(patient_obj_cluster) <- "RNA"
Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother"))

# T and NK cell cluster(5,2,4,3,15,19,14,0)
dir_DE_Child_Mother <- paste0("/data/sharing/Davey/CVID_PatientYang/DE_Patient_Child_Mother_Cluster")
dir.create(dir_DE_Child_Mother)

# DE analysis for each cluster
for (i in c(0,2,3,4,5,14,15,19)){
        patient_obj_cluster <- subset(patient_obj, subset=integrated_snn_res.0.6==i)
        DefaultAssay(patient_obj_cluster) <- "RNA"
        Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
        patient_obj_cluster <- subset(patient_obj_cluster,idents=c("Child","Mother"))
 
        FindDE_cluster <- FindMarkers(patient_obj_cluster, ident.1="Child", ident.2="Mother", test.use = "wilcox", assay="RNA", slot="data", logfc.threshold = 0.25, min.pct = 0.1)
        FindDE_cluster$cluster <- paste0("cluster_", i)
        FindDE_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_cluster))

# filter RPLXX (Ribosomal protein)
        tmp <- grep("^RP[L|S]",rownames(FindDE_cluster), perl=T)
        if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
        tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
        if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
# filter p val adjust
        FindDE_cluster <- FindDE_cluster[FindDE_cluster$p_val < 0.05, ]

        output=paste0(dir_DE_Child_Mother, "/", "DE_Child_Mother.cluster_", i, ".pVal_0.05.csv")
        write.csv(FindDE_cluster, file=output, quote=F)

        # top10 DE VlnPlot visualization
        pdf(file=paste0(dir_DE_Child_Mother, "/VlnPlot_DE_Child_Mother.cluster_", i, "_Top10.pdf"), width=16, height=8)
        print(VlnPlot(patient_obj_cluster, assay="SCT", pt.size=0, group.by="Source", idents=c("Child","Mother"), ncol=5, features = rownames(FindDE_cluster)[1:10]))
        dev.off()

        # DotPlot
        library(viridis)
        pdf(file=paste0(dir_DE_Child_Mother, "/DotPlots_DE_Child_Mother.cluster_",i,"_Top10.pdf"), width=15, height=3.5)
        print(DotPlot(patient_obj_cluster, assay="RNA", group.by="Source", idents=c("Child","Mother"), features = rownames(FindDE_cluster)[1:10], cols=c("green","red")) + scale_color_viridis(direction=1))
        dev.off()

        # Child-vs-Mother average expression scater plot
        Idents(patient_obj_cluster) <- patient_obj_cluster@meta.data[,"Source"]
        child_mother_overall_avg_expr <- log1p(AverageExpression(patient_obj_cluster, verbose=F)$RNA)

# filter RPLXX (Ribosomal protein)
        tmp <- grep("^RP[L|S]",rownames(child_mother_overall_avg_expr), perl=T)
        if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}
        tmp <- grep("^MT-",rownames(child_mother_overall_avg_expr), perl=T)
        if(length(tmp) > 0){child_mother_overall_avg_expr <- child_mother_overall_avg_expr[-tmp, ]}
        colnames(child_mother_overall_avg_expr) <- c("Child", "Mother")

        top10_up <- FindDE_cluster[order(FindDE_cluster$avg_logFC,decreasing=T),][1:10,]
        top10_down <- FindDE_cluster[order(FindDE_cluster$avg_logFC),][1:10,]
        top20 <- c(rownames(top10_up),rownames(top10_down))

        library(ggplot2)
        p <- ggplot(child_mother_overall_avg_expr, aes(Child, Mother)) + geom_point(size=2) + ggtitle( paste0("Child vs. Mother average expression (Cluster",i,")" )) + theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC > 0))),], col="red") + geom_point(data = child_mother_overall_avg_expr[rownames(subset(FindDE_cluster,(avg_logFC <= 0))),], col="green") + geom_abline(slope=1,intercept=0.25,col="gray",lty=2) + geom_abline(slope=1,intercept=-0.25,col="gray",lty=2)
        p <- LabelPoints(plot = p, points = top20, size= 5, repel = TRUE, max.overlaps = 20)
        pdf(file=paste0(dir_DE_Child_Mother, "/ScatterPlot_Child_vs_Mother_AverageExpression.cluster_",i,".pdf"), width=9, height=9)
        print(p)
        dev.off()
}

