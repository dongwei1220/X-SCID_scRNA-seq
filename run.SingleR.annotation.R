###################################################################################################
# SingleR for cell type annotation                                                                #
###################################################################################################
# https://osca.bioconductor.org/cell-type-annotation.html
rm(list=ls())

# load required packages
library(Seurat)
library(cowplot)
library(ggplot2)

library(SingleR)
library(scater)
library("BiocParallel")
library("SingleCellExperiment")
library(pheatmap)

# https://satijalab.org/seurat/v3.0/future_vignette.html
library(future)
# it is tricky to set options
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 30 * 1024^3)

# set working directory
setwd("/data/dongwei/projects/X-SCID_scRNA-seq/")

################################################################################################################
# run SingleR for cell type annotation
# use NovershternHematopoieticData reference data
################################################################################################################
# load 10X RDS and convert to SingleCellExperiment object
data_integrated <- readRDS("All_patient_and_healthy_donors_BMMC_PBMC.rds")
DefaultAssay(data_integrated) <- "RNA"
Idents(object = data_integrated) <- 'integrated_snn_res.0.6'

# convert to SingleCellExperiment object
data_integrated.sce <- as.SingleCellExperiment(data_integrated)

# load SingleR reference(NovershternHematopoieticData)
# main clusters: 16, fine clusters: 38
ref <- readRDS(file="./SingleR_ref_DB/NovershternHematopoieticData.RDS")
ref
#class: SummarizedExperiment
#dim: 13276 211
#metadata(0):
#assays(1): logcounts
#rownames(13276): 13CDNA73 15E1.2 ... ZZEF1 ZZZ3
#rowData names(0):
#colnames(211): GSM609632 GSM609633 ... GSM609841 GSM609842
#colData names(2): label.main label.fine
table(ref$label.main)
#        B cells       Basophils    CD4+ T cells    CD8+ T cells            CMPs
#             29               6              21              24               4
#Dendritic cells     Eosinophils Erythroid cells            GMPs    Granulocytes
#             10               5              33               4              13
#           HSCs  Megakaryocytes            MEPs       Monocytes        NK cells
#             14              12               9               9              14
#     NK T cells
#              4
length(table(ref$label.fine))
#[1] 38
length(table(ref$label.main))
#[1] 16

# intersect common genes
common <- intersect(rownames(data_integrated.sce), rownames(ref))
ref <- ref[common,]
data_integrated.sce <- data_integrated.sce[common,]

# normalize count
data_integrated.sce <- logNormCounts(data_integrated.sce)

# assign clusters to SingleR object
data_integrated.sce$cluster <- data_integrated@meta.data$integrated_snn_res.0.6

# setting parallel
param <- SnowParam(workers = 20, type = "SOCK")

###### cluster-based annotation #####
## annotation of clusters(for main cluster label)
pred_cluster <- SingleR(test = data_integrated.sce, ref = ref, labels = ref$label.main, 
                        method="cluster", clusters=data_integrated.sce$cluster,
                        de.method="wilcox", BPPARAM=param)

anno_cluster <- data.frame(Cluster=rownames(pred_cluster), Annotation=pred_cluster$pruned.labels)
write.table(anno_cluster, file="SingleR.30_clusters.annotation_main.xls", sep="\t", row.names=F, quote=F, col.names=T)

## annotation of clusters(for fine cluster label)
pred_cluster_fine <- SingleR(test = data_integrated.sce, ref = ref, labels = ref$label.fine, 
                             method="cluster", clusters=data_integrated.sce$cluster,
                             de.method="wilcox", BPPARAM=param)

anno_cluster_fine <- data.frame(Cluster=rownames(pred_cluster_fine), Annotation=pred_cluster_fine$pruned.labels)
write.table(anno_cluster_fine, file="SingleR.30_clusters.annotation_fine.xls", sep="\t", row.names=F, quote=F, col.names=T)

##### annotation of individual cells #####
pred_fine <- SingleR(test = data_integrated.sce, ref = ref, labels = ref$label.fine, de.method="wilcox", BPPARAM=param)
table(pred_fine$labels)

cell_labels <- data.frame(Cell=pred_fine@rownames, Label=pred_fine$pruned.labels)
cell_meta <- data.frame(Cell=colnames(data_integrated), Sample=data_integrated@meta.data$sample, Cluster=Idents(data_integrated))
#check order
table(ifelse(cell_labels$Cell == rownames(cell_meta), "Yes", "No"))

cell_labels$Sample <- cell_meta$Sample
cell_labels$Cluster <- cell_meta$Cluster
write.table(cell_labels, file="SingleR.30_clusters.annotation_fine.individual_cell_infor.NovershternHematopoieticData.xls", sep="\t", row.names=F, quote=F, col.names=T)

# cell cluster number
cell_cluster_num <- table(cell_labels$Label,cell_labels$Cluster)
head(cell_cluster_num)
write.table(cell_cluster_num, file="SingleR.30_clusters.annotation_fine.individual_cell_infor.NovershternHematopoieticData_num.xls",sep="\t", row.names=T, quote=F, col.names=T)

# cell cluster percent
cell_cluster_per <- prop.table(cell_cluster_num, 2)
head(cell_cluster_per)
write.table(cell_cluster_per, file="SingleR.30_clusters.annotation_fine.individual_cell_infor.NovershternHematopoieticData_num_per.xls",sep="\t", row.names=T, quote=F, col.names=T)

#######################################################################################################################
library("RColorBrewer")

tab <- read.csv(file="SingleR.30_clusters.annotation_fine.individual_cell_infor.NovershternHematopoieticData.per.ordered.csv", header=T, row.names=1, check.names=F)
head(tab)
#colnames(tab) <- c("cell_type", paste0("cluster_", 0:29) )
tab_matrix <- as.matrix(tab)

jBrewColors_1 <- c("#E31A1C","#FFD700","#771122","#777711","#1F78B4","#68228B","#AAAA44",
                 "#60CC52","#771155","#DDDD77","#774411","#AA7744","#AA4455","#117744",
                 "#000080","#44AA77","#AA4488","#DDAA77")
jBrewColors_2 <- brewer.pal(n = 12, name = "Set3")

#jBrewColors <- c(brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 8, name = "Dark2"))
jBrewColors <-c(jBrewColors_1, jBrewColors_2)

# barplot visualization
pdf("SingleR.30_clusters.annotation_fine.individual_cell_infor.NovershternHematopoieticData_num_per.barplot.pdf",height=8,width=13)
barplot(tab_matrix, col=jBrewColors)
plot(1, type = "n", axes = FALSE, ann = FALSE)
legend("center", legend=rownames(tab_matrix), text.font=1, ncol=3, fill=jBrewColors)
dev.off()
