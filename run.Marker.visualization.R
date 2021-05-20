###################################################################################################
# Cell type marker gene visualization                                                             #
###################################################################################################

# load required packages
library(Seurat)
library(cowplot)
library(ggplot2)

# https://satijalab.org/seurat/v3.0/future_vignette.html
library(future)
# it is tricky to set options
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 30 * 1024^3)

# set working directory
setwd("/data/dongwei/projects/X-SCID_scRNA-seq/")

########################################################################################################################
#load 10X RDS
data_integrated <- readRDS("All_patient_and_healthy_donors_BMMC_PBMC.rds")

# change default assay to "SCT" for marker gene visualization
DefaultAssay(data_integrated) <- "SCT"

Idents(data_integrated) <- "integrated_snn_res.0.6"
levels(data_integrated) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29")

### cell type marker gene visualization
feature_input_list = c("CD34","CD164","HLF","CRHBP","ADGRG6","PCDH9","DUSP1","CD14","LYZ","MS4A7","MS4A6A","ANXA2","PPBP","PF4","GP9","MME","EBF1","PAX5","ITGA2B","PLEK","VWF","MPO","ELANE","CTSG","PRTN3","CSF1R","HBG1","HBG2","GATA1","KLF1","CA1","TFR2","HBB","CD3D","CD3E","CD4","CXCR3","CXCR5","IL2RA","CCR6","FOXP3","TNFRSF18","CD44","CD8A","CD8B","SELL","CCR7","IL7R","KLRG1","GZMA","GZMB","GZMK","PRF1","S100A4","MS4A1","IGHD","IGHM","IGHA1","IGHA2","IGHG1","IGHG3","CD19","CD27","CD79A","CD79B","CD22","FCER2","FCGR3A","FCGR3B","GNLY","NKG7","FCER1A","CST3","MPEG1","IRF8","SPIB","CCR2","IGKC","CD244","MKI67")

# Feature plot visualization
dir_FeaturePlot <- paste0("FeaturePlot_All_CellType_Markers")
dir.create(dir_FeaturePlot)
for (i in 1:length(feature_input_list)){
        feature_input <- feature_input_list[i]
        pdf(file=paste0(dir_FeaturePlot, "/FeaturePlots_", feature_input, "_ClusterPlot.pdf"), width=14, height=14)
        print(FeaturePlot(data_integrated, slot="data", features = feature_input, max.cutoff = "q95", pt.size=0.5, cols = c("grey", "red"), label=F, order=F) + xlim(-14,12) + ylim(-11,18))
        dev.off()
        png(file=paste0(dir_FeaturePlot, "/FeaturePlots_", feature_input, "_ClusterPlot.png"), width=800, height=800)
        print(FeaturePlot(data_integrated, slot="data", features = feature_input, max.cutoff = "q95", pt.size=0.5, cols = c("grey", "red"), label=F, order=F) + xlim(-14,12) + ylim(-11,18))
        dev.off()

}

# Dot plot visualization
library(ggplot2)
library(viridis)

HSCmarks <- c("CD34","CRHBP")
Erymarks <- c("GATA1","KLF1","CA1","HBB")
Platelet <- c("PF4","ITGA2B","PPBP")
Granmarks <- c("MPO","ELANE","CTSG","PRTN3")
Monomarks <- c("CD14","LYZ","MS4A7","MS4A6A","ANXA2","FCGR3A")
DCmarks <- c("FCER1A","CST3","MPEG1","IGKC")
Bcellmarks <- c("CD79A","CD79B","MS4A1","IGHD","IGHM","IGHA1","IGHG1")
Tcellmarks <- c("CD3D","CD3E","CD4","IL7R","CCR7","FOXP3","CD8A","CD8B","SELL","KLRG1","GZMA","GZMB","GZMK")
NKmarks <- c("GNLY","NKG7")

features <- list("HSC/CMP" = HSCmarks, "Erythrocytes" = Erymarks, "Megakaryocytes/Platelet" = Platelet, "Granulocytes" = Granmarks, "Monocytes" = Monomarks, "DC" = DCmarks, "B Cells" = Bcellmarks, "T Cells" = Tcellmarks, "NK Cells" = NKmarks)
pdf("Integrated_DotPlot_All_cluster_markers.pdf",height=7,width=15)
levels(data_integrated) <- c("13", "25", "6", "17", "26", "16", "9", "1", "11", "7", "20", "22", "18", "27", "23", "24", "12", "10", "8", "28", "21", "2", "4", "19", "5", "15", "3", "14", "29", "0")
DotPlot(object = data_integrated, features=features) + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

#######################################################################################################################
##### DE analysis for all clusters #####
# FindAllMarkers for DE analysis
DefaultAssay(data_integrated) <- "RNA"

dir_DE_All_Cluster <- paste0("DE_All_Cluster")
dir.create(dir_DE_All_Cluster)

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_integrated, only.pos = T, logfc.threshold = 1, min.pct = 0.25)

FindDE_all_cluster$cluster <- paste0("cluster_", FindDE_all_cluster$cluster)
FindDE_all_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_all_cluster))

# filtering p val adjust
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.01, ]

# filtering RPLXX (Ribosomal protein)
# MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

output=paste0(dir_DE_All_Cluster, "/", "DE_all.cluster.logFC_1_pValAdj_0.01.csv")
write.csv(FindDE_all_cluster, file=output, quote=F)

# DE gene visualization
DefaultAssay(data_integrated) <- "SCT"
Idents(data_integrated) <- "integrated_snn_res.0.6"
levels(data_integrated) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
"21","22","23","24","25","26","27","28","29")

library(dplyr)
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste0(dir_DE_All_Cluster, "/", "DE_all.cluster.top5.heatmap.pdf"), width=30, height=18)
p <- DoHeatmap(data_integrated, cells=1:3000, features = as.vector(top5$gene), group.colors=cluster_colors) + scale_fill_viridis()
print(p)
dev.off()

#######################################################################################################################
# Visualize selected gene expression in each sample
DefaultAssay(data_integrated) <- "SCT"

obj.list <- SplitObject(data_integrated, split.by = "sample")

sample_name <- names(obj.list)

# FeaturePlot Visualization
feature_input_list = c("XIST", "IL2", "IL2RA","IL2RB","IL2RG","IL4","IL4R","IL7","IL7R","IL9R","IL15","IL15RA","IL21R","JAK1","JAK2","JAK3","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","FOS","FOSB","JUN","JUNB","JUND","MYC","BCL2")

for (i in 1:length(feature_input_list)){
	feature_input <- feature_input_list[i]
	dir_FeaturePlot <- paste0("FeaturePlot_", feature_input)
	dir.create(dir_FeaturePlot)

	for (i in 1:10){
		pdf(file=paste0(dir_FeaturePlot, "/FeaturePlots_", feature_input, "_ClusterPlot_", sample_name[i], ".pdf"), width=10, height=10)
		print(FeaturePlot(obj.list[i][[1]], slot="data", features = feature_input, max.cutoff = "q95", cols = c("grey", "red"), label=F, order=T) + xlim(-14,12) + ylim(-11,18) + labs(title = sample_name[i]))
		dev.off()
		png(file=paste0(dir_FeaturePlot, "/FeaturePlots_", feature_input, "_ClusterPlot_", sample_name[i], ".png"), width=800, height=800)
		print(FeaturePlot(obj.list[i][[1]], slot="data", features = feature_input, max.cutoff = "q95", cols = c("grey", "red"), label=F, order=T) + xlim(-14,12) + ylim(-11,18) + labs(title = sample_name[i]))
		dev.off()
	}
}
