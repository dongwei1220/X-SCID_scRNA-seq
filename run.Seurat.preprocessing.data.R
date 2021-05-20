###################################################################################################
# Pre-processing X-SCID patient BMMCs scRNA-seq data                                              #
###################################################################################################
rm(list=ls())

# load required packages
library(Seurat)
library(cowplot)
library(ggplot2)

# https://satijalab.org/seurat/v3.0/future_vignette.html
library(future)
# it is tricky to set options
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 30 * 1024^3)

#workflow
#https://github.com/satijalab/seurat/issues/1836
#1. Create Seurat object
#2. QC by filtering out cells based on percent.mito and nFeature_RNA
#3. SCT normalize each dataset specifying the parameter vars.to.regress = percent.mito
#4. Integrate all datasets
#5. Run PCA, UMAP, FindClusters, FindNeighbors (on default assay which is "integrated")
#6. Change default assay to "SCT" to generate FeaturePlots and perform differential expression analysis

########################################################################################################################
# define function for loading 10x data
ReadData10X_Patient <- function(x){
  # rep1
  data_1 <- Read10X(data.dir = "/data/dongwei/projects/X-SCID_scRNA-seq/PatientBMMC_rep1_out/outs/filtered_feature_bc_matrix")
  object_1 <- CreateSeuratObject(counts = data_1, min.cells = 3, min.features = 200)
  # rep2
  data_2 <- Read10X(data.dir = "/data/dongwei/projects/X-SCID_scRNA-seq/PatientBMMC_rep2_out/outs/filtered_feature_bc_matrix")
  object_2 <- CreateSeuratObject(counts = data_2, min.cells = 3, min.features = 200)

 # combine two replicates
 combined <- merge(object_1, y = c(object_2),
                   add.cell.ids = c("rep1", "rep2"),
                   project = "X-SCID")
  combined[["sample"]] <- "Patient"
  combined[["chemistry"]] <- "V3"
  return(combined)
}

# read data
Patient_BMMC <- ReadData10X_Patient(1)
Patient_list <- c(Patient_BMMC)

# QC and SCTransform normalization
for (i in 1:length(Patient_list)) {
  Patient_list[[i]] <- PercentageFeatureSet(Patient_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  Patient_list[[i]] <- PercentageFeatureSet(Patient_list[[i]], pattern = "^RP[SL]", col.name = "percent.ribo")
  # QC plot1
  pdf(paste0("/data/dongwei/projects/X-SCID_scRNA-seq/Patient_BMMC_",i,"_QCPlot1.pdf"),width=12,height=8)
  print(VlnPlot(Patient_list[[i]], features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"), pt.size = 0.1, ncol=4))
  dev.off()

  Patient_list[[i]] <- CellCycleScoring(object = Patient_list[[i]], g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
  # QC plot2
  pdf(paste0("/data/dongwei/projects/X-SCID_scRNA-seq/Patient_BMMC_",i,"_QCPlot2.pdf"),width=10,height=7)
  print(VlnPlot(Patient_list[[i]], features = c("S.Score","G2M.Score"), pt.size = 0.1))
  dev.off()

  # data filtering
  Patient_list[[i]] <- subset(Patient_list[[i]], subset = percent.mt < 10)
  
  # SCTransform normalization
  Patient_list[[i]] <- SCTransform(Patient_list[[i]], vars.to.regress = "percent.mt")
}

saveRDS(Patient_list,"/data/dongwei/projects/X-SCID_scRNA-seq/Patient_BMMC.rds")
