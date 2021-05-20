###################################################################################################
# single cell geneset enrichment                                                                  #
###################################################################################################
# AUCell
# https://www.bioconductor.org/packages/release/bioc/html/AUCell.html

# load required packages
library("AUCell")
library("GSEABase")
library("NMF")
library("Seurat")

# set working directory
setwd("/data/dongwei/projects/X-SCID_scRNA-seq/")

patient_healthy_obj <- readRDS("integrated_patient_healthy_obj.rds")

# get T cells subset
patient_healthy_Tcells <- subset(patient_healthy_obj, idents=c(2,3,4,5,14,15,19))

# remove other categray
patient_healthy_Tcells <- subset(patient_healthy_Tcells,subset=Source!="Other")

exprMatrix <- as.matrix(patient_healthy_Tcells[['RNA']]@data)

cellsUmap <- patient_healthy_Tcells$umap@cell.embeddings

pdf("AUCell.patient_healthy.Tcell_cells_rankings.pdf")
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 20)
dev.off()

# Define function for preparing genesets
Prepare_GeneSets <- function(entrez_list, geneSet_names){

	geneSets <- c()

	for (i in 1:length(geneSet_names)){
		entrez_id <- entrez_list[names(entrez_list) %in% geneSet_names[i]]
		entrez_id <- GeneSet(entrez_id[[1]], setName=geneSet_names[i])

		geneIdType(entrez_id) <- EntrezIdentifier()
		geneIds <- geneIds(mapIdentifiers(entrez_id, SymbolIdentifier(annotation="org.Hs.eg")))

		geneSet <- GSEABase::GeneSet(geneIds, setName=geneSet_names[i])
		geneSets <- c(geneSets, geneSet)
	}
	names(geneSets) <- geneSet_names
	return(geneSets)
}


Hs_c5_bp <- readRDS("/data/public/GSEA_dataset/Hs.c5.bp.v7.1.entrez.rds")

target <- "_T_cell"
geneSet_names <- names(Hs_c5_bp)[grep(target, names(Hs_c5_bp), ignore.case=T)]

dir_AUCellPlot <- paste0("AUCell_UMAP_Tcells", target)
dir.create(dir_AUCellPlot)

geneSets <- Prepare_GeneSets(Hs_c5_bp, geneSet_names)

geneSets <- GeneSetCollection(geneSets)

geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 

# calculate AUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05, nCores = 20)

set.seed(123)
pdf("AUCell.patient_healthy.Tcells.cells_assignment.pdf")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=20, assign=TRUE)
dev.off()

selectedThresholds <- getThresholdSelected(cells_assignment)

#Exploring the cell-assignment (table & heatmap)
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])

set.seed(123)
miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]

library(NMF)
assignmentMat_anno <- patient_healthy_Tcells@meta.data[rownames(patient_healthy_Tcells@meta.data) %in% colnames(assignmentMat), "Source"]

pdf("AUCell.patient_healthy.Tcells.cells.assignmentMat.pdf",height=10,width=15,onefile=F)
aheatmap(assignmentMat, annCol=assignmentMat_anno, scale="none", color="black", legend=FALSE)
dev.off()

########################################################################
selectedThresholds <- getThresholdSelected(cells_assignment)
selectedThresholds
#target_GO <- "GO_T_cell_proliferation_involved_in_immune_response"
#selectedThresholds <- selectedThresholds[names(selectedThresholds) == target_GO]

library(viridis)
for(geneSetName in names(selectedThresholds))
{
  print(paste0("processing ", geneSetName))
  nBreaks <- 100

  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    

	#cellsUmap_df <- data.frame(UMAP_1=cellsUmap[,1], UMAP_2=cellsUmap[,2])
	cellsUmap_df <- data.frame(UMAP_1=cellsUmap[,1], UMAP_2=cellsUmap[,2], Source=patient_healthy_Tcells$Source)
	cellsUmap_df$auc <- min(aucSplit[[1]], aucSplit[[2]])
	cellsUmap_df[rownames(cellsUmap_df) %in% names(aucSplit[[2]]),]$auc <- aucSplit[[2]]

	cutoff_min <- quantile(aucSplit[[2]], probs=seq(0,1,0.05))[[2]]
	cutoff_max <- quantile(aucSplit[[2]], probs=seq(0,1,0.05))[[20]]

	cellsUmap_df$auc[cellsUmap_df$auc < cutoff_min] <- cutoff_min
	cellsUmap_df$auc[cellsUmap_df$auc > cutoff_max] <- cutoff_max

	cellsUmap_sub <- cellsUmap[rownames(cellsUmap) %in% names(aucSplit[[2]][aucSplit[[2]] >= cutoff_max]), ]

#  	plot_colour <- c("lightgrey", brewer.pal(n = 8, name = "Oranges"))[1:9]
#  	plot_colour <- c("lightgrey", brewer.pal(n = 9, name = "YlGnBu"))[1:9]
#  	plot_colour <- c("lightgrey", brewer.pal(8, "Spectral"))[1:9]
    plot_colour <- colorRampPalette(c("gray","orange","red"))(500)

  	pp <- ggplot(data=cellsUmap_df) + 
  		geom_point(aes(x=UMAP_1, y=UMAP_2, color=auc, fill=auc), alpha=1, size=0.6) + 
  		geom_point(data=cellsUmap_df[rownames(cellsUmap_df) %in% rownames(cellsUmap_sub), ], aes(x=UMAP_1, y=UMAP_2, color=auc, fill=auc), alpha=1, size=0.8) +
  		scale_color_gradientn(colors=plot_colour) + 
  		scale_fill_gradientn(colors=plot_colour) +
  		#scale_color_viridis(option="B") +
  		#scale_fill_viridis(option="B") + 
  		xlim(-14,12) + ylim(-11,18) +
  		theme_bw() + facet_wrap(.~Source,nrow=1) +
  		ggtitle(geneSetName) +
  		theme(legend.position = "none", panel.grid.major=element_blank(), panel.grid.minor=element_blank())

   pdf(paste0(dir_AUCellPlot, "/", geneSetName, ".pdf"), height=7,width=19)
   print(pp)
   dev.off()

  }
}
