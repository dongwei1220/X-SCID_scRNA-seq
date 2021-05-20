###################################################################################################
# GO and KEGG enrichment for DEs                                                                  #
###################################################################################################

# load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# set working directory
dirRoot <- "/data/dongwei/projects/X-SCID_scRNA-seq/")

#######################################################################################################################
# Patient-vs-Healthy T cells DE enrichment
de_genes <- read.csv(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Healthy.Tcells.pValAdj_0.05.csv"),row.names=1)

### patient up gene ###
patient_up <- rownames(de_genes[de_genes$avg_logFC > 0,])

FCgenelist <- de_genes[de_genes$avg_logFC > 0,]$avg_logFC
names(FCgenelist) <- patient_up
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)

# 转换ID
patient_up = bitr(patient_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

### GO Enrichment
eggo <- enrichGO(patient_up$ENTREZID, 
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 ont = "BP", readable = T,
                 pvalueCutoff = 0.01, 
                 qvalueCutoff = 0.05,
                 pAdjustMethod = "BH")

eggosimp <- simplify(eggo,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")

write.csv(eggosimp@result,paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.csv"))

# data visualization
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.barplot.pdf"),height=8,width=10)
barplot(eggosimp,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.dotplot.pdf"),height=8,width=10)
dotplot(eggosimp, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.cnetplot.pdf"),height=8,width=12.5)
cnetplot(eggosimp,circular=F, colorEdge =T, showCategory=5, foldChange=FCgenelist)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(eggosimp)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.GO_enrichment.heatplot.pdf"),height=6,width=12)
heatplot(eggosimp,foldChange=FCgenelist,showCategory=10)
dev.off()

### KEGG Enrichment
egkegg <- enrichKEGG(patient_up$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)
egkegg <- setReadable(egkegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(egkegg@result,paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.KEGG_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.KEGG_enrichment.barplot.pdf"),height=8,width=8)
barplot(egkegg,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.KEGG_enrichment.dotplot.pdf"),height=8,width=8)
dotplot(egkegg, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.KEGG_enrichment.cnetplot.pdf"),height=8,width=10)
cnetplot(egkegg,circular=T, colorEdge =T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Patient_Up.Tcells.KEGG_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(egkegg)
dev.off()

### healthy up gene ###
healthy_up <- rownames(de_genes[de_genes$avg_logFC < 0,])

FCgenelist <- abs(de_genes[de_genes$avg_logFC < 0,]$avg_logFC)
names(FCgenelist) <- healthy_up
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)

# 转换ID
healthy_up = bitr(healthy_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

### GO Enrichment
eggo2 <- enrichGO(healthy_up$ENTREZID, 
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont = "BP", readable = T,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  pAdjustMethod = "BH")
				  
eggosimp2 <- simplify(eggo2,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")

write.csv(eggosimp2@result,paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.barplot.pdf"),height=8,width=9)
barplot(eggosimp2,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.dotplot.pdf"),height=8,width=9)
dotplot(eggosimp2,showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.cnetplot.pdf"),height=9,width=12.5)
cnetplot(eggosimp2,circular=F, foldChange=FCgenelist, colorEdge =T, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(eggosimp2)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.GO_enrichment.heatplot.pdf"),height=6,width=12)
heatplot(eggosimp2,foldChange=FCgenelist,showCategory=10)
dev.off()

### KEGG Enrichment
egkegg2 <- enrichKEGG(healthy_up$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
egkegg2 <- setReadable(egkegg2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(egkegg2@result,paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.KEGG_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.KEGG_enrichment.barplot.pdf"),height=8,width=8)
barplot(egkegg2,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.KEGG_enrichment.dotplot.pdf"),height=8,width=8)
dotplot(egkegg2,showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.KEGG_enrichment.cnetplot.pdf"),height=8,width=10)
cnetplot(egkegg2,circular=T, colorEdge =T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_T_cells/DE_Healthy_Up.Tcells.KEGG_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(egkegg2)
dev.off()

#######################################################################################################################
# Patient-vs-Healthy NK cells DE enrichment
de_genes <- read.csv(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Healthy.NKcells.pValAdj_0.05.csv"),row.names=1)

### patient up gene ###
patient_up <- rownames(de_genes[de_genes$avg_logFC > 0,])

FCgenelist <- de_genes[de_genes$avg_logFC > 0,]$avg_logFC
names(FCgenelist) <- patient_up
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order

# 转换ID
patient_up = bitr(patient_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

### GO Enrichment
eggo <- enrichGO(patient_up$ENTREZID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 ont = "BP", readable = T,
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 pAdjustMethod = "BH")

eggosimp <- simplify(eggo,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")

write.csv(eggosimp@result,paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.csv"))

# data visualization
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.barplot.pdf"),height=8,width=10)
barplot(eggosimp,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.dotplot.pdf"),height=8,width=10)
dotplot(eggosimp, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.cnetplot.pdf"),height=8,width=11)
cnetplot(eggosimp,circular=F, foldChange=FCgenelist, colorEdge =T,showCategory=5)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(eggosimp)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.GO_enrichment.heatplot.pdf"),height=6,width=12)
heatplot(eggosimp,foldChange=FCgenelist,showCategory=10)
dev.off()

### KEGG Enrichment
egkegg <- enrichKEGG(patient_up$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)
egkegg <- setReadable(egkegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(egkegg@result,paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.KEGG_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.KEGG_enrichment.barplot.pdf"),height=8,width=8)
barplot(egkegg,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.KEGG_enrichment.dotplot.pdf"),height=8,width=8)
dotplot(egkegg, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.KEGG_enrichment.cnetplot.pdf"),height=8,width=10)
cnetplot(egkegg,circular=T, colorEdge =T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Patient_Up.NKcells.KEGG_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(egkegg)
dev.off()

### healthy up gene ###
healthy_up <- rownames(de_genes[de_genes$avg_logFC < 0,])

FCgenelist <- abs(de_genes[de_genes$avg_logFC < 0,]$avg_logFC)
names(FCgenelist) <- healthy_up
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)

# 转换ID
healthy_up = bitr(healthy_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

### GO Enrichment
eggo2 <- enrichGO(healthy_up$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont = "BP", readable = T,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  pAdjustMethod = "BH")

eggosimp2 <- simplify(eggo2,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")

write.csv(eggosimp2@result,paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.barplot.pdf"),height=8,width=9)
barplot(eggosimp2,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.dotplot.pdf"),height=8,width=9)
dotplot(eggosimp2,showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.cnetplot.pdf"),height=9,width=12)
cnetplot(eggosimp2,circular=F,foldChange=FCgenelist, colorEdge =T, showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(eggosimp2)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.GO_enrichment.heatplot.pdf"),height=6,width=12)
heatplot(eggosimp2,foldChange=FCgenelist,showCategory=10)
dev.off()

### KEGG Enrichment
egkegg2 <- enrichKEGG(healthy_up$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
egkegg2 <- setReadable(egkegg2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(egkegg2@result,paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.KEGG_enrichment.csv"))

pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.KEGG_enrichment.barplot.pdf"),height=8,width=8)
barplot(egkegg2,showCategory=10,drop=T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.KEGG_enrichment.dotplot.pdf"),height=8,width=8)
dotplot(egkegg2,showCategory=10)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.KEGG_enrichment.cnetplot.pdf"),height=8,width=11)
cnetplot(egkegg2,circular=T, colorEdge =T)
dev.off()
pdf(paste0(dirRoot,"DE_Patient_Healthy_NK_cells/DE_Healthy_Up.NKcells.KEGG_enrichment.emapplot.pdf"),height=8,width=10)
emapplot(egkegg2)
dev.off()

