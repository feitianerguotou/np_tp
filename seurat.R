library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(MySeuratWrappers)


NP_1 <- Read10X("singlecellr/1_Cellranger_result/NP_2/")
NP_2 <- Read10X("singlecellr/1_Cellranger_result/NP_3/")
NP_3 <- Read10X("singlecellr/1_Cellranger_result/gse162122/061919/")
NP_4 <- Read10X("singlecellr/1_Cellranger_result/gse162122/11564/")
NP_5 <- Read10X("singlecellr/1_Cellranger_result/gse162122/11911/")
NP_6 <- Read10X("singlecellr/1_Cellranger_result/gse162122/12640/")
NP_7 <- Read10X("singlecellr/1_Cellranger_result/gse162122/12745/")

TNL_1 <- Read10X("singlecellr/1_Cellranger_result/TNL_1/")
TNL_2 <- Read10X("singlecellr/1_Cellranger_result/TNL_2/")
TNL_3 <- Read10X("singlecellr/1_Cellranger_result/TNL_3/")

NP_1object <- CreateSeuratObject(counts = NP_1, project = "NP_1",min.cells = 3, min.features = 200)
NP_2object <- CreateSeuratObject(counts = NP_2, project = "NP_2",min.cells = 3, min.features = 200)
NP_3object <- CreateSeuratObject(counts = NP_3, project = "NP_3",min.cells = 3, min.features = 200)
NP_4object <- CreateSeuratObject(counts = NP_4, project = "NP_4",min.cells = 3, min.features = 200)
NP_5object <- CreateSeuratObject(counts = NP_5, project = "NP_5",min.cells = 3, min.features = 200)
NP_6object <- CreateSeuratObject(counts = NP_6, project = "NP_6",min.cells = 3, min.features = 200)
NP_7object <- CreateSeuratObject(counts = NP_7, project = "NP_7",min.cells = 3, min.features = 200)

TNL_1object <- CreateSeuratObject(counts = TNL_1, project = "TNL_1",min.cells = 3, min.features = 200)
TNL_2object <- CreateSeuratObject(counts = TNL_2, project = "TNL_2",min.cells = 3, min.features = 200)
TNL_3object <- CreateSeuratObject(counts = TNL_3, project = "TNL_3",min.cells = 3, min.features = 200)

NP_1object[["percent.mt"]] <- PercentageFeatureSet(NP_1object, pattern = "^MT-")
NP_2object[["percent.mt"]] <- PercentageFeatureSet(NP_2object, pattern = "^MT-")
NP_3object[["percent.mt"]] <- PercentageFeatureSet(NP_3object, pattern = "^MT-")
NP_4object[["percent.mt"]] <- PercentageFeatureSet(NP_4object, pattern = "^MT-")
NP_5object[["percent.mt"]] <- PercentageFeatureSet(NP_5object, pattern = "^MT-")
NP_6object[["percent.mt"]] <- PercentageFeatureSet(NP_6object, pattern = "^MT-")
NP_7object[["percent.mt"]] <- PercentageFeatureSet(NP_7object, pattern = "^MT-")
TNL_1object[["percent.mt"]] <- PercentageFeatureSet(TNL_1object, pattern = "^MT-")
TNL_2object[["percent.mt"]] <- PercentageFeatureSet(TNL_2object, pattern = "^MT-")
TNL_3object[["percent.mt"]] <- PercentageFeatureSet(TNL_3object, pattern = "^MT-")

VlnPlot(NP_7object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(TNL_1object, feature1 = "nFeature_RNA", feature2 = "percent.mt")

NP_1object <- subset(NP_1object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_2object <- subset(NP_2object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_3object <- subset(NP_3object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_4object <- subset(NP_4object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_5object <- subset(NP_5object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_6object <- subset(NP_6object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
NP_7object <- subset(NP_7object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

TNL_1object <- subset(TNL_1object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_2object <- subset(TNL_2object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_3object <- subset(TNL_3object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

NP_TNL <- c(NP_1object,NP_2object,NP_3object,NP_4object,NP_5object,NP_6object,NP_7object,TNL_1object,TNL_2object,TNL_3object)
names(NP_TNL) <- c('NP_1','NP_2','NP_3','NP_4','NP_5','NP_6','NP_7','TNL_1','TNL_2','TNL_3')
NP_TNL <- lapply(NP_TNL, function(x) {x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)})

features <- SelectIntegrationFeatures(object.list = NP_TNL)
anchors <- FindIntegrationAnchors(object.list = NP_TNL, anchor.features = features)
NP_TNL <- IntegrateData(anchorset = anchors)

DefaultAssay(NP_TNL) <- "integrated"
NP_TNL <- ScaleData(NP_TNL, verbose = FALSE)
NP_TNL <- RunPCA(NP_TNL, verbose = FALSE)
NP_TNL <- RunUMAP(NP_TNL, reduction = "pca", dims = 1:10)
NP_TNL <- FindNeighbors(NP_TNL, reduction = "pca", dims = 1:10)

NP_TNL <- FindClusters(NP_TNL, resolution = 1)
NP_TNL@meta.data$groups[NP_TNL@meta.data$orig.ident=="TNL_3"]=c("TNL")
NP_TNL@meta.data$groups[NP_TNL@meta.data$groups=="TNL"]=c("TP")
saveRDS(NP_TNL, file = "singlecellr/til_tnl/NP_TNL2022.rds")#####
DimPlot(NP_TNL, reduction = "umap",split.by = "orig.ident")
Idents(NP_TNL) <- "orig.ident"
Idents(NP_TNL) <- "integrated_snn_res.1"
p1 <- DimPlot(NP_TNL, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
table(Idents(NP_TNL))
table(NP_TNL$groups)

DimPlot(NP_TNL, cells.highlight = CellsByIdentities(NP_TNL, idents = "35"))
plot <- DimPlot(NP_TNL, reduction = "umap")
NP_TNL <- CellSelector(plot=plot,object = NP_TNL,ident = "61")
plot <- DimPlot(NP_TNL, reduction = "umap")
NP_TNL <- CellSelector(plot=plot,object = NP_TNL,ident = "62")
DimPlot(NP_TNL, reduction = "umap")
NP_TNL@active.ident


genes_to_check3 = c('UBE2C', 'CDC20', 'TOP2A', 'STMN1','CD28', 'UCP2')#other
genes_to_check3 <- c('LUM','COL1A2','CD14','C1QB','TAGLN','ACTA2','TRAC','CD3E','TFF3','CCL21','VWF','PECAM1','KLRD1','KLRB1','TPSB2','TPSAB1','CD79A','CD79B','FCGR3B','NAMPT','HBB','HBA2','PAPPA2','PAPPA')
genes_to_check3 <- c('VWF', 'PECAM1','HLA-G','KRT7','DCN', 'LUM', 'TAGLN', 'ACTA2', 'HBB', 'HBA1','TRAC', 'CD3D','CD14', 'KLRD1','NCAM1', 'CD34')#all
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')#T
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('KRT7','HLA-G')
genes_to_check3 = c('LUM','CD14','TAGLN','TRAC','TFF3','VWF','KLRD1','TPSB2','CD79A','FCGR3B','HBB','PAPPA2')
VlnPlot(NP_TNL,features = genes_to_check3, pt.size = 0,ncol = 5)
MySeuratWrappers::VlnPlot(NP_TNL,features = genes_to_check3,stacked=T,direction = "vertical",pt.size=0)+NoLegend()
MySeuratWrappers::VlnPlot(NP_TNL,features = top10$gene[51:60],stacked=T,direction = "vertical",pt.size=0)+NoLegend()
MySeuratWrappers::MultiFeaturePlot(NP_TNL, features = mapgene, reduction = "umap",min.cutoff = 0,ncol =8)
FeaturePlot(NP_TNL, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "groups")
DefaultAssay(NP_TNL) <- "SCT"
NP_TNL_markers <- FindAllMarkers(NP_TNL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_markers %>% group_by(cluster) %>% top_n(n = 10)
write.csv(top10,file="singlecellr/np_top10.csv")


SMC=c(4,7,9,11,12,14,18)
Monocytic=c(1,6,20,21,25,27,31,34)
Epithelial_like_cells=c(19)
Endothelial_cells=c(2,3,15,17,32,62,63)
Fibroblasts=c(0,8)
LEC=c(13)
Red_blood_cells=c(22)
Neutrophils=c(30)
NK_cells=c(16)
T_cells=c(5,24)
B_cells=c(23)
Mast_cells=c(61)
Unknown=c(10,26,28,29,33,35)

current.cluster.ids <- c(SMC,Monocytic,Neutrophils,T_cells,B_cells,NK_cells,Mast_cells,
                         Epithelial_like_cells,
                         Endothelial_cells,
                         Fibroblasts,
                         LEC,Red_blood_cells,Unknown)

new.cluster.ids <- c(rep("SMC",length(SMC)),
                     rep("Monocytic",length(Monocytic)),
                     rep("Neutrophil",length(Neutrophils)),
                     rep("T_cells",length(T_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("Epithelial_like_cells",length(Epithelial_like_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("Fibroblasts",length(Fibroblasts)),
                     rep("LEC",length(LEC)),
                     rep("Red_blood_cells",length(Red_blood_cells)),
                     rep("Unknown",length(Unknown)))
NP_TNL@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(NP_TNL@active.ident)), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(NP_TNL, reduction = "umap",group.by = "celltype")

NP_TNL <- readRDS("singlecellr/np_tnl/np_tnl202210.rds")#####
saveRDS(NP_TNL, file = "singlecellr/np_tnl/np_tnl202210.rds")#####

NP_TNL_clear <- subset(NP_TNL, celltype == "Unknown", invert = T)
NP_TNL_clear <- readRDS("singlecellr/np_tnl202210/NP_TNL_clear.rds")#####
saveRDS(NP_TNL_clear, file = "singlecellr/np_tnl202210/NP_TNL_clear.rds")#####

DimPlot(NP_TNL_clear, reduction = "umap",group.by = "celltype",label = T,repel = T)
DefaultAssay(NP_TNL_clear) <- "SCT"#integrated SCT
Idents(NP_TNL_clear) <- "celltype"
NP_TNL_clear.markers <- FindAllMarkers(NP_TNL_clear, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_clear.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]

write.csv(top10,file="singlecellr/np_tnl202210/NP_TNL_clear.topmarkers.csv")
write.csv(NP_TNL_clear.markers,file="singlecellr/np_tnl202210/NP_TNL_clear.markers.csv")
DefaultAssay(NP_TNL_12) <- "integrated"

NP_TNL_clear_300 <- subset(NP_TNL_clear, downsample = 300)
DefaultAssay(NP_TNL_clear_300) <- "integrated" #integrated
DoHeatmap(NP_TNL_clear_300, features = top10_order$gene, angle = 25,size = 4,group.by = "celltype")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))

mapgene <- c('LUM','CD14','TAGLN','TRAC','TFF3','VWF','KLRD1','TPSB2','CD79A','FCGR3B','HBB','PAPPA2')
mapgene <- c('CD79A','VWF','KRT7','LUM','TFF3','TPSB2')
mapgene <- c('CD14','FCGR3B','KLRD1','HBB','ACTA2','CD3D')
mapgene <- c('CD79A','VWF','KRT7','LUM','TFF3','TPSB2','CD14','FCGR3B','KLRD1','HBB','ACTA2','CD3D')
mapgene <- c('ARNT2','RUNX1','MAFB','KLF5','POU2F2','NPAS2','MAFF','MAFK','PROX1','SPI1')
mapgene <- c('LUM','COL1A2','CD14','C1QB','TAGLN','ACTA2','CD2','CD3E','TFF3','CCL21','VWF','PECAM1','KLRD1','KLRB1','TPSB2','TPSAB1','CD79A','CD79B','FCGR3B','NAMPT','HBB','HBA2','KRT7','PAPPA')
MySeuratWrappers::VlnPlot(NP_TNL_clear,features = mapgene,stacked=T,direction = "vertical",pt.size=0, group.by = "celltype")+NoLegend()+theme(axis.text.y = element_text(size = 5))
FeaturePlot(NP_TNL_clear, features = mapgene, reduction = "umap",min.cutoff = 0)
DimPlot(NP_TNL_clear, reduction = "umap",group.by = "celltype",split.by = "groups")
MySeuratWrappers::MultiFeaturePlot(NP_TNL_clear, features = mapgene, reduction = "umap",min.cutoff = 0,ncol =8)

table(Idents(NP_TNL_clear), NP_TNL_clear$groups)
prop.table(table(Idents(NP_TNL_clear), NP_TNL_clear$groups))

#找差异基因
NP_TNL_clear$celltype.group <- paste(Idents(NP_TNL_clear),NP_TNL_clear$groups,sep = "_")
Idents(NP_TNL_clear) <- "celltype.group"
# Endothelial_cells Fibroblasts SMC Monocytic T_cells LEC NK_cells B_cells Epithelial_like_cells
DefaultAssay(NP_TNL_clear) <- "SCT"
FindMarkers(NP_TNL_clear, ident.1 = "Epithelial_like_cells_NP", ident.2 = "Epithelial_like_cells_TP", test.use = "wilcox",verbose = FALSE, min.pct = 0.5) %>% write.csv(file="singlecellr/np_tnl202210/Epithelial_like_cells_np_tp.csv")

Idents(NP_TNL_clear) <- "celltype"
mapgene <- c('ACTG2', 'GJA4', 'MYLK', 'CNN1','FTL', 'TIMP1','TAGLN', 'DES','ATP1A1', 'STMN1')#SMC
mapgene <- c('MYLK', 'OXTR','GJA1','PTGS2','PTGER3','ACTA2','PTGDS','GJA4')#shousuo
mapgene <- c('ENO1', 'COL3A1', 'COL5A2', 'IGFBP2','RPS27', 'IFI6','MDK', 'ADIRF','COL4A1', 'THBS1')#Fibroblasts
mapgene <- c('ALDOA', 'CRIP1', 'CCL4', 'SPP1','CCL3', 'ISG15','IFI6', 'HLA-A','CD14', 'CD163')#Monocytic
mapgene <- c('FN1', 'CCN2', 'C1QB', 'C1QA','IGFBP7')#Endothelial_cells
mapgene <- c('LYVE1', 'IFITM3', 'IGFBP6', 'GSN','AKAP12')#LEC
mapgene <- c('WFDC2', 'ESR1', 'SCGB2A1', 'MMP7','SOX17','BRI3','TIMP3', 'SERPINE2', 'PAPPA2','XAGE2')#Epithelial_like_cells
mapgene <- c('C1QB', 'C1QA', 'CCL3', 'TYROBP','RNASE1', 'CXCL14','IFITM1', 'FBLN1','SFRP1', 'PI16')#T_cells
mapgene <- c('C1QB', 'C1QA', 'IFITM3','PTPRCAP', 'IFITM1')#NK_cells
mapgene <- c('HLA-E', 'B2M', 'KLF2','DCN', 'CCDC80')#B_cells
mapgene <- c('LYVE1', 'IFITM3', 'JUND', 'IFITM1','IGFBP6', 'GSN','IFI6', 'FBLN1','DCN', 'LUM')#LEC
VlnPlot(NP_TNL_clear, features = mapgene, pt.size = 0,ncol = 10, split.by = "groups",idents = c("Monocytic"), y.lab = "", x.lab = "")
MySeuratWrappers::VlnPlot(NP_TNL_clear, features = mapgene, pt.size = 0,ncol =10, split.by = "groups",idents = c("Endothelial_cells"), y.lab = "", x.lab = "")   



#CellChat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(mindr)
library(NMF)
library(circlize)
library(Seurat)
library(SeuratData)

NP <- subset(NP_TNL_clear, groups == "NP")
NP <- SCTransform(NP)
DefaultAssay(NP) <- "SCT"
TP <- subset(NP_TNL_clear, groups == "TP")
TP <- SCTransform(TP)

cellchat_np <- createCellChat(object = NP, group.by = "celltype", assay = "RNA")
cellchat_tp <- createCellChat(object = TP, group.by = "celltype", assay = "RNA")
groupSize_np <- as.numeric(table(cellchat_np@idents))
groupSize_np <- as.numeric(table(cellchat_tp@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

CellChatDB_interaction <-CellChatDB$interaction
write.csv(CellChatDB_interaction, file = "singlecellr/np_tnl/CellChatDB_interaction.csv")
cellchat_np@DB <- CellChatDB.use
cellchat_tp@DB <- CellChatDB.use

cellchat_np <- subsetData(cellchat_np) 
cellchat_np <- identifyOverExpressedGenes(cellchat_np)
cellchat_np <- identifyOverExpressedInteractions(cellchat_np)
cellchat_np <- projectData(cellchat_np, PPI.human)  
cellchat_np <- computeCommunProb(cellchat_np)
cellchat_np <- filterCommunication(cellchat_np, min.cells = 10)
cellchat_np <- computeCommunProbPathway(cellchat_np)
cellchat_np <- aggregateNet(cellchat_np)

cellchat_np <- netAnalysis_computeCentrality(cellchat_np, slot.name = "netP")
cellchat_np@netP$pathways
cellchat_np@LR$LRsig$pathway_name
cellchat_np@LR$LRsig$antagonist

nPatterns = 3
cellchat_np <- identifyCommunicationPatterns(cellchat_np, pattern = "outgoing", k = nPatterns)
netVisual_heatmap(cellchat_np, signaling = c("CCL"), color.heatmap = "Reds")

all_inter_tn <- subsetCommunication(cellchat_np)
smc_tn <- subsetCommunication(cellchat_np, sources.use = c(8))

cellchat_tp <- subsetData(cellchat_tp) 
cellchat_tp <- identifyOverExpressedGenes(cellchat_tp)
cellchat_tp <- identifyOverExpressedInteractions(cellchat_tp)
cellchat_tp <- projectData(cellchat_tp, PPI.human)  
cellchat_tp <- computeCommunProb(cellchat_tp)
cellchat_tp <- filterCommunication(cellchat_tp, min.cells = 10)
cellchat_tp <- computeCommunProbPathway(cellchat_tp)
cellchat_tp <- aggregateNet(cellchat_tp)

cellchat_tp <- netAnalysis_computeCentrality(cellchat_tp, slot.name = "netP")
cellchat_tp@netP$pathways
cellchat_tp@LR$LRsig$pathway_name
cellchat_tp@LR$LRsig$antagonist

nPatterns = 3
cellchat_tp <- identifyCommunicationPatterns(cellchat_tp, pattern = "outgoing", k = nPatterns)
netVisual_heatmap(cellchat_tp, signaling = c("CCL"), color.heatmap = "Reds")

saveRDS(cellchat_tp, file = "singlecellr/np_tnl/cellchat_tp.rds")#####
saveRDS(cellchat_np, file = "singlecellr/np_tnl/cellchat_np.rds")#####
cellchat_tp <- readRDS("singlecellr/np_tnl/cellchat_tp.rds")
cellchat_np <- readRDS("singlecellr/np_tnl/cellchat_np.rds")

#merge
np_object.list <- list(NP = cellchat_np, TP = cellchat_tp)

cellchat_np_merge <- mergeCellChat(np_object.list, add.names = names(np_object.list))
cellchat_np_merge
gg1 <- compareInteractions(cellchat_np_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_np_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
compareInteractions(cellchat_np_merge, show.legend = F, group = c(1,2),size.text = 30)


netVisual_heatmap(cellchat_np)
netVisual_heatmap(cellchat_tp)
netVisual_heatmap(cellchat_np_merge, measure = "weight")
netVisual_heatmap(cellchat_np_merge)
netVisual_circle(cellchat_np@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_np@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat_tp@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_tp@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellchat_tp@net$weight)
netVisual_circle(cellchat_np@net$count)
#选5个 #红色]（或[蓝色]边表示信号在第二个数据集中增加或[减少]）
cellchat_tp_11<- subsetCellChat(object = cellchat_tp, idents.use = c("Red_blood_cells"),invert = T)
cellchat_np_11<- subsetCellChat(object = cellchat_np, idents.use = c("Red_blood_cells"),invert = T)
netVisual_heatmap(cellchat_tp_11)
netVisual_heatmap(cellchat_np_11)
np_object.list <- list(NP = cellchat_np_11, TP = cellchat_tp_11)

cellchat_np_merge <- mergeCellChat(np_object.list, add.names = names(np_object.list))
cellchat_np_merge
gg1 <- compareInteractions(cellchat_np_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_np_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

netVisual_diffInteraction(cellchat_np_merge, weight.scale = T,label.edge = T)
netVisual_diffInteraction(cellchat_np_merge, weight.scale = T, measure = "weight",label.edge = T)
netVisual_heatmap(cellchat_np_merge)
saveRDS(cellchat_np_merge, file = "singlecellr/np_tnl/cellchat_np_11_merge.rds")
cellchat_np_5_merge <- readRDS("singlecellr/np_tnl/cellchat_np_5_merge.rds")

rankSimilarity(cellchat_np_merge, type = "functional")
gg1 <- rankNet(cellchat_np_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_np_merge, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
cellchat_np_11 <- netAnalysis_computeCentrality(cellchat_np_11, slot.name = "netP")
cellchat_tp_11 <- netAnalysis_computeCentrality(cellchat_tp_11, slot.name = "netP")

pathway.union <- union(cellchat_np_11@netP$pathways, cellchat_tp_11@netP$pathways)
netAnalysis_signalingRole_heatmap(cellchat_np_11, width = 5, height = 16,
                                  pattern = "all")
netAnalysis_signalingRole_heatmap(cellchat_tp_11, width = 5, height = 16,
                                  pattern = "outgoing", signaling = pathway.union)



library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets 
netAnalysis_signalingRole_heatmap(cellchat_np_11, width = 5, height = 12,
                                  pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat_tp_4, width = 5, height = 18,
                                  pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat_np_4, width = 5, height = 12,
                                  pattern = "outgoing",color.heatmap = "GnBu")
netAnalysis_signalingRole_heatmap(cellchat_tp_4, width = 5, height = 18,
                                  pattern = "outgoing",color.heatmap = "GnBu")
netVisual_bubble(cellchat_np_merge, comparison = c(1, 2), angle.x = 45, thresh = 0.01)
netVisual_bubble(cellchat_np_merge, comparison = c(1, 2), max.dataset = 2, angle.x = 45, remove.isolate = T)

cellchat_np_merge <- identifyOverExpressedGenes(cellchat_np_merge, 
                                                pos.dataset = "TP", features.name = "TP", 
                                                only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
net <- netMappingDEG(cellchat_np_merge, features.name = "TP")
net.up <- subsetCommunication(cellchat_np_merge, net = net, datasets = "TP",ligand.logFC = 1, receptor.logFC = 1)
net.down <- subsetCommunication(cellchat_np_merge, net = net, datasets = "NP",ligand.logFC = -1, receptor.logFC = -1)
write.csv(net,"singlecellr/np_tnl/cellchat_np_merge202210des.csv")
netVisual_chord_gene(cellchat_tp_11, targets.use = "SMC",
                     slot.name = 'net', net = net.up ,lab.cex = 0.8, small.gap = 3.5)
netVisual_chord_gene(cellchat_tp_11, sources.use = "SMC",
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5)
netVisual_chord_gene(cellchat_np_11, targets.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5 )
netVisual_chord_gene(cellchat_np_11, sources.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5 )


netVisual_bubble(cellchat_np_merge, pairLR.use = use.up, targets.use = "SMC", comparison = c(1, 2),angle.x = 90, remove.isolate = T)

pairLR.use.down = net.down[, "interaction_name", drop = F]




#SMC
NP_TNL_SMC <- subset(NP_TNL, celltype == "SMC")
NP_TNL_SMC <- ScaleData(NP_TNL_SMC, verbose = FALSE)
NP_TNL_SMC <- RunPCA(NP_TNL_SMC, npcs = 30, verbose = FALSE)
NP_TNL_SMC <- RunUMAP(NP_TNL_SMC, reduction = "pca", dims = 1:30)
NP_TNL_SMC <- FindNeighbors(NP_TNL_SMC, reduction = "pca", dims = 1:30)

p1 <- DimPlot(NP_TNL_SMC, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_SMC, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DefaultAssay(NP_TNL_SMC) <- "integrated"
NP_TNL_SMC <- FindClusters(NP_TNL_SMC, resolution = 0.1)#分群数
DimPlot(NP_TNL_SMC, reduction = "umap", group.by = "groups",)
library(clustree)
clustree(NP_TNL_SMC@meta.data, prefix = "integrated_snn_res.")

NP_TNL_SMC <- FindClusters(NP_TNL_SMC, resolution = 0.5)

genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'ACTA1','GJA1')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA','OXTR', 'GJA1')#代谢
genes_to_check3 <- c('TPM2', 'TES','VIM','ZEB2','UTRN', 'MTATP6P1', 'TMSB4X', 'XIST', 'VCAN', 'TFRC','TYROBP')#all
genes_to_check3 = c('UBE2C', 'CDC20', 'TOP2A', 'STMN1','CD28', 'UCP2')#other
genes_to_check3 = c('CCN1', 'ACTG2', 'ADIRF', 'SNCG','STEAP4', 'LHFPL6','C7', 'LUM')#marker
genes_to_check3 <- c('OXTR', 'GJA1','ELANE', 'IFNG','SERPINE1')
VlnPlot(NP_TNL_SMC,features = genes_to_check3, pt.size = 0,ncol = 2)
FeaturePlot(NP_TNL_SMC, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
DoHeatmap(object = NP_TNL_SMC, features = top10$gene)
DoHeatmap(NP_TNL_SMC, features = top10$gene, group.by = "SMC_type", angle = 25)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#标注
SMC_1=c(0,4)
SMC_2=c(1)
SMC_3=c(2)
SMC_4=c(3)

current.cluster.ids <- c(SMC_1,
                         SMC_2,
                         SMC_3,
                         SMC_4)

new.cluster.ids <- c(rep("SMC_1",length(SMC_1)),
                     rep("SMC_2",length(SMC_2)),
                     rep("SMC_3",length(SMC_3)),
                     rep("SMC_4",length(SMC_4)))

NP_TNL_SMC@meta.data$SMC_type <- plyr::mapvalues(x = as.integer(as.character(NP_TNL_SMC@meta.data$integrated_snn_res.0.1)), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(NP_TNL_SMC, reduction = "umap",group.by = "SMC_type", split.by = "groups")
DimPlot(NP_TNL_SMC, reduction = "umap",group.by = "SMC_type")

table(Idents(NP_TNL_SMC))
table(NP_TNL_SMC$SMC_type,NP_TNL_SMC$groups)

Idents(NP_TNL_SMC) <- "SMC_type"
DefaultAssay(NP_TNL_SMC) <- "SCT"
DefaultAssay(NP_TNL_SMC) <- "integrated"
NP_TNL_SMC_markers <- FindAllMarkers(NP_TNL_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_SMC_markers %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(NP_TNL_SMC, features = top10_order$gene,slot = "scale.data", angle = 20,label = F,group.by = "SMC_type")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

NP_TNL_SMC_300 <- subset(NP_TNL_SMC, downsample = 300)
DefaultAssay(NP_TNL_SMC_300) <- "SCT" #integrated
DefaultAssay(NP_TNL_SMC_300) <- "integrated" #integrated
DoHeatmap(NP_TNL_SMC_300, features = top10_order$gene, angle = 25,size = 4,group.by = "SMC_type")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))



write.csv(NP_TNL_SMC_markers,file="singlecellr/NP_TNL_SMC_markers202210.csv")
VlnPlot(NP_TNL_SMC,features = top10$gene, pt.size = 0,ncol = 2)
SMC_marker <- c('STEAP4', 'THY1','GGT5','DCN','CLDN1', 'BNC2', 'DES', 'ACTG2', 'TES', 'RERGL', 'NTRK2', 'SNCG')
SMC_marker <- c('OXTR','GJA1','ELANE','IFNG','PTGFR','PTGS2')
VlnPlot(NP_TNL_SMC,features = SMC_marker, pt.size = 0,ncol = 3,group.by = "SMC_type")
MySeuratWrappers::VlnPlot(NP_TNL_SMC,features = SMC_marker,stacked=F,pt.size=0,group.by = "SMC_type",x.lab = "",y.lab = "",ncol = 3)
MySeuratWrappers::MultiFeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,ncol =8)
VlnPlot(NP_TNL_SMC,features = SMC_marker, pt.size = 0,ncol = 3,group.by = "groups")
FeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,keep.scale = "all",split.by = "groups",ncol = 4)
FeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,keep.scale = "all",ncol = 4)
saveRDS(NP_TNL_SMC, file = "singlecellr/np_tnl/NP_TNL_SMC20221025.rds")#####
NP_TNL_SMC <- readRDS("singlecellr/np_tnl/NP_TNL_SMC20221025.rds")#####

#smc空间
NP_v <- readRDS("singlecellr/np_tnl/NP_v.rds")#####
TP_v <- readRDS("singlecellr/np_tnl/TP_v.rds")#####

NP_TNL_SMC_NP <- subset(NP_TNL_SMC, groups == "NP")
NP_TNL_SMC_TP <- subset(NP_TNL_SMC, groups == "TNL")

NP_TNL_SMC_NP <- SCTransform(NP_TNL_SMC_NP)
DefaultAssay(NP_TNL_SMC_NP) <- "SCT"
NP_TNL_SMC_TP <- SCTransform(NP_TNL_SMC_TP)
DefaultAssay(NP_TNL_SMC_TP) <- "SCT"
DefaultAssay(NP_v) <- "SCT"
DefaultAssay(TP_v) <- "SCT"
anchors_SMC_NP <- FindTransferAnchors(reference = NP_TNL_SMC_NP, query = NP_v, normalization.method = "SCT")
predictions.assay_SMC_NP <- TransferData(anchorset = anchors_SMC_NP, refdata = NP_TNL_SMC_NP$SMC_type, prediction.assay = TRUE, 
                                         weight.reduction = NP_v[["pca"]], dims = 1:30)
NP_v_SMC <- NP_v
NP_v_SMC[["predictions"]] <- predictions.assay_SMC_NP
DefaultAssay(NP_v_SMC) <- "predictions"
SpatialFeaturePlot(NP_v_SMC, features = c("SMC-1","SMC-2","SMC-3","SMC-4"), alpha = c(0.1, 1), min.cutoff = 0)

anchors_SMC_TP <- FindTransferAnchors(reference = NP_TNL_SMC_TP, query = TP_v, normalization.method = "SCT")
predictions.assay_SMC_TP <- TransferData(anchorset = anchors_SMC_TP, refdata = NP_TNL_SMC_TP$SMC_type, prediction.assay = TRUE, 
                                         weight.reduction = TP_v[["pca"]], dims = 1:30)
TP_v_SMC <- TP_v
TP_v_SMC[["predictions"]] <- predictions.assay_SMC_TP
DefaultAssay(TP_v_SMC) <- "predictions"
SpatialFeaturePlot(TP_v_SMC, features = c("SMC-1","SMC-2","SMC-3","SMC-4"), alpha = c(0.1, 1), min.cutoff = 0)

DefaultAssay(NP_TNL_SMC) <- "integrated"
DoHeatmap(object = NP_TNL_SMC, features = top10$gene)

#smc空间
NP_v <- readRDS("singlecellr/np_tnl/NP_v.rds")#####
TP_v <- readRDS("singlecellr/np_tnl/TP_v.rds")#####

NP_TNL_NP <- subset(NP_TNL_clear, groups == "NP")
NP_TNL_TP <- subset(NP_TNL_clear, groups == "TP")

NP_TNL_NP <- SCTransform(NP_TNL_NP)
DefaultAssay(NP_TNL_NP) <- "SCT"
NP_TNL_TP <- SCTransform(NP_TNL_TP)
DefaultAssay(NP_TNL_TP) <- "SCT"
DefaultAssay(NP_v) <- "SCT"
DefaultAssay(TP_v) <- "SCT"
anchors_NP <- FindTransferAnchors(reference = NP_TNL_NP, query = NP_v, normalization.method = "SCT")
predictions.assay_NP <- TransferData(anchorset = anchors_NP, refdata = NP_TNL_NP$celltype, prediction.assay = TRUE, 
                                         weight.reduction = NP_v[["pca"]], dims = 1:30)

seleccells <- c("SMC","LEC","Endothelial-cells","Fibroblasts","Monocytic","T-cells","B-cells","Mast-cells","Neutrophil","NK-cells","Epithelial-like-cells","Red-blood-cells")

NP_v[["predictions"]] <- predictions.assay_NP
DefaultAssay(NP_v) <- "predictions"
SpatialFeaturePlot(NP_v, features = seleccells, alpha = c(0.1, 1), min.cutoff = 0)

anchors_TP <- FindTransferAnchors(reference = NP_TNL_TP, query = TP_v, normalization.method = "SCT")
predictions.assay_TP <- TransferData(anchorset = anchors_TP, refdata = NP_TNL_TP$celltype, prediction.assay = TRUE, 
                                         weight.reduction = TP_v[["pca"]], dims = 1:30)

TP_v[["predictions"]] <- predictions.assay_TP
DefaultAssay(TP_v) <- "predictions"
SpatialFeaturePlot(TP_v, features = seleccells, alpha = c(0.1, 1), min.cutoff = 0)

DefaultAssay(TP_v) <- "SCT"
DefaultAssay(NP_v) <- "SCT"

genes_to_check3 <- c('LUM','CD14','TAGLN','CD3E','TFF3','VWF','PECAM1','KLRD1','TPSB2','CD79A','FCGR3B','PAPPA')
genes_to_check3 <- c('CD79A','VWF','KRT7','LUM','TFF3','TPSB2','CD14','FCGR3B','KLRD1','HBB','ACTA2','CD3D')
SpatialFeaturePlot(TP_v, features = genes_to_check3, alpha = c(0.1, 1), min.cutoff = 0)
SpatialFeaturePlot(NP_v, features = genes_to_check3, alpha = c(0.1, 1), min.cutoff = 0)
FeaturePlot(NP_v, features = genes_to_check3, reduction = "umap",min.cutoff = 0.1,keep.scale ="all")
FeaturePlot(TP_v, features = genes_to_check3, reduction = "umap",min.cutoff = 0.1,keep.scale ="all")

#Fibroblasts
NP_TNL_Fibroblasts <- subset(NP_TNL, celltype == "Fibroblasts")
NP_TNL_Fibroblasts <- ScaleData(NP_TNL_Fibroblasts, verbose = FALSE)
NP_TNL_Fibroblasts <- RunPCA(NP_TNL_Fibroblasts, npcs = 30, verbose = FALSE)
NP_TNL_Fibroblasts <- RunUMAP(NP_TNL_Fibroblasts, reduction = "pca", dims = 1:10)
NP_TNL_Fibroblasts <- FindNeighbors(NP_TNL_Fibroblasts, reduction = "pca", dims = 1:10)
Idents(NP_TNL_Fibroblasts) <- "groups"
Idents(NP_TNL_Fibroblasts) <- "integrated_snn_res.0.1"
Idents(NP_TNL_Fibroblasts) <- "fibor"
DefaultAssay(NP_TNL_Fibroblasts) <- "integrated"
p1 <- DimPlot(NP_TNL_Fibroblasts, reduction = "umap", group.by = "integrated_snn_res.0.5")
p2 <- DimPlot(NP_TNL_Fibroblasts, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
NP_TNL_Fibroblasts@meta.data$fibor 
NP_TNL_Fibroblasts@meta.data$fibor[NP_TNL_Fibroblasts@meta.data$integrated_snn_res.0.1=="2"]=c("MyoFibroblasts")
NP_TNL_Fibroblasts <- FindClusters(NP_TNL_Fibroblasts, resolution = 0.1)
DimPlot(NP_TNL_Fibroblasts, reduction = "umap", split.by = "groups",group.by = "fibor")
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'MYH11','GJA1')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('APOD', 'MMP2', 'SPOCK1', 'PRG4','TPPP3', 'CD55','ACTA2','MYH11', 'ACTG2')
genes_to_check3 = c('SERPINE2', 'IGFBP6', 'RGCC', 'NR4A1','EGR1', 'NDRG2','MYH11','ACTA2', 'MYLK')
DefaultAssay(NP_TNL_Fibroblasts) <- "SCT"
VlnPlot(NP_TNL_Fibroblasts,features = genes_to_check3, pt.size = 0,ncol = 3)
FeaturePlot(NP_TNL_Fibroblasts, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all",)
FeaturePlot(NP_TNL_Fibroblasts, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all",split.by = "groups")

MySeuratWrappers::VlnPlot(NP_TNL_Fibroblasts, features = genes_to_check3, pt.size = 0,group.by = "fibor",ncol =3, y.lab = "", x.lab = "",stacked = F)

NP_TNL_Fibroblasts_markers <- FindAllMarkers(NP_TNL_Fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_Fibroblasts_markers %>% group_by(cluster) %>% top_n(n = 10, wt=avg_log2FC)
DoHeatmap(object = NP_TNL_Fibroblasts, features = top10$gene,angle = 30)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

NP_TNL_Fibroblasts_300 <- subset(NP_TNL_Fibroblasts, downsample = 300)
DefaultAssay(NP_TNL_Fibroblasts_300) <- "integrated" #integrated
DoHeatmap(NP_TNL_Fibroblasts_300, features = top10$gene, angle = 25,size = 4,group.by = "fibor")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))
saveRDS(NP_TNL_Fibroblasts, file = "singlecellr/np_tnl/NP_TNL_Fibroblasts20221025.rds")#####
NP_TNL_Fibroblasts <- readRDS("singlecellr/np_tnl/NP_TNL_Fibroblasts20221025.rds")#####

#单核细胞
NP_TNL_monocytic <- subset(NP_TNL, celltype == "Monocytic")
NP_TNL_monocytic <- ScaleData(NP_TNL_monocytic, verbose = FALSE)
NP_TNL_monocytic <- RunPCA(NP_TNL_monocytic, npcs = 30, verbose = FALSE)
NP_TNL_monocytic <- RunUMAP(NP_TNL_monocytic, reduction = "pca", dims = 1:10)
NP_TNL_monocytic <- FindNeighbors(NP_TNL_monocytic, reduction = "pca", dims = 1:10)
Idents(NP_TNL_monocytic) <- "groups"
DefaultAssay(NP_TNL_monocytic) <- "integrated"
DefaultAssay(NP_TNL_monocytic) <- "SCT"
p1 <- DimPlot(NP_TNL_monocytic, reduction = "umap", group.by = "integrated_snn_res.0.1")
p2 <- DimPlot(NP_TNL_monocytic, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

NP_TNL_monocytic <- FindClusters(NP_TNL_monocytic, resolution = 0.3)
DimPlot(NP_TNL_monocytic, reduction = "umap", group.by = "integrated_snn_res.0.3")
genes_to_check3 = c('CD163', 'CD86', 'CD68', 'CD14', 'CD209','IL1B', 'MRC1','S100A9','S100A8','CIITA','ITGAX','FCGR3A')#巨噬
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('CD14','CD163', 'CD86','MRC1', 'FCN1','IL1B')
genes_to_check3 = c('CD163', 'CD86', 'CD14', 'PTPRC','IL1B', 'MRC1','S100A9','S100A8','LYZ','TYROBP')#巨噬
genes_to_check3 = c('FCGR3A', 'NLRP3', 'CD14', 'FCN1','IL1B', 'C1QC','SPP1')#巨噬
genes_to_check3 = c('CD40', 'CD80', 'CD86','LILRA4', 'CD1C', 'BATF3')#DC
genes_to_check3 = c('CLEC9A', 'XCR1', 'CADM1','CD1C', 'FCER1A', 'CLEC10A','LAMP3', 'CD80', 'CCR7')#DC2
genes_to_check3 = c('C1QB', 'SPP1','IL1B','FCN1','COL1A1','CCL21','TAGLN','IGFBP5','CD1C','FCER1A','STMN1','MKI67','HBA1','HBA2')
genes_to_check3 = c('IL1B', 'FCN1','CD163','MRC1','GZMB','GZMA','TAGLN','IGFBP5','AQP1','CLDN5','CCL21','AKAP12')
genes_to_check3 = c('IL1B', 'CD163','GZMA','TAGLN','AQP1','CCL21')
DefaultAssay(NP_TNL_monocytic) <- "SCT"
VlnPlot(NP_TNL_monocytic,features = genes_to_check3, pt.size = 0,ncol = 4)
FeaturePlot(NP_TNL_monocytic, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all")
FeaturePlot(NP_TNL_monocytic, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all",split.by = "groups")
saveRDS(NP_TNL_monocytic, file = "singlecellr/np_tnl/NP_TNL_monocytic.rds")

NP_TNL_monocytic_markers <- FindAllMarkers(NP_TNL_monocytic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_monocytic_markers %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = NP_TNL_monocytic, features = top10$gene)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#标注
write.csv(NP_TNL_monocytic_markers,file="singlecellr/NP_TNL_monocytic_markers20221028.csv")

genes_to_check3 = c('CD163', 'CD86', 'CD14', 'CD209','IL1B', 'MRC1','S100A9','S100A8','CD68','MMP19')#巨噬
genes_to_check3 = c('C1QA', 'MRC1','S100A8','FCN1','CALD1','RGS5','MMRN1','CAVIN2','MKI67','TOP2A')#
VlnPlot(NP_TNL_monocytic,features = genes_to_check3, pt.size = 0,ncol = 2, combine = TRUE)
FeaturePlot(NP_TNL_monocytic, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
FeaturePlot(NP_TNL_monocytic, features = genes_to_check3, reduction = "umap",min.cutoff = 0,pt.size = 0.5,split.by = "groups")
MySeuratWrappers::VlnPlot(NP_TNL_monocytic, features = genes_to_check3, pt.size = 0,ncol =2, y.lab = "", x.lab = "",stacked = F,group.by = "monocytic")

#标注
Idents(NP_TNL_monocytic) <- "integrated_snn_res.0.1"
M1=c(3)
M2=c(0,1,2,4)
Mono_TAGLN=c(5)
Mono_CCL21=c(8)
Mono_AQP1=c(6)
Mono_GZMA=c(7)

current.cluster.ids <- c(M1,M2,Mono_TAGLN,Mono_CCL21,Mono_AQP1,Mono_GZMA)

new.cluster.ids <- c(rep("M1",length(M1)),
                     rep("M2",length(M2)),
                     rep("Mono_TAGLN",length(TAGLN_Monocytes)),
                     rep("Mono_CCL21",length(CCL21_Monocytes)),
                     rep("Mono_AQP1",length(AQP1_Monocytes)),
                     rep("Mono_GZMA",length(GZMA_Monocytes)))

NP_TNL_monocytic@meta.data$monocytic <- plyr::mapvalues(x = as.integer(as.character(NP_TNL_monocytic@meta.data$integrated_snn_res.0.3)), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(NP_TNL_monocytic, reduction = "umap",group.by = "monocytic", split.by = "groups")
DimPlot(NP_TNL_monocytic, reduction = "umap")

NP_TNL_monocytic_300 <- subset(NP_TNL_monocytic, downsample = 100)
DefaultAssay(NP_TNL_monocytic_300) <- "integrated" #integrated
DoHeatmap(NP_TNL_monocytic_300, features = top10_order$gene, angle = 25,size = 4,group.by = "monocytic")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))


Idents(NP_TNL_monocytic) <- "monocytic"
saveRDS(NP_TNL_monocytic, file = "singlecellr/np_tnl/NP_TNL_monocytic20221028.rds")
NP_TNL_monocytic <- readRDS("singlecellr/np_tnl/NP_TNL_monocytic20221028.rds")#####

#NP_TNL_monocytic空间
NP_TNL_monocytic_NP <- subset(NP_TNL_monocytic, groups == "NP")
NP_TNL_monocytic_TP <- subset(NP_TNL_monocytic, groups == "TP")

NP_TNL_monocytic_NP <- SCTransform(NP_TNL_monocytic_NP)
DefaultAssay(NP_TNL_monocytic_NP) <- "SCT"
NP_TNL_monocytic_TP <- SCTransform(NP_TNL_monocytic_TP)
DefaultAssay(NP_TNL_monocytic_TP) <- "SCT"
DefaultAssay(NP_v) <- "SCT"
DefaultAssay(TP_v) <- "SCT"



#配体受体 
NP_v@reductions$spatial = NP_v@reductions$umap
NP_v@reductions$spatial@key = 'spatial_'
NP_v@reductions$spatial@cell.embeddings = as.matrix(NP_v@images$slice1@coordinates[,c(3,2)])
NP_v@reductions$spatial@cell.embeddings[,1] = -NP_v@reductions$spatial@cell.embeddings[,1]
colnames(NP_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(NP_v, features = c("CD3E","CD3D"), blend = T,cols = c('lightgrey','blue3','firebrick3'),reduction = 'spatial',combine = T)+coord_flip()+NoLegend()&NoAxes()

TP_v@reductions$spatial = TP_v@reductions$umap
TP_v@reductions$spatial@key = 'spatial_'
TP_v@reductions$spatial@cell.embeddings = as.matrix(TP_v@images$slice1@coordinates[,c(3,2)])
TP_v@reductions$spatial@cell.embeddings[,1] = -TP_v@reductions$spatial@cell.embeddings[,1]
colnames(TP_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(TP_v, features = c("CD3E","CD3D"), blend = T,cols = c('lightgrey','blue3','firebrick3'),reduction = 'spatial',combine = T)+coord_flip()+NoLegend()&NoAxes()



#Tcells
DefaultAssay(NP_TNL) <- "integrated"
Idents(NP_TNL) <- NP_TNL$celltype
NP_TNL_TNK <- subset(NP_TNL, idents = c("T_cells"))
NP_TNL_TNK <- ScaleData(NP_TNL_TNK, verbose = FALSE)
NP_TNL_TNK <- RunPCA(NP_TNL_TNK, verbose = FALSE)
ElbowPlot(NP_TNL_TNK)
NP_TNL_TNK <- RunUMAP(NP_TNL_TNK, reduction = "pca", dims = 1:5)
NP_TNL_TNK <- FindNeighbors(NP_TNL_TNK, reduction = "pca", dims = 1:5)
Idents(NP_TNL_TNK) <- "celltype"
DimPlot(NP_TNL_TNK, reduction = "umap",group.by = "groups")

p1 <- DimPlot(NP_TNL_TNK, reduction = "umap", group.by = "groups",split.by = "groups")
p2 <- DimPlot(NP_TNL_TNK, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(NP_TNL_TNK, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")
DefaultAssay(NP_TNL_TNK) <- "integrated"
NP_TNL_TNK <- FindClusters(object = NP_TNL_TNK, resolution = 0.1)

TNL_TIL_TNK@meta.data$Ttype[TNL_TIL_TNK@meta.data$integrated_snn_res.0.2=="8"]=c("CD8+ T_cells")
TNL_TIL_TNK$Ttype.group.new <- paste(Idents(TNL_TIL_TNK),TNL_TIL_TNK$groups2,sep = "_")
table(Idents(TNL_TIL_TNK))
table(TNL_TIL_TNK$Macrotype,TNL_TIL_macro$groups)
Idents(NP_TNL_TNK) <- "Ttype"
Idents(NP_TNL_TNK) <- "Ttype.group.new"
Idents(NP_TNL_TNK) <- "integrated_snn_res.0.2"
DefaultAssay(NP_TNL_TNK) <- "SCT"
genes_to_check3 = c('CD3E', 'CD3D', 'CD2')#T
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5')#免疫抑制
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('CD4','CD8A', 'BTLA','CD69', 'CCR7','CCL5','KLRD1','C1QB', 'GSN','CXCL14', 'GNLY','FCGR3A', 'CD2','TOP2A')
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CCL4','TNFSF10','IL7R','CCL3','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL4','TNFSF10','IL7R','CCL3','CCL5','XCL1','XCL2','IL2RB')#T NK炎性
genes_to_check3 = c('CD1C', 'LUM','IL1B','MRC1','MKI67','HBA1','TAGLN')#top
genes_to_check3 = c('CCR7', 'IL7R','FOXP3','CXCL13','ZFP36','GZMK','IFIT1','IFNG','LAG3','AREG','TRDC','CD7')#T marker
genes_to_check3 = c('CCR7','IL7R','CD8A','TAGLN','KLRD1','CD4','MKI67','KIT')#top
genes_to_check3 = c('CD4','CD8A','ACTA2','KLRD1','MAFB','MKI67','KIT','FOXP3')#top
genes_to_check3 = c('CD38','IFNG','ICOS','CCR6','ID2','CD2')#top
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67')
genes_to_check3 = c('C1QB', 'SPP1','IL1B','FCN1','COL1A1','CCL21','TAGLN','IGFBP5','CD1C','FCER1A','STMN1','MKI67','HBA1','HBA2')
genes_to_check3 = c('C1QB', 'MRC1','IL1B','FCN1','COL1A1','LUM','TAGLN','RGS5','CD1C','FCER1A','MKI67','HBA1')
FeaturePlot(NP_TNL_TNK, features = genes_to_check3, reduction = "umap",min.cutoff = 0,keep.scale = "all",split.by = "groups")
MySeuratWrappers::MultiFeaturePlot(NP_TNL_TNK,features = genes_to_check3,ncol=4)
saveRDS(NP_TNL_TNK, file = "singlecellr/np_tnl/NP_TNL_TNK20221028.rds")#####
NP_TNL_TNK <-readRDS(file = "singlecellr/np_tnl/NP_TNL_TNK20221028.rds")#####

VlnPlot(NP_TNL_TNK, features = genes_to_check3, pt.size = 0,ncol = 4)
VlnPlot(NP_TNL_TNK, features = genes_to_check3, pt.size = 0,ncol = 4)
MySeuratWrappers::VlnPlot(NP_TNL_TNK,features = genes_to_check3,stacked=T,pt.size=0,
                          cols = color7)
color7 <- c("#F8766D","#F8766D", "#CD9600","#CD9600", "#7CAE00","#7CAE00", "#00BE67","#00BE67", "#00BFC4","#00BFC4", "#00A9FF","#00A9FF", "#C77CFF", "#C77CFF", "#FF61CC", "#FF61CC")
NP_TNL_TNK.markers <- FindAllMarkers(NP_TNL_TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_TNK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
set.seed(42)
TNL_TIL_TNK_80 <- subset(TNL_TIL_TNK, downsample = 80)
DefaultAssay(TNL_TIL_TNK_80) <- "integrated"
DoHeatmap(TNL_TIL_TNK_80, features = top10_order$gene, angle = 15,group.by = "Ttype")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

DoHeatmap(object = TNL_TIL_TNK, features = genes_to_check3)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13_TNK.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_13_macro.csv")
write.csv(TNL_TIL_TNK.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_new_TNK.csv")

#(Endothelial_cells)
NP_TNL_endo <- subset(NP_TNL, celltype == c("Endothelial_cells"))
NP_TNL_endo <- ScaleData(NP_TNL_endo, verbose = FALSE)
NP_TNL_endo <- RunPCA(NP_TNL_endo, verbose = FALSE)
ElbowPlot(NP_TNL_endo)
NP_TNL_endo <- RunUMAP(NP_TNL_endo, reduction = "pca", dims = 1:30)
NP_TNL_endo <- FindNeighbors(NP_TNL_endo, reduction = "pca", dims = 1:30)
Idents(NP_TNL_endo) <- "celltype"
DimPlot(NP_TNL_endo, reduction = "umap",split.by = "groups")

p1 <- DimPlot(NP_TNL_endo, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")
DefaultAssay(NP_TNL_endo) <- "integrated"
NP_TNL_endo <- FindClusters(object = NP_TNL_endo, resolution = 0.1)

NP_TNL_endo@meta.data$endo_type[NP_TNL_endo@meta.data$integrated_snn_res.0.1=="2"]=c("Arterial")
NP_TNL_endo$endo_type.group.new <- paste(Idents(NP_TNL_endo),NP_TNL_endo$groups,sep = "_")
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "endo_type.group.new")
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "endo_type.group.new")
table(Idents(TNL_TIL_B))
table(TNL_TIL_endo$endo_type,TNL_TIL_endo$groups)
Idents(NP_TNL_endo) <- "integrated_snn_res.0.2"
Idents(NP_TNL_endo) <- "endo_type"
DefaultAssay(TNL_TIL_endo) <- "SCT"
genes_to_check3 = c('PROX1','CCL21', 'PROX1', 'HEY1','CXCL12', 'IGFBP3','CD36', 'CA4','ACKR1')#ymphatic ECs (LECs; CCL21, PROX1).arteries (HEY1, IGFBP3), capillaries (CD36, CA4), veins (ACKR1)
genes_to_check3 = c('ACKR1', 'IGFBP3')#ymphatic 
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('APLNR', 'INSR', 'ESM1', 'KDR','VWA1', 'COL4A1')#angiogenic
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('VWF','THBD','EDN1', 'SELP','OLR1')#血管内皮损伤marker
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('MGP','BMX', 'PRSS23','SLC6A2')
genes_to_check3 = c('TNFSF10')#endo炎性因子2
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(NP_TNL_endo, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_endo,features = genes_to_check3,ncol=2)
VlnPlot(NP_TNL_endo, features = genes_to_check3, pt.size = 0, group.by = "endo_type")
saveRDS(NP_TNL_endo, file = "singlecellr/np_tnl/np_tnl_endo20221113.rds")#####
TNL_TIL_endo <-readRDS(file = "singlecellr/np_tnl/np_tnl_endo.rds")#####

NP_TNL_endo.markers <- FindAllMarkers(NP_TNL_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_endo, features = top10$gene,group.by = "NP_TNL_endo")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#Venous Arterial
Idents(NP_TNL_endo) <- "endo_type.group.new"
FindMarkers(NP_TNL_endo, ident.1 = "Arterial_1TNL", ident.2 = "Arterial_2TIL",test.use = "wilcox", verbose = FALSE, min.pct = 0.5) %>% write.csv(file="singlecellr/til_tnl/endo_cells_TNL_TIL_Arterial.csv")

#(Epithelial_like_cells)
NP_TNL_epi <- subset(NP_TNL, celltype == c("Trophectoderm_cells"))
NP_TNL_epi <- ScaleData(NP_TNL_epi, verbose = FALSE)
NP_TNL_epi <- RunPCA(NP_TNL_epi, verbose = FALSE)
ElbowPlot(NP_TNL_epi)
NP_TNL_epi <- RunUMAP(NP_TNL_epi, reduction = "pca", dims = 1:10)
NP_TNL_epi <- FindNeighbors(NP_TNL_epi, reduction = "pca", dims = 1:10)

DimPlot(NP_TNL_epi, reduction = "umap",split.by = "groups")

p1 <- DimPlot(NP_TNL_epi, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_epi, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
Idents(NP_TNL_epi) <- "groups"
DefaultAssay(NP_TNL_epi) <- "integrated"
NP_TNL_epi <- FindClusters(object = NP_TNL_epi, resolution = 0.1)
DimPlot(NP_TNL_epi, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")

NP_TNL_epi.markers <- FindAllMarkers(NP_TNL_epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


DefaultAssay(NP_TNL_epi) <- "SCT"
genes_to_check3 = c('HLA-G', 'KRT18', 'KRT8', 'GATA3','KRT17', 'KRT7')#滋养
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'DCN','LUM', 'OXTR', 'GJA1')#smc
genes_to_check3 = c('EPCAM','CD24','CCL21', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('EPCAM', 'CD24', 'KRT7', 'HLA-G')#滋养
genes_to_check3 = c('PARP1', 'MET', 'CDH1', 'HLA-G','PAPPA2', 'MMP11', 'ERVFRD-1', 'CYP19A1','CGA')#3种滋养
#villous cytotrophoblast cells [CTB, markers: PARP1 (15), MET (16),CDH1 (14)]
#extravillous trophoblast cells [EVT, markers:HLA-G (15), PAPPA2 (16), MMP11]
#syncytiotrophoblastcells [STB, markers: ERVFRD-1 (15), CYP19A1 (30) and CGA]
genes_to_check3 = c('EPCAM', 'CD24', 'KRT7','PARP1', 'MET', 'CDH1', 'HLA-G','PAPPA2', 'MMP11')#3种滋养
FeaturePlot(NP_TNL_epi, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "groups",keep.scale = "all")
VlnPlot(NP_TNL_epi, features = genes_to_check3, pt.size = 0, group.by = "endo_type")
saveRDS(NP_TNL_epi, file = "singlecellr/np_tnl/NP_TNL_epi20221113.rds")#####
NP_TNL_epi <-readRDS(file = "singlecellr/np_tnl/np_tnl_endo.rds")#####
