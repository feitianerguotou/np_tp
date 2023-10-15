library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(MySeuratWrappers)
library(SPOTlight)
library(DoubletFinder)
library(Libra)

NP_1 <- Read10X("single_cells/nptp/NP_2/")
NP_2 <- Read10X("single_cells/nptp/NP_3/")
NP_3 <- Read10X("single_cells/nptp/gse162122/061919/")
NP_4 <- Read10X("single_cells/nptp/gse162122/11564/")
NP_5 <- Read10X("single_cells/nptp/gse162122/11911/")
NP_6 <- Read10X("single_cells/nptp/gse162122/12640/")
NP_7 <- Read10X("single_cells/nptp/gse162122/12745/")

TNL_1 <- Read10X("single_cells/nptp/TNL_1/")
TNL_2 <- Read10X("single_cells/nptp/TNL_2/")
TNL_3 <- Read10X("single_cells/nptp/TNL_3/")
TNL_4 <- Read10X("single_cells/nptp/TNL_4/")
TNL_5 <- Read10X("single_cells/nptp/TNL_5/")

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
TNL_4object <- CreateSeuratObject(counts = TNL_4, project = "TNL_4",min.cells = 3, min.features = 200)
TNL_5object <- CreateSeuratObject(counts = TNL_5, project = "TNL_5",min.cells = 3, min.features = 200)

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
TNL_4object[["percent.mt"]] <- PercentageFeatureSet(TNL_4object, pattern = "^MT-")
TNL_5object[["percent.mt"]] <- PercentageFeatureSet(TNL_5object, pattern = "^MT-")

VlnPlot(NP_TNL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
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
TNL_4object <- subset(TNL_4object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_5object <- subset(TNL_5object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

NP_TNL <- c(NP_1object,NP_2object,NP_3object,NP_4object,NP_5object,NP_6object,NP_7object,TNL_1object,TNL_2object,TNL_3object,TNL_4object,TNL_5object)
names(NP_TNL) <- c('NP_1','NP_2','NP_3','NP_4','NP_5','NP_6','NP_7','TNL_1','TNL_2','TNL_3','TNL_4','TNL_5')

NP_TNL <- lapply(NP_TNL, function(x) {x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)})
NP_TNL <- lapply(NP_TNL, function(x) {x <- ScaleData(x, verbose = FALSE)})
NP_TNL <- lapply(NP_TNL, function(x) {x <- RunPCA(x, verbose = FALSE)})

Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:20, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  DoubletRate=ncol(data)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}

NP_TNL2 <- lapply(NP_TNL, function(x) {x <- Find_doublet(x)})

NP_TNL <- readRDS("single_cells/nptp/NP_TNL_doublefinder0905.rds")#####
saveRDS(NP_TNL2, file = "single_cells/nptp/NP_TNL2_doublefinder0905.rds")#####

NP_TNL3 <- lapply(NP_TNL2, function(x) {x <-subset(x,subset=doublet_info=="Singlet")})

features <- SelectIntegrationFeatures(object.list = NP_TNL3)
anchors <- FindIntegrationAnchors(object.list = NP_TNL3, anchor.features = features)
NP_TNL3 <- IntegrateData(anchorset = anchors)

saveRDS(NP_TNL3, file = "single_cells/nptp/NP_TNL3_doublefinder0905.rds")#####

DefaultAssay(NP_TNL3) <- "integrated"
NP_TNL3 <- ScaleData(NP_TNL3, verbose = FALSE)
NP_TNL3 <- RunPCA(NP_TNL3, verbose = FALSE)
ElbowPlot(NP_TNL3)
NP_TNL3 <- RunUMAP(NP_TNL3, reduction = "pca", dims = 1:20)
NP_TNL3 <- FindNeighbors(NP_TNL3, reduction = "pca", dims = 1:20)
DimPlot(NP_TNL3, reduction = "umap",split.by = "orig.ident")

NP_TNL3 <- FindClusters(NP_TNL3, resolution = 1)

DimPlot(NP_TNL3, reduction = "umap",cells.highlight = CellsByIdentities(NP_TNL3, idents = c(21:24)),cols.highlight = c("darkblue", "darkred","yellow","orange"))

plot <- DimPlot(NP_TNL3, reduction = "umap")
NP_TNL3 <- CellSelector(plot=plot,object = NP_TNL3,ident = "34")

DimPlot(NP_TNL3, reduction = "umap",label = T)

saveRDS(NP_TNL3, file = "single_cells/nptp/NP_TNL_double_clean.rds")#####
NP_TNL3 <- readRDS("single_cells/nptp/NP_TNL_double_clean.rds")#####

Idents(NP_TNL3) <- "integrated_snn_res.1"

genes_to_check3 = c('UBE2C', 'CDC20', 'TOP2A', 'STMN1','CD28', 'UCP2')#other
genes_to_check3 <- c('LUM','COL1A2','CD14','C1QB','TAGLN','ACTA2','TRAC','CD3E','TFF3','CCL21','VWF','PECAM1','KLRD1','KLRB1','TPSB2','TPSAB1','CD79A','CD79B','FCGR3B','NAMPT','HBB','HBA2','KRT7','PAPPA')
genes_to_check3 <- c('VWF', 'PECAM1','HLA-G','KRT7','DCN', 'LUM', 'TAGLN', 'ACTA2', 'HBB', 'HBA1','TRAC', 'CD3D','CD14', 'KLRD1','NCAM1', 'CD34')#all
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')#T
genes_to_check3 = c('CD163', 'CD86', 'CD14', 'PTPRC','IL1B', 'MRC1','S100A9','S100A8')#
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('KRT7','HLA-G')
genes_to_check3 = c('LUM','CD14','TAGLN','TRAC','TFF3','VWF','KLRD1','TPSB2','CD79A','FCGR3B','HBB','PAPPA2')
genes_to_check3 <- c('CD79A','VWF','KRT7','LUM','TFF3','TPSB2','CD14','FCGR3B','KLRD1','HBB','ACTA2','CD3D')
genes_to_check3 = c('DCN', 'LUM', 'TAGLN', 'ACTA2')
genes_to_check3 = c('CD14', 'TAGLN')
genes_to_check3 = c('DLK1', 'EGFL6', 'PRL', 'PTGDS')
VlnPlot(NP_TNL_dfclean,features = genes_to_check3, pt.size = 0,ncol = 1)
MySeuratWrappers::VlnPlot(NP_TNL,features = genes_to_check3,stacked=T,direction = "vertical",pt.size=0)+NoLegend()
MySeuratWrappers::VlnPlot(NP_TNL,features = top10$gene[51:60],stacked=T,direction = "vertical",pt.size=0)+NoLegend()
MySeuratWrappers::MultiFeaturePlot(NP_TNL, features = mapgene, reduction = "umap",min.cutoff = 0,ncol =8)
FeaturePlot(NP_TNL3, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol =2)
FeaturePlot(NP_TNL3, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "orig.ident")
DefaultAssay(NP_TNL3) <- "SCT"
NP_TNL_markers <- FindAllMarkers(NP_TNL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_markers %>% group_by(cluster) %>% top_n(n = 10)
write.csv(top10,file="singlecellr/np_top10.csv")

SMC=c(2,8,9,11,15,16,17,18)
Monocyte=c(0,19,20,21)
Epithelial_like_cell=c(5,13,33)
Endothelial_cell=c(1,3,14,22,24)
Fibroblast=c(4,10,31,32)
LEC=c(6)
Red_blood_cell=c(25,26)
Neutrophil=c(34)
NK_cell=c(12)
T_cell=c(7,23)
B_cell=c(30)
Mast_cell=c(29)
Unknown=c(27,28)

current.cluster.ids <- c(SMC,Monocyte,Neutrophil,NK_cell,T_cell,B_cell,Mast_cell,Unknown,
                         Endothelial_cell,
                         Fibroblast,
                         Epithelial_like_cell,
                         LEC,Red_blood_cell)

new.cluster.ids <- c(rep("SMC",length(SMC)),rep("Monocyte",length(Monocyte)), rep("Neutrophil",length(Neutrophil)),
                     rep("NK_cell",length(NK_cell)),
                     rep("T_cell",length(T_cell)),rep("B_cell",length(B_cell)),
                     rep("Mast_cell",length(Mast_cell)),rep("Unknown",length(Unknown)),
                     rep("Endothelial_cell",length(Endothelial_cell)),
                     rep("Fibroblast",length(Fibroblast)),
                     rep("Epithelial_like_cell",length(Epithelial_like_cell)),
                     rep("LEC",length(LEC)),
                     rep("Red_blood_cell",length(Red_blood_cell)))

NP_TNL3@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(NP_TNL3@active.ident)), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(NP_TNL3, reduction = "umap",label = T,group.by = "celltype")
DimPlot(NP_TNL3, reduction = "umap",label = T,group.by = "groups")
NP_TNL3@meta.data$groups[NP_TNL3@meta.data$orig.ident=="TNL_5"]=c("TP")

abc <- grep("pANN_",colnames(NP_TNL3@meta.data))
NP_TNL3@meta.data <- NP_TNL3@meta.data[,-abc]

NP_TNL3_clear <- subset(NP_TNL3, celltype == "Unknown", invert = T)

saveRDS(NP_TNL3_clear, file = "single_cells/nptp/NP_TNL_double_clean_unknown.rds")#####
NP_TNL3_clear <- readRDS("single_cells/nptp/NP_TNL_double_clean_unknown.rds")#####

DimPlot(NP_TNL3_clear, reduction = "umap",label = T,group.by = "celltype",split.by = "groups")
DimPlot(NP_TNL3_clear, reduction = "umap",label = T,group.by = "groups")

Idents(NP_TNL3_clear) <- "celltype"
DefaultAssay(NP_TNL3_clear) <- "SCT"
genes_to_check3 <- c('CD79A','VWF','KRT7','LUM','TFF3','TPSB2','CD14','FCGR3B','KLRD1','HBB','ACTA2','CD3D')
MySeuratWrappers::VlnPlot(NP_TNL3_clear,features = genes_to_check3,stacked=T,direction = "vertical",pt.size=0,group.by = "celltype")+NoLegend()+theme(axis.text.y = element_text(size = 5))

NP_TNL3_clear.markers <- FindAllMarkers(NP_TNL3_clear, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')
top10 <- NP_TNL3_clear.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]

NP_TNL3_clear_300 <- subset(NP_TNL3_clear, downsample = 300)
DefaultAssay(NP_TNL3_clear_300) <- "integrated" #integrated
DoHeatmap(NP_TNL3_clear_300, features = top10_order$gene, angle = 25,size = 4,group.by = "celltype")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))

table(Idents(NP_TNL3_clear), NP_TNL3_clear$groups)

NP_TNL3_clear$celltype.groups <- paste(Idents(NP_TNL3_clear),NP_TNL3_clear$groups,sep = "_")
Idents(NP_TNL3_clear) <- "celltype.groups"
# Endothelial_cell Fibroblast SMC Monocyte T_cell LEC NK_cell B_cell Epithelial_like_cell Mast_cell Neutrophil

FindMarkers(NP_TNL3_clear, ident.1 = "SMC_TP", ident.2 = "SMC_NP", verbose = FALSE, assay = 'RNA',slot = 'counts',test.use = "MAST",
            logfc.threshold = 0.01,only.pos =F,min.pct = 0.25) %>% write.csv(
              file="single_cells/nptp/SMC_des_TP_NP_mast.csv")

Idents(NP_TNL3_clear) <- "celltype"
mapgene <- c('TAGLN', 'DES', 'ACTG2', 'CNN1','MYLK', 'STMN1','GJA4', 'FTL','TIMP1', 'ATP1A1')#SMC
mapgene <- c('MYLK', 'OXTR','GJA1','PTGS2','PTGER3','ACTA2','PTGDS','GJA4')#shousuo
mapgene <- c('MDK', 'ADIRF', 'RPS27', 'ENO1','COL3A1', 'COL5A2', 'IGFBP2', 'IFI6','COL4A1', 'THBS1')#Fibroblasts
mapgene <- c('ALDOA', 'CRIP1', 'HLA-A','CD14', 'CD163', 'CCL4', 'SPP1','CCL3', 'ISG15','IFI6')#Monocyte
mapgene <- c('C1QB', 'C1QA', 'CCL3', 'TYROBP','RNASE1', 'CXCL14','IFITM1', 'FBLN1','SFRP1', 'PI16')#T_cells
mapgene <- c('NME2', 'ZFP36','FN1', 'CCN2', 'MARCKS', 'C1QA','IGFBP7','FN1', 'CCN2', 'C1QB')#Endothelial_cells
mapgene <- c('LYVE1', 'IFITM3', 'IGFBP6', 'GSN','AKAP12')#LEC
mapgene <- c('WFDC2', 'ESR1', 'SCGB2A1', 'MMP7','SOX17','BRI3','TIMP3', 'SERPINE2', 'PAPPA2','XAGE2')#Epithelial_like_cells

mapgene <- c('PTPRCAP', 'IFITM1','C1QB', 'C1QA', 'IFITM3','PTPRCAP', 'IFITM1','C1QB', 'C1QA', 'IFITM3')#NK_cells
mapgene <- c('EEF1G', 'PTPRCAP', 'IFI44L','HLA-DRB6', 'SPP1','FGL2', 'SLAMF7', 'ISG15','EEF1G', 'CPVL')#B_cells
mapgene <- c('SDPR', 'ATP5E', 'GNB2L1', 'SELENOW','LYVE1', 'GSN','IFI6', 'FBLN1','DCN', 'LUM')#LEC
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5','IDO1', 'CTLA4', 'HAVCR2', 'LAG3','IL10', 'CXCL13','SIRPG', 'CD27')#

VlnPlot(NP_TNL3_clear, features = mapgene, pt.size = 0,ncol = 10, split.by = "groups",idents = c("B_cell"), y.lab = "", x.lab = "")
MySeuratWrappers::VlnPlot(NP_TNL3_clear, features = mapgene, pt.size = 0,ncol =10, split.by = "groups", y.lab = "", x.lab = "")   
DoHeatmap(object=subset(NP_TNL3_clear, celltype == c("B_cell")), features = mapgene, angle = 25,size = 4,group.by = "groups",slot = "scale.data")
MySeuratWrappers::VlnPlot(NP_TNL3_clear,group.by = "celltype.groups",
                          features = genes_to_check3,stacked=T,sort = F,direction = "vertical",pt.size=0,axis.ticks.size = 1,axis.text.size = 12,axis.title.size=15)+NoLegend()+ theme(axis.text.y = element_text(size = 5))
DotPlot(NP_TNL3_clear,group.by = "celltype.groups",features = genes_to_check3,cols = c("lightgrey", "blue"))+theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ggplot2:::coord_flip()

#pseudobulk_DESeq2
head(NP_TNL3_clear@meta.data)
DE = run_de(NP_TNL3_clear, de_family = 'pseudobulk', de_method = 'DESeq2', de_type = 'LRT',n_threads = 20,cell_type_col = "celltype",replicate_col = "orig.ident",label_col = "groups")
head(DE)
write.csv(DE,file="single_cells/nptp/NP_TNL3_clear_deseq2.csv")

library(scRNAtoolVis)
markers <- read_csv(file="single_cells/nptp/NP_TNL3_clear_deseq2_2.csv")#adjp 001

jjVolcano(diffData = markers)
jjVolcano(diffData = markers,tile.col=mycolor2,log2FC.cutoff = 0.58492501,
          topGeneN = 3, legend.position = c(0.95,0.1))
mycolor <- c('#808000','Orange','#8A2BE2','#9ACD32','#D55E00','#F5DEB3','#CC79A7','#008B8B','#7FFFD4','#FF0000')
mycolor2 <- c('#808000','Orange','#8A2BE2','#9ACD32','#FF00FF','#D55E00','#F5DEB3','#CC79A7','#FFC0CB')

NP_v <- Load10X_Spatial("singlecellr/spaceranger_result/A1/",filename = "filtered_feature_bc_matrix.h5", assay = "Spatial",
                        slice = "slice1",
                        filter.matrix = TRUE)
TP_v <- Load10X_Spatial("singlecellr/spaceranger_result/D1/",filename = "filtered_feature_bc_matrix.h5", assay = "Spatial",
                        slice = "slice1",
                        filter.matrix = TRUE)
Idents(NP_v) <- "orig.ident"
Idents(TP_v) <- "orig.ident"

VlnPlot(TP_v, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
VlnPlot(NP_v, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p1 <- SpatialFeaturePlot(NP_v, features = "nCount_Spatial") 
p2 <- SpatialFeaturePlot(TP_v, features = "nCount_Spatial") 
p1+p2
NP_v <- SCTransform(NP_v, assay = "Spatial", verbose = FALSE)
TP_v <- SCTransform(TP_v, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(NP_v, features = c("CD163", "KRT8","TAGLN"))
SpatialFeaturePlot(TP_v, features = c("KRT7", "KRT8","PTPRC","KRT18"),alpha = c(0.1, 1))

NP_v <- RunPCA(NP_v, assay = "SCT", verbose = FALSE)
NP_v <- FindNeighbors(NP_v, reduction = "pca", dims = 1:30)
NP_v <- FindClusters(NP_v, verbose = FALSE, resolution = 0.8)
NP_v <- RunUMAP(NP_v, reduction = "pca", dims = 1:30)

TP_v <- RunPCA(TP_v, assay = "SCT", verbose = FALSE)
TP_v <- FindNeighbors(TP_v, reduction = "pca", dims = 1:30)
TP_v <- FindClusters(TP_v, verbose = FALS, Eresolution = 0.8)
TP_v <- RunUMAP(TP_v, reduction = "pca", dims = 1:30)

p1 <- DimPlot(TP_v, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(TP_v, label = TRUE, label.size = 3)
p1 + p2

p1 <- DimPlot(NP_v, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(NP_v, label = TRUE, label.size = 3)
p1 + p2

NP_v@reductions$spatial = NP_v@reductions$umap
NP_v@reductions$spatial@key = 'spatial_'
NP_v@reductions$spatial@cell.embeddings = as.matrix(NP_v@images$slice1@coordinates[,c(3,2)])
NP_v@reductions$spatial@cell.embeddings[,1] = -NP_v@reductions$spatial@cell.embeddings[,1]
colnames(NP_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(NP_v, features = c("CD86","CD163"), blend = T,pt.size = 1,cols = c('lightgrey','blue','red'),reduction = 'spatial',combine = T)+coord_flip()+NoLegend()&NoAxes()

TP_v@reductions$spatial = TP_v@reductions$umap
TP_v@reductions$spatial@key = 'spatial_'
TP_v@reductions$spatial@cell.embeddings = as.matrix(TP_v@images$slice1@coordinates[,c(3,2)])
TP_v@reductions$spatial@cell.embeddings[,1] = -TP_v@reductions$spatial@cell.embeddings[,1]
colnames(TP_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(TP_v, features = c("CD86","CD163"), blend = T,pt.size = 1, cols = c('lightgrey','blue','red'),reduction = 'spatial',combine = T)+coord_flip()+NoLegend()&NoAxes()

NP_v <- readRDS("single_cells/nptp/NP_v.rds")#####
TP_v <- readRDS("single_cells/nptp/TP_v.rds")#####
SpatialDimPlot(NP_v)

library(SPOTlight)
library(scater)
library(scran)
library(future)
set.seed(123)
plan("multisession", workers = 20)
options(future.globals.maxSize = 4000000000)

markers_all <- Seurat::FindAllMarkers(NP_TNL3_clear, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1,test.use = "roc", verbose = FALSE, assay = 'RNA',slot = 'counts')
markers_top10 <- markers_all %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
write.csv(markers_all,file="singlecellr/til_tnl/marker_TNL_TIL_13_roc.csv")
markers_all <- read.csv("singlecellr/til_tnl/marker_TNL_TIL_13_roc.csv")

sc_p <- DimPlot(TIL,reduction = 'umap',label = T,group.by = 'celltype.new')
st_p <- Seurat::SpatialDimPlot(TIL_v,label = T)
sc_p+st_p

set.seed(123)
NP_TNL3_clear_3000 <- subset(NP_TNL3_clear, downsample = 3000)
TIL_3000 <- subset(TNL, downsample = 3000)

spotlight_NP <- SPOTlight(x=NP_TNL3_clear_3000,
                           y = NP_v@assays$Spatial@counts,
                           groups = as.character(NP_TNL3_clear_3000$celltype),
                           mgs = markers_top10,
                           weight_id = "myAUC",
                           group_id = "cluster",
                           gene_id = "gene")
saveRDS(spotlight_NP, file = "single_cells/nptp/spotlight_NP.rds")

head(mat <- spotlight_NP$mat)[, seq_len(3)]
mod <- spotlight_NP$NMF

plotTopicProfiles(
  x = mod,
  y = NP_TNL3_clear_3000$celltype,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = TIL_300$celltype.new,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)


library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
plotCorrelationMatrix(mat)
plotSpatialScatterpie(
  x = TP_v,
  y = mat)
plotCorrelationMatrix(mat)
plotInteractions(mat, "network")

NP_TNL3_clear_200 <- subset(NP_TNL3_clear, downsample = 200)
NP_TNL3_clear_200 <- subset(NP_TNL3_clear_200, celltype == "Red_blood_cell", invert = T)
DefaultAssay(NP_TNL3_clear_200) <- "RNA"
sce <- as.SingleCellExperiment(NP_TNL3_clear_200)
sce <- logNormCounts(sce)
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))

dec <- modelGeneVar(sce, subset.row = genes)
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sce) <- colData(sce)$celltype
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
head(mgs_df)

idx <- split(seq(ncol(sce)), sce$celltype)

n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

#TP_v
res2 <- SPOTlight(
  x = sce,
  y = TP_v@assays$Spatial@counts,
  groups = sce$celltype,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

decon_mtrx <- res2$mat
colnames(decon_mtrx) <- paste(colnames(decon_mtrx),"Spotlight",sep = "_")
TP_v@meta.data <- cbind(TP_v@meta.data, decon_mtrx)
TP_v[["SPOTlight"]] <- CreateAssayObject(t(res2$mat))
DefaultAssay(TP_v) <- "SPOTlight"

head(mat <- res2$mat)[, seq_len(length(unique(sce$celltype)))]
mod <- res2$NMF

res.data <- (mat <- res2$mat)[, seq_len(length(unique(sce$celltype)))]
plotTopicProfiles(
  x = mod,
  y = sce$celltype,  #
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

ct <- colnames(mat)
mat[mat < 0.1] <- 0

paletteMartin <- c(
  "#000000", "#004949", "#009292", "#db6d00", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  "#920000", "#924900", "#24ff24", "#ffff6d","#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))


names(pal) <- ct

plotSpatialScatterpie(
  x = TP_v,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))

celltypes = rownames(TP_v)
SpatialFeaturePlot(TP_v, features = celltypes, 
                   pt.size.factor = 1.6, 
                   ncol = 4, 
                   crop = TRUE)

plotSpatialScatterpie(
  x = TP_v,
  y = mat,
  cell_types = colnames(mat),
  img = T, 
  scatterpie_alpha = 1,
  pie_scale = 0.4, 
  # Rotate the image 90 degrees counterclockwise
  degrees = -90,
  # Pivot the image on its x axis
  axis = "h") +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))


saveRDS(res2, file = "single_cells/nptp/spotlight_TP.rds")
res2 <- readRDS("single_cells/nptp/spotlight_NP.rds")


#NP_v
res <- SPOTlight(
  x = sce, 
  y = NP_v@assays$Spatial@counts,
  groups = as.character(sce$celltype),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

saveRDS(res, file = "single_cells/nptp/spotlight_NP.rds")
res <- readRDS("single_cells/nptp/spotlight_NP.rds")

decon_mtrx <- res$mat
colnames(decon_mtrx) <- paste(colnames(decon_mtrx),"Spotlight",sep = "_")
TP_v@meta.data <- cbind(TP_v@meta.data, decon_mtrx)
TP_v[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(TP_v) <- "SPOTlight"

head(mat <- res$mat)[, seq_len(length(unique(sce$celltype)))]
mod <- res$NMF

res.data <- (mat <- res$mat)[, seq_len(length(unique(sce$celltype)))]
plotTopicProfiles(
  x = mod,
  y = sce$celltype,  #
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

ct <- colnames(mat)
mat[mat < 0.1] <- 0

paletteMartin <- c(
  "#000000", "#004949", "#009292", "#db6d00", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  "#920000", "#924900", "#24ff24", "#ffff6d","#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))

names(pal) <- ct

plotSpatialScatterpie(
  x = NP_v,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))

celltypes = rownames(TP_v)
SpatialFeaturePlot(TP_v, features = celltypes, 
                   pt.size.factor = 1.6, 
                   ncol = 4, 
                   crop = TRUE)


#CellChat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(mindr)
library(NMF)
library(circlize)

NP <- subset(NP_TNL3_clear, groups == "NP")
NP <- subset(NP, celltype == "Red_blood_cell", invert = T)
NP <- SCTransform(NP)
DefaultAssay(NP) <- "SCT"
TP <- subset(NP_TNL3_clear, groups == "TP")
TP <- subset(TP, celltype == "Red_blood_cell", invert = T)
TP <- SCTransform(TP)
DefaultAssay(TP) <- "SCT"

cellchat_np <- createCellChat(object = NP, group.by = "celltype", assay = "RNA")
cellchat_tp <- createCellChat(object = TP, group.by = "celltype", assay = "RNA")
groupSize_np <- as.numeric(table(cellchat_np@idents))
groupSize_np <- as.numeric(table(cellchat_tp@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

CellChatDB_interaction <-CellChatDB$interaction
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

saveRDS(cellchat_tp, file = "single_cells/nptp/cellchat_tp.rds")#####
saveRDS(cellchat_np, file = "single_cells/nptp/cellchat_np.rds")#####
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


netVisual_diffInteraction(cellchat_np_merge, weight.scale = T,label.edge = T)

cellchat_np_merge <- identifyOverExpressedGenes(cellchat_np_merge, 
                                                pos.dataset = "TP", features.name = "TP", 
                                                only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
net <- netMappingDEG(cellchat_np_merge, features.name = "TP")
net.up <- subsetCommunication(cellchat_np_merge, net = net, datasets = "TP",ligand.logFC = 1, receptor.logFC = 1)
net.down <- subsetCommunication(cellchat_np_merge, net = net, datasets = "NP",ligand.logFC = -1, receptor.logFC = -1)
write.csv(net,"singlecellr/np_tnl/cellchat_np_merge202210des.csv")
netVisual_chord_gene(cellchat_tp, targets.use = "SMC",
                     slot.name = 'net', net = net.up ,lab.cex = 0.8, small.gap = 3.5)
netVisual_chord_gene(cellchat_tp, sources.use = "SMC",
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5)
netVisual_chord_gene(cellchat_np, targets.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5 )
netVisual_chord_gene(cellchat_np, sources.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.5, small.gap = 3.5,legend.pos.x = 1,legend.pos.y =1 )




#SMC
NP_TNL_SMC <- subset(NP_TNL3_clear, celltype == "SMC")
NP_TNL_SMC <- ScaleData(NP_TNL_SMC, verbose = FALSE)
NP_TNL_SMC <- RunPCA(NP_TNL_SMC, npcs = 30, verbose = FALSE)
NP_TNL_SMC <- RunUMAP(NP_TNL_SMC, reduction = "pca", dims = 1:20)
NP_TNL_SMC <- FindNeighbors(NP_TNL_SMC, reduction = "pca", dims = 1:20)

p1 <- DimPlot(NP_TNL_SMC, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_SMC, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DefaultAssay(NP_TNL_SMC) <- "integrated"
NP_TNL_SMC <- FindClusters(NP_TNL_SMC, resolution = 0.1)#
DimPlot(NP_TNL_SMC, reduction = "umap")
library(clustree)
clustree(NP_TNL_SMC@meta.data, prefix = "integrated_snn_res.")

NP_TNL_SMC <- FindClusters(NP_TNL_SMC, resolution = 0.5)
DefaultAssay(NP_TNL_SMC) <- "SCT"
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'ACTA1','GJA1')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA','OXTR', 'GJA1')#
genes_to_check3 <- c('TPM2', 'TES','VIM','ZEB2','UTRN', 'MTATP6P1', 'TMSB4X', 'XIST', 'VCAN', 'TFRC','TYROBP')#all
genes_to_check3 = c('UBE2C', 'CDC20', 'TOP2A', 'STMN1','CD28', 'UCP2')#other
genes_to_check3 = c('CCN1', 'ACTG2', 'ADIRF', 'SNCG','STEAP4', 'LHFPL6','C7', 'LUM')#marker
genes_to_check3 <- c('OXTR', 'GJA1','ELANE', 'IFNG','SERPINE1')
VlnPlot(NP_TNL_SMC,features = genes_to_check3, pt.size = 0,ncol = 2)
FeaturePlot(NP_TNL_SMC, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
DoHeatmap(object = NP_TNL_SMC, features = top10$gene)
DoHeatmap(NP_TNL_SMC, features = top10$gene, group.by = "SMC_type", angle = 25)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#标注
SMC_1=c(0)
SMC_2=c(1)
SMC_3=c(2)
SMC_4=c(3)
SMC_5=c(4)
current.cluster.ids <- c(SMC_1,
                         SMC_2,
                         SMC_3,
                         SMC_4,SMC_5)

new.cluster.ids <- c(rep("SMC_1",length(SMC_1)),
                     rep("SMC_2",length(SMC_2)),
                     rep("SMC_3",length(SMC_3)),
                     rep("SMC_4",length(SMC_4)),
                     rep("SMC_5",length(SMC_5)))

NP_TNL_SMC@meta.data$SMC_type <- plyr::mapvalues(x = as.integer(as.character(NP_TNL_SMC@meta.data$integrated_snn_res.0.1)), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(NP_TNL_SMC, reduction = "umap",group.by = "SMC_type", split.by = "groups")
DimPlot(NP_TNL_SMC, reduction = "umap",group.by = "SMC_type")

table(Idents(NP_TNL_SMC))
table(NP_TNL_SMC$SMC_type,NP_TNL_SMC$groups)

Idents(NP_TNL_SMC) <- "SMC_type"
DefaultAssay(NP_TNL_SMC) <- "SCT"
DefaultAssay(NP_TNL_SMC) <- "integrated"
NP_TNL_SMC.markers <- FindAllMarkers(NP_TNL_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')

top10 <- NP_TNL_SMC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(NP_TNL_SMC, features = top10_order$gene,slot = "scale.data", angle = 20,label = F,group.by = "SMC_type")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

NP_TNL_SMC_300 <- subset(NP_TNL_SMC, downsample = 300)
DefaultAssay(NP_TNL_SMC_300) <- "SCT" #integrated
DefaultAssay(NP_TNL_SMC_300) <- "RNA" #integrated
DoHeatmap(NP_TNL_SMC_300, features = top10_order$gene, angle = 25,size = 4,group.by = "SMC_type")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))



write.csv(NP_TNL_SMC_markers,file="singlecellr/NP_TNL_SMC_markers202210.csv")
VlnPlot(NP_TNL_SMC,features = top10$gene, pt.size = 0,ncol = 2)
SMC_marker <- c('STEAP4', 'THY1','RGS16','DES', 'ACTG2', 'SLMAP','BNC2','CLDN1', 'ISOC1', 'RERGL', 'NTRK2', 'SNCG','SPP1', 'SERPINE2', 'AOC1')
SMC_marker <- c('STEAP4', 'THY1','DES', 'SLMAP','BNC2','CLDN1', 'RERGL', 'NTRK2','SPP1', 'SERPINE2')
SMC_marker <- c('OXTR','GJA1','ELANE','IFNG','PTGFR','PTGS2')
VlnPlot(NP_TNL_SMC,features = SMC_marker, pt.size = 0,ncol = 3,group.by = "SMC_type")
MySeuratWrappers::VlnPlot(NP_TNL_SMC,features = SMC_marker,stacked=T,pt.size=0,group.by = "SMC_type",x.lab = "",y.lab = "",ncol = 3)
MySeuratWrappers::MultiFeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,ncol =8)
VlnPlot(NP_TNL_SMC,features = SMC_marker, pt.size = 0,ncol = 3)
FeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,keep.scale = "all",split.by = "groups",ncol = 4)
FeaturePlot(NP_TNL_SMC, features = SMC_marker, reduction = "umap",min.cutoff = 0,keep.scale = "all",ncol = 4)
DotPlot(NP_TNL_SMC,group.by = "SMC_type.groups",features = SMC_marker,cols = c("lightgrey", "blue"))+theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ggplot2:::coord_flip()

NP_TNL_SMC$SMC_type.groups <- paste(Idents(NP_TNL_SMC),NP_TNL_SMC$groups,sep = "_")
saveRDS(NP_TNL_SMC, file = "single_cells/nptp/NP_TNL_SMC2023.rds")#####
NP_TNL_SMC <- readRDS("single_cells/nptp/NP_TNL_SMC2023.rds")#####

NP_v <- readRDS("single_cells/nptp/NP_v.rds")#####
TP_v <- readRDS("single_cells/nptp/TP_v.rds")#####

NP_TNL_SMC_NP <- subset(NP_TNL_SMC, groups == "NP")
NP_TNL_SMC_TP <- subset(NP_TNL_SMC, groups == "TP")

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
SpatialFeaturePlot(NP_v_SMC, features = c("SMC-1","SMC-2","SMC-3","SMC-4","SMC-5"), alpha = c(0.1, 1), min.cutoff = 0, ncol = 5)

anchors_SMC_TP <- FindTransferAnchors(reference = NP_TNL_SMC_TP, query = TP_v, normalization.method = "SCT")
predictions.assay_SMC_TP <- TransferData(anchorset = anchors_SMC_TP, refdata = NP_TNL_SMC_TP$SMC_type, prediction.assay = TRUE, 
                                         weight.reduction = TP_v[["pca"]], dims = 1:30)
TP_v_SMC <- TP_v
TP_v_SMC[["predictions"]] <- predictions.assay_SMC_TP
DefaultAssay(TP_v_SMC) <- "predictions"
SpatialFeaturePlot(TP_v_SMC, features = c("SMC-1","SMC-2","SMC-3","SMC-4","SMC-5"), alpha = c(0.1, 1), min.cutoff = 0, ncol = 5)

DefaultAssay(NP_TNL_SMC) <- "integrated"
DoHeatmap(object = NP_TNL_SMC, features = top10$gene)



#Fibroblasts
NP_TNL_Fibroblasts <- subset(NP_TNL3_clear, celltype == "Fibroblast")
NP_TNL_Fibroblasts <- ScaleData(NP_TNL_Fibroblasts, verbose = FALSE)
NP_TNL_Fibroblasts <- RunPCA(NP_TNL_Fibroblasts, npcs = 30, verbose = FALSE)
NP_TNL_Fibroblasts <- RunUMAP(NP_TNL_Fibroblasts, reduction = "pca", dims = 1:20)
NP_TNL_Fibroblasts <- FindNeighbors(NP_TNL_Fibroblasts, reduction = "pca", dims = 1:20)

Idents(NP_TNL_Fibroblasts) <- "groups"
Idents(NP_TNL_Fibroblasts) <- "integrated_snn_res.0.1"
Idents(NP_TNL_Fibroblasts) <- "fibor"
DefaultAssay(NP_TNL_Fibroblasts) <- "integrated"
DimPlot(NP_TNL_Fibroblasts, reduction = "umap", group.by = "integrated_snn_res.0.2",split.by = "groups")

NP_TNL_Fibroblasts@meta.data$fibor 
NP_TNL_Fibroblasts@meta.data$fibor[NP_TNL_Fibroblasts@meta.data$integrated_snn_res.0.1=="3"]=c("DSC")
NP_TNL_Fibroblasts <- FindClusters(NP_TNL_Fibroblasts, resolution = 0.2)
DimPlot(NP_TNL_Fibroblasts, reduction = "umap", split.by = "groups",group.by = "fibor")

genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'MYH11','GJA1')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('APOD', 'MMP2', 'SPOCK1', 'PRG4','TPPP3', 'CD55','ACTA2','MYH11', 'ACTG2')
genes_to_check3 = c('SERPINE2', 'IGFBP6', 'RGCC', 'NR4A1','EGR1', 'NDRG2','MYH11','ACTA2', 'MYLK')
genes_to_check3 = c('DLK1','EGFL6', 'PRL','PTGDS','MFAP5','PCOLCE2', 'NDRG2','MTND1P23')
DefaultAssay(NP_TNL_Fibroblasts) <- "SCT"
VlnPlot(NP_TNL_Fibroblasts,features = genes_to_check3, pt.size = 0,ncol = 2)
FeaturePlot(NP_TNL_Fibroblasts, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all",)
FeaturePlot(NP_TNL_Fibroblasts, features = genes_to_check3, reduction = "umap",min.cutoff = 0, keep.scale = "all",split.by = "groups")
DotPlot(NP_TNL_Fibroblasts, features = genes_to_check3, group.by = "fibor",col.min=0)+ theme(axis.text.x=element_text(hjust = 1,vjust=0.5))+theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#加框

MySeuratWrappers::VlnPlot(NP_TNL_Fibroblasts, features = genes_to_check3, pt.size = 0,group.by = "fibor",ncol =2, y.lab = "", x.lab = "",stacked = T)+NoLegend()+theme(axis.text.y = element_text(size = 5))

NP_TNL_Fibroblasts_markers <- FindAllMarkers(NP_TNL_Fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NP_TNL_Fibroblasts_markers <- FindAllMarkers(NP_TNL_Fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')

top10 <- NP_TNL_Fibroblasts_markers %>% group_by(cluster) %>% top_n(n = 10, wt=avg_log2FC)
DoHeatmap(object = NP_TNL_Fibroblasts, features = top10$gene,angle = 30)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

NP_TNL_Fibroblasts_300 <- subset(NP_TNL_Fibroblasts, downsample = 300)
DefaultAssay(NP_TNL_Fibroblasts_300) <- "SCT" #integrated
DoHeatmap(NP_TNL_Fibroblasts_300, features = top10$gene, angle = 25,size = 4,group.by = "fibor")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(NP_TNL_Fibroblasts_300, features = top10_order$gene,slot = "scale.data", angle = 20,label = F,group.by = "fibor")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

saveRDS(NP_TNL_Fibroblasts, file = "single_cells/nptp/NP_TNL_Fibroblasts2023.rds")#####
NP_TNL_Fibroblasts <- readRDS("single_cells/nptp/NP_TNL_Fibroblasts2023.rds")#####

#(Endothelial_cells)
NP_TNL_endo <- subset(NP_TNL3_clear, celltype == c("Endothelial_cell"))
NP_TNL_endo <- ScaleData(NP_TNL_endo, verbose = FALSE)
NP_TNL_endo <- RunPCA(NP_TNL_endo, verbose = FALSE)
ElbowPlot(NP_TNL_endo)
NP_TNL_endo <- RunUMAP(NP_TNL_endo, reduction = "pca", dims = 1:20)
NP_TNL_endo <- FindNeighbors(NP_TNL_endo, reduction = "pca", dims = 1:20)
Idents(NP_TNL_endo) <- "celltype"
Idents(NP_TNL_endo) <- "fibor"
DefaultAssay(NP_TNL_endo) <- "integrated"
DimPlot(NP_TNL_endo, reduction = "umap",split.by = "groups")

p1 <- DimPlot(NP_TNL_endo, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")
DefaultAssay(NP_TNL_endo) <- "integrated"
NP_TNL_endo <- FindClusters(object = NP_TNL_endo, resolution = 0.1)

NP_TNL_endo@meta.data$endo_type[NP_TNL_endo@meta.data$integrated_snn_res.0.1=="1"]=c("Arterial")
NP_TNL_endo$endo_type.group.new <- paste(Idents(NP_TNL_endo),NP_TNL_endo$groups,sep = "_")
DimPlot(NP_TNL_endo, reduction = "umap", repel = TRUE,group.by = "endo_type",split.by = "groups")
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "endo_type.group.new")
table(Idents(TNL_TIL_B))
table(TNL_TIL_endo$endo_type,TNL_TIL_endo$groups)
Idents(NP_TNL_endo) <- "integrated_snn_res.0.2"
Idents(NP_TNL_endo) <- "endo_type"
DefaultAssay(TNL_TIL_endo) <- "SCT"
genes_to_check3 = c('PROX1','CCL21', 'HEY1','CXCL12', 'IGFBP3','CD36', 'CA4','ACKR1')#ymphatic ECs (LECs; CCL21, PROX1).arteries (HEY1, IGFBP3), capillaries (CD36, CA4), veins (ACKR1)
genes_to_check3 = c('ACKR1', 'IGFBP3','CD36')#
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
MySeuratWrappers::VlnPlot(NP_TNL_endo, features = genes_to_check3, pt.size = 0,group.by = "endo_type",ncol =2, y.lab = "", x.lab = "",stacked = T,split.by = "groups")+NoLegend()+theme(axis.text.y = element_text(size = 8))
DotPlot(NP_TNL_endo, features = genes_to_check3, group.by = "endo_type",col.min=0)+ theme(axis.text.x=element_text(hjust = 1,vjust=0.5))+theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#

VlnPlot(NP_TNL_endo, features = genes_to_check3, pt.size = 0, group.by = "endo_type")
saveRDS(NP_TNL_endo, file = "single_cells/nptp/np_tnl_endo2023.rds")#####
NP_TNL_endo <-readRDS(file = "single_cells/nptp/np_tnl_endo2023.rds")#####


NP_TNL_endo.markers <- FindAllMarkers(NP_TNL_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NP_TNL_endo_markers <- FindAllMarkers(NP_TNL_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')

top10 <- NP_TNL_endo_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_endo, features = top10$gene,group.by = "NP_TNL_endo")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#Venous Arterial

NP_TNL_endo$celltype.groups <- paste(Idents(NP_TNL_endo),NP_TNL_endo$groups,sep = "_")
Idents(NP_TNL_endo) <- "celltype.groups"
DimPlot(NP_TNL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "celltype.groups")
MySeuratWrappers::VlnPlot(NP_TNL_endo, features = genes_to_check3, pt.size = 0,group.by = "celltype.groups",ncol =2, y.lab = "", x.lab = "",stacked = T,
                          cols = c("#F8766D", "#F8768D","#00CFC4","#00BFC4"))+NoLegend()+theme(axis.text.y = element_text(size = 8))
DE = run_de(NP_TNL_endo, de_family = 'pseudobulk', de_method = 'DESeq2', de_type = 'LRT',n_threads = 20,cell_type_col = "endo_type",replicate_col = "orig.ident",label_col = "groups")


#(Epithelial_like_cells)
NP_TNL_epi <- subset(NP_TNL3_clear, celltype == c("Epithelial_like_cell"))
NP_TNL_epi <- ScaleData(NP_TNL_epi, verbose = FALSE)
NP_TNL_epi <- RunPCA(NP_TNL_epi, verbose = FALSE)
ElbowPlot(NP_TNL_epi)
NP_TNL_epi <- RunUMAP(NP_TNL_epi, reduction = "pca", dims = 1:20)
NP_TNL_epi <- FindNeighbors(NP_TNL_epi, reduction = "pca", dims = 1:20)

DimPlot(NP_TNL_epi, reduction = "umap",split.by = "groups")

p1 <- DimPlot(NP_TNL_epi, reduction = "umap", group.by = "groups")
p2 <- DimPlot(NP_TNL_epi, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
Idents(NP_TNL_epi) <- "integrated_snn_res.0.1"
DefaultAssay(NP_TNL_epi) <- "integrated"
DefaultAssay(NP_TNL_epi) <- "RNA"
NP_TNL_epi <- FindClusters(object = NP_TNL_epi, resolution = 0.1)
DimPlot(NP_TNL_epi, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")
NP_TNL_epi@meta.data$epi[NP_TNL_epi@meta.data$integrated_snn_res.0.1=="0"]=c("EVT")


NP_TNL_epi.markers <- FindAllMarkers(NP_TNL_epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NP_TNL_epi.markers <- FindAllMarkers(NP_TNL_epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')



top10 <- NP_TNL_epi.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


DefaultAssay(NP_TNL_epi) <- "SCT"
genes_to_check3 = c('HLA-G', 'KRT18', 'KRT8', 'GATA3','KRT17', 'KRT7')#
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'DCN','LUM', 'OXTR', 'GJA1')#smc
genes_to_check3 = c('EPCAM','CD24','CCL21', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#
genes_to_check3 = c('EPCAM', 'CD24', 'KRT7', 'HLA-G')#
genes_to_check3 = c('PARP1', 'MET', 'CDH1','PAGE4', 'HLA-G','PAPPA2', 'MMP11', 'ERVFRD-1', 'CYP19A1','CGA')#
genes_to_check3 = c('KRT7','KRT14', 'PAPPA','PARP1','ERVFRD-1')#
#villous cytotrophoblast cells [CTB, markers: PARP1 (15), MET (16),CDH1 (14)]
#extravillous trophoblast cells [EVT, markers:HLA-G (15), PAPPA2 (16), MMP11]
#syncytiotrophoblastcells [STB, markers: ERVFRD-1 (15), CYP19A1 (30) and CGA]
genes_to_check3 = c('EPCAM', 'CD24', 'KRT7','PARP1', 'MET', 'CDH1', 'HLA-G','PAPPA2', 'MMP11')#

Idents(NP_TNL_epi) <- "epi"
FeaturePlot(NP_TNL_epi, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
VlnPlot(NP_TNL_epi, features = genes_to_check3, pt.size = 0)
MySeuratWrappers::VlnPlot(NP_TNL_epi, features = genes_to_check3, pt.size = 0,group.by = "epi",ncol =2, y.lab = "", x.lab = "",stacked = T)+NoLegend()+theme(axis.text.y = element_text(size = 8))
DotPlot(NP_TNL_epi, features = genes_to_check3, group.by = "epi",col.min=0)+ theme(axis.text.x=element_text(hjust = 1,vjust=0.5))+theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#
NP_TNL_epi.markers <- FindAllMarkers(NP_TNL_epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')
top10 <- NP_TNL_epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]

NP_TNL_epi_300 <- subset(NP_TNL_epi, downsample = 300)
DefaultAssay(NP_TNL_epi_300) <- "SCT" #integrated
DoHeatmap(NP_TNL_epi_300, features = top10_order$gene, angle = 25,size = 4,group.by = "epi")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))

saveRDS(NP_TNL_epi, file = "single_cells/nptp/NP_TNL_epi2023.rds")#####
NP_TNL_epi <-readRDS(file = "single_cells/nptp/NP_TNL_epi2023.rds")#####

NP_TNL_monocytic <- subset(NP_TNL3_clear, celltype == "Monocyte")
NP_TNL_monocytic <- ScaleData(NP_TNL_monocytic, verbose = FALSE)
NP_TNL_monocytic <- RunPCA(NP_TNL_monocytic, npcs = 30, verbose = FALSE)
NP_TNL_monocytic <- RunUMAP(NP_TNL_monocytic, reduction = "pca", dims = 1:20)
NP_TNL_monocytic <- FindNeighbors(NP_TNL_monocytic, reduction = "pca", dims = 1:20)
NP_TNL_monocytic <- FindClusters(NP_TNL_monocytic, resolution = 0.4)
Idents(NP_TNL_monocytic) <- "integrated_snn_res.0.1"
DefaultAssay(NP_TNL_monocytic) <- "integrated"
DefaultAssay(NP_TNL_monocytic) <- "SCT"
DimPlot(NP_TNL_monocytic, reduction = "umap", group.by = "integrated_snn_res.0.1")
DimPlot(NP_TNL_monocytic, reduction = "umap", group.by = "macro",split.by = "groups")

NP_TNL_monocytic@meta.data$macro[NP_TNL_monocytic@meta.data$integrated_snn_res.0.1=="2"]=c("TAGLN+")

NP_TNL_monocytic <- PrepSCTFindMarkers(NP_TNL_monocytic)
NP_TNL_monocytic.markers <- FindAllMarkers(NP_TNL_monocytic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NP_TNL_monocytic.markers <- FindAllMarkers(NP_TNL_monocytic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST", verbose = FALSE, assay = 'RNA',slot = 'counts')
top10 <- NP_TNL_monocytic.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

Idents(NP_TNL_monocytic) <- "macro"
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#
genes_to_check3 = c('CD163', 'CD86', 'CD14', 'CD68','CD209','IL1B', 'MRC1','S100A9','S100A8','CIITA','ITGAX','FCGR3A','SPP1', 'TREM2','FOLR2','FN1','FCN1','TAGLN')#
genes_to_check3 = c('SPP1', 'TREM2','FOLR2','FN1','FCN1','TAGLN')#
genes_to_check3 = c('CD14', 'FCN1','FOLR2','TAGLN')
genes_to_check3 = c('CD14','CD163', 'MRC1','CD209','CD86','IL1B','FOLR2','TREM2','SPP1')#
FeaturePlot(NP_TNL_monocytic, features = top10$gene, reduction = "umap",min.cutoff = 0,ncol=6)
FeaturePlot(NP_TNL_monocytic, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol=3)
VlnPlot(NP_TNL_monocytic,features = genes_to_check3, pt.size = 0,ncol = 6, combine = TRUE)
MySeuratWrappers::VlnPlot(NP_TNL_monocytic, features = genes_to_check3, pt.size = 0,ncol =2, y.lab = "", x.lab = "",stacked = T)+NoLegend()+theme(axis.text.y = element_text(size = 8))
saveRDS(NP_TNL_monocytic, file = "single_cells/nptp/NP_TNL_monocytic2023.rds")#####
NP_TNL_monocytic <- readRDS(file = "single_cells/nptp/NP_TNL_monocytic2023.rds")#####
#("CD163", "CD86"))
SpatialFeaturePlot(NP_v, features = c("CD68"))+ theme(legend.position="right")
SpatialFeaturePlot(TP_v, features = c("CD68"))+ theme(legend.position="right")



#Tcells
NP_TNL_TNK <- subset(NP_TNL3_clear, idents = c("T_cell"))
NP_TNL_TNK <- ScaleData(NP_TNL_TNK, verbose = FALSE)
NP_TNL_TNK <- RunPCA(NP_TNL_TNK, verbose = FALSE)
ElbowPlot(NP_TNL_TNK)
NP_TNL_TNK <- RunUMAP(NP_TNL_TNK, reduction = "pca", dims = 1:20)
NP_TNL_TNK <- FindNeighbors(NP_TNL_TNK, reduction = "pca", dims = 1:20)
Idents(NP_TNL_TNK) <- "celltype"
DimPlot(NP_TNL_TNK, reduction = "umap",group.by = "groups")

DimPlot(NP_TNL_TNK, reduction = "umap", group.by = "groups",split.by = "groups")

DefaultAssay(NP_TNL_TNK) <- "integrated"
NP_TNL_TNK <- FindClusters(object = NP_TNL_TNK, resolution = 0.1)
DimPlot(NP_TNL_TNK, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.1")


genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A')
genes_to_check3 = c('CD4','CD8A','ACTA2','KLRD1','MAFB','MKI67','KIT','FOXP3')#top
FeaturePlot(NP_TNL_TNK, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "groups")

DefaultAssay(NP_TNL_TNK) <- "SCT"
NP_TNL_TNK <- PrepSCTFindMarkers(NP_TNL_TNK)
NP_TNL_TNK.markers <- FindAllMarkers(NP_TNL_TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- NP_TNL_TNK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

FeaturePlot(NP_TNL_TNK, features = top10$gene, reduction = "umap",min.cutoff = 0,ncol=6)


library(monocle)
devtools::load_all("/home/data/vip01/R/x86_64-pc-linux-gnu-library/4.2/monocle/monocle/")
Idents(NP_TNL3) <- "celltype"
NP_TNL_SMC <- subset(NP_TNL3, idents = c("SMC"))

NP_TNL_SMC <- readRDS("single_cells/nptp/NP_TNL_SMC2023.rds")#####
Idents(NP_TNL_SMC) <- "SMC_type"
data <- GetAssayData(NP_TNL_SMC, assay = 'RNA', slot = 'counts')
cell_metadata <- NP_TNL_SMC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)


data_matrix <- as(as.matrix(NP_TNL_SMC@assays$RNA@data), 'sparseMatrix')
feature_ann <- data.frame(gene_id=rownames(data_matrix),gene_short_name=rownames(data_matrix))
rownames(feature_ann)<-rownames(data_matrix)
data_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-NP_TNL_SMC@meta.data
rownames(sample_ann)<-colnames(data_matrix)
data_pd <- new("AnnotatedDataFrame", data =sample_ann)
data_cds <- newCellDataSet(data_matrix,phenoData =data_pd, featureData =data_fd, expressionFamily=negbinomial.size())



#2
expr_matrix <- as(as.matrix(NP_TNL_SMC@assays$RNA@counts), 'sparseMatrix')
p_data <- NP_TNL_SMC@meta.data

p_data$celltype <- NP_TNL_SMC@active.ident

f_data <- data.frame(gene_short_name = row.names(expr_matrix),row.names = row.names(expr_matrix))

pd <- new('AnnotatedDataFrame', data = p_data)

fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) 

Idents(NP_TNL_SMC) <- "groups"
NP_TNL_SMC <- PrepSCTFindMarkers(NP_TNL_SMC)
deg.cluster <- FindAllMarkers(NP_TNL_SMC)

express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 *dispersion_fit)$gene_id
cds <- setOrderingFilter(cds,disp.genes)


plot_cell_trajectory(cds,cell_size = 0.1,color_by = "Pseudotime")


saveRDS(cds, file = "single_cells/nptp/SMC_data2_cds.rds")#######
cds <- readRDS(file = "single_cells/nptp/SMC_data2_cds.rds")#######

p1 <- plot_cell_trajectory(cds,cell_size = 0.3,color_by = "groups")
p2 <- plot_cell_trajectory(cds,cell_size = 0.3,color_by = "SMC_type")
p3 <- plot_cell_trajectory(cds,cell_size = 0.3,color_by = "Pseudotime")
p1|p2|p3

plot_complex_cell_trajectory(cds, x = 1, y = 2, 
                             color_by = "groups")
library(ggpubr)
df <- pData(cds)
ggplot(df, aes(Pseudotime, colour = SMC_type, fill=SMC_type)) + geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()


unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
ordering_genes = unsup_clustering_genes$gene_id
plot_pseudotime = plot_pseudotime_heatmap(cds[ordering_genes, ], num_clusters = 3,cores = 2, return_heatmap = T, show_rownames = F)

clusters <- cutree(plot_pseudotime$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.csv(clustering,file="single_cells/nptp/smc-pse-cluster")
test_genes <- c('ACTG2', 'CCL2', 'CCL19', 'CCL21','DES', 'EGR1','JUN', 'JUNB', 'MYLK', 'MYH11','SERPINE1', 'SLC2A3')
plot_genes_in_pseudotime(cds[test_genes,],
                         color_by = "groups",
                         cell_size=0.3,
                         ncol = 4)

NP_TNL_SMC <- readRDS("single_cells/nptp/NP_TNL_SMC2023.rds")#####
Idents(NP_TNL_SMC) <- "SMC_type"






