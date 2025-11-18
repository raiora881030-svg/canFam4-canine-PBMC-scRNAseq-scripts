##############################################
# Script: Canine PBMC scRNA-seq Analysis
# Author: Myung-Chul Kim (mck@knu.ac.kr)
#
# Description:
#   This script contains the analysis workflow used 
#   to process and visualize single-cell RNA-seq 
#   data from healthy dogs.
#
#   Included sections:
#     - Loading required packages
#     - Loading Seurat objects
#     - Subset-specific analyses (CD4, CD8/NKT, B cells, Myeloid)
#     - Visualization (UMAP/t-SNE, FeaturePlot, VlnPlot)
#     - Marker gene analysis and heatmaps
#     - CanFam3 vs CanFam4 comparison
#
# Notes:
#   - Absolute file paths have been removed for privacy.
#   - Script preserves original exploratory flow.
#   - Update working directory paths before use.
#
# Last Updated: 2025-XX-XX
##############################################

# mapping the healthy circulating leuokcytes in 6 dogs
# setwd("PATH_TO_PROJECT")
# setwd("PATH_TO_PROJECT")
# setwd("PATH_TO_PROJECT")

#
# --------------------------------------------
# Load required packages
# --------------------------------------------
library(GSVA)
library(ggrepel)  
library(monocle3)
library(Seurat)
library(Nebulosa)
library(CellChat)
library(ggpubr)
library(scDblFinder)
library(scRNAseq)
library(dplyr)
library(SingleR)
library(ggplot2)
library(dittoSeq)
library(escape)
library(scales)
library(RColorBrewer)
library(viridis)
library(harmony)
library(scRepertoire)
library(RColorBrewer)

# load the rds Seuerat object canfam3

# --------------------------------------------
# Load Seurat objects
# --------------------------------------------
Dog.combined.3 <- readRDS("Dog.combined.healthy.V5.3.7.rds")
Dog.combined.singlet.3 <- readRDS("Dog.combined.singlet.final.rds")
Dog.combined.doublet <- readRDS("Dog.combined.doublet.rds")
Dog.combined.singlet.CD4 <- readRDS("Dog.combined.singlet.CD4.rds")
Dog.combined.singlet.CD8NKT <- readRDS("Dog.combined.singlet.CD8NKT.rds")
Dog.combined.singlet.B <- readRDS("Dog.combined.singlet.B.rds")
Dog.combined.singlet.myeloid <- readRDS("Dog.combined.singlet.myeloid.rds")

# canfam4
Dog.combined <- readRDS("PATH_TO_FILE")
Dog.combined.singlet <- readRDS("PATH_TO_FILE")
Dog.combined.doublet <- readRDS("PATH_TO_FILE")
Dog.combined.singlet.T <- readRDS("PATH_TO_FILE")
Dog.combined.singlet.myeloid <- readRDS("PATH_TO_FILE")
Dog.combined.singlet.CD4 <- readRDS("PATH_TO_FILE")
Dog.combined.singlet.B <- readRDS("PATH_TO_FILE")
Dog.combined.singlet.CD8NKT <- readRDS("PATH_TO_FILE")

saveRDS(Dog.combined.singlet.CD4, "PATH_TO_FILE")
saveRDS(Dog.combined.singlet.B, "PATH_TO_FILE")

# age 
unique(Dog.combined.singlet$Subset5) 


#
cellchat.Dog <- readRDS("Data_Healthy/cellchat/cellchat.Dog.rds")
cellchat.Dog <- readRDS("Data_Healthy/cellchat/cellchat.Dog.label.rds.rds")
cellchat.Dog <- readRDS("Data_Healthy/cellchat/cellchat.Dog.label.T.rds")

#
ES.immune <- readRDS("Data_Healthy/ES.immune.singlet.CD8NKT.rds")
Dog.combined.singlet.CD8NKT <-AddMetaData(Dog.combined.singlet.CD8NKT, ES.immune)

ES.immune <- readRDS("Data_Healthy/ES.immune.singlet.CD4.rds")
Dog.combined.singlet.CD4 <-AddMetaData(Dog.combined.singlet.CD4, ES.immune)
View(Dog.combined.singlet.myeloid[[]])

ES.immune2 <- readRDS("Data_Healthy/ES.immune.singlet.myeloid.rds")
Dog.combined.singlet.myeloid <-AddMetaData(Dog.combined.singlet.myeloid, ES.immune2)

#
VlnPlot(Dog.combined.singlet.CD4, features = "CTLA4",sort = T) 
VlnPlot(Dog.combined.singlet.CD8NKT, features = "CTLA4",sort = T) 

#
VlnPlot(Dog.combined.singlet.CD4, features = "Exhaustion", group.by = "PD1", cols = c("#440154FF","#FDE725FF")) #+  stat_compare_means("t.test") # atomic vector issue
p <-VlnPlot(Dog.combined.singlet.CD4, features = "Exhaustion", group.by = "PD1", cols = c("#440154FF","#FDE725FF"))
p + stat_compare_means(method = "t.test")
ggsave("Data_Healthy/CD4/PD1+_CD4_Exhaustion.png", width = 4, height = 4)

VlnPlot(Dog.combined.singlet.CD8NKT, 
        features = c("Proinflammatory","Antiinflammatory","Exhaustion","T_Cell_Terminal_Differentiation"),
        group.by = "PD1", 
        cols = c("#440154FF","#FDE725FF"))
p <-VlnPlot(Dog.combined.singlet.CD8NKT, 
            features = "Anergy", 
            group.by = "PD1", 
            cols = c("#440154FF","#FDE725FF"))
p + stat_compare_means(method = "t.test")
ggsave("Data_Healthy/CD8/PD1+_CD8_AdaptiveImmunity.png", width = 3, height = 3)

pairwise <- pairwise.wilcox.test(df$value, df$group, p.adjust.method = "BH")


# Exhaustion status of Lag3
VlnPlot(Dog.combined.singlet.CD8NKT, features = "Exhaustion", group.by = "PD1_LAG3")
VlnPlot(subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos")), 
        features = "Exhaustion", group.by = "PD1_LAG3") +  stat_compare_means("t.test") # atomic vector issue

p <- VlnPlot(subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos")),
             features = "T_Cell_Terminal_Differentiation", group.by = "PD1_LAG3", c("#440154FF","#FDE725FF"))
p + stat_compare_means(method = "t.test")
ggsave("Data_Healthy/Lag3Exhaustion/PD1_TIM3_CD4_Exhaustion.png", width = 3, height = 3)

VlnPlot(subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos")), 
        features = "Exhaustion", group.by = "PD1_CTLA4") +  stat_compare_means("t.test") # atomic vector issue

p <- VlnPlot(subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos")),
             features = "Exhaustion", group.by = "PD1_CTLA4", c("#440154FF","#FDE725FF"))
p + stat_compare_means(method = "t.test")
ggsave("Data_Healthy/Lag3Exhaustion/PD1_CTLA4_CD8T_Exhaustion.png", width = 3, height = 3)

#
#
FeaturePlot(Dog.combined.singlet.CD4, features = "CCR4", order = T)
Dog.combined.singlet.CD4.Treg <- subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(7,17))

DefaultAssay(Dog.combined.singlet.CD4.Treg) <- "integrated"
Dog.combined.singlet.CD4.Treg <-ScaleData(Dog.combined.singlet.CD4.Treg.CD4.Treg, verbose = FALSE)
Dog.combined.singlet.CD4.Treg <-RunPCA(Dog.combined.singlet.CD4.Treg, npcs = 30, verbose = FALSE)
Dog.combined.singlet.CD4.Treg <-RunUMAP(Dog.combined.singlet.CD4.Treg, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4.Treg <- RunTSNE(Dog.combined.singlet.CD4.Treg, dims = 1:30, reduction = "pca")
Dog.combined.singlet.CD4.Treg <- FindNeighbors(Dog.combined.singlet.CD4.Treg, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4.Treg <- FindClusters(Dog.combined.singlet.CD4.Treg, resolution = 1)
DefaultAssay(Dog.combined.singlet.CD4.Treg) <- "RNA"
Dog.combined.singlet.CD4.Treg <-ScaleData(Dog.combined.singlet.CD4.Treg, verbose = FALSE)
Dog.combined.singlet.CD4.Treg <- JoinLayers(Dog.combined.singlet.CD4.Treg)
DefaultAssay(Dog.combined.singlet.CD4.Treg) <- "RNA"

DimPlot(Dog.combined.singlet.CD4.Treg, reduction = "tsne")
VlnPlot(Dog.combined.singlet.CD4.Treg, features = c("CCR4"))
FeaturePlot(Dog.combined.singlet.CD4.Treg, features = c("CCR4"), order = T)

Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters
Dog.markers <- FindAllMarkers(Dog.combined.singlet, only.pos = TRUE, min.pct = 0.36, logfc.threshold = 0.36) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top10, "Data_Healthy/CD4/CCR4_DEG_Top10.csv",append = "F",sep = ",")
ggsave("Data_Healthy/Tsne_Labeled_Screening_3.5.png", width=7, height=5)

write.table(top10, "grant/Top10.csv",append = "F",sep = ",")

DoHeatmap(subset(Dog.combined.singlet.CD4.Treg, downsample = 50),
          features = top5$gene, 
          #disp.max = 2,
          # disp.min = -2,
          size = 3, angle = 50, cells = 1:100)+
  theme(text = element_text(size=8)) + 
  scale_fill_viridis(option = "E") + 
  labs(x = "Viridis E", y = NULL) #+ NoLegend()

#+ scale_fill_viridis(100)
#+ NoLegend()
ggsave("Data_Healthy/Heatmap_Labeled_Singlet_legend.png", width=8, height=7)

Dog.combined.singlet.CD4.Treg$ccr4 <- Dog.combined.singlet.CD4.Treg@assays$RNA$data["CCR4",]
Dog.combined.singlet.CD4.Treg$CCR4_Sep <- ifelse(Dog.combined.singlet.CD4.Treg$ccr4>0, "Pos","Neg")
FeaturePlot(Dog.combined.singlet.CD4.Treg, features = c("CCR4"), split.by = "CCR4_Sep", order = T)

Dog.combined.singlet.CD4$CTLA4_Sep <- ifelse(Dog.combined.singlet.CD4$ctla4>0, "Pos","Neg")
FeaturePlot(Dog.combined.singlet.CD4, features = c("CTLA4"), split.by = "CTLA4_Sep", order = T)
Idents(Dog.combined.singlet.CD4) <- Dog.combined.singlet.CD4$CTLA4_Sep
CD4_CTLA4_marker <- FindMarkers(Dog.combined.singlet.CD4, ident.1 = "Pos", ident.2 = "Neg")
write.table(CD4_CTLA4_marker, "Data_Healthy/CD4_CTLA4_marker.csv", append = F, sep =",")

Dog.combined.singlet.CD8NKT$ctla4 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["CTLA4",]
Dog.combined.singlet.CD8NKT$CTLA4_Sep <- ifelse(Dog.combined.singlet.CD8NKT$ctla4>0, "Pos","Neg")
FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("CTLA4"), split.by = "CTLA4_Sep", order = T)
Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$CTLA4_Sep
CD8NKT_CTLA4_marker <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = "Pos", ident.2 = "Neg")
write.table(CD8NKT_CTLA4_marker, "CD8_CTLA4_marker.csv", append = F, sep =",")

Dog.combined.singlet.T$ctla4 <- Dog.combined.singlet.T@assays$RNA$data["CTLA4",]
Dog.combined.singlet.T$CTLA4_Sep <- ifelse(Dog.combined.singlet.T$ctla4>0, "Pos","Neg")
FeaturePlot(Dog.combined.singlet.T, features = c("CTLA4"), split.by = "CTLA4_Sep", order = T)
Idents(Dog.combined.singlet.T) <- Dog.combined.singlet.T$CTLA4_Sep
T_CTLA4_marker <- FindMarkers(Dog.combined.singlet.T, ident.1 = "Pos", ident.2 = "Neg")
write.table(T_CTLA4_marker, "T_CTLA4_marker.csv", append = F, sep =",")


Idents(Dog.combined.singlet.CD4.Treg) <- Dog.combined.singlet.CD4.Treg$CCR4_Sep
CCR4marker <- FindMarkers(Dog.combined.singlet.CD4.Treg, ident.1 = "Pos", ident.2 = "Neg")
write.table(CCR4marker, "CCR4marker.csv", append = F, sep =",")

Dog.combined.singlet.CD4$ccr4 <- Dog.combined.singlet.CD4@assays$RNA$data["CCR4",]
Dog.combined.singlet.CD4$CCR4_Sep <- ifelse(Dog.combined.singlet.CD4$ccr4>0, "Pos","Neg")
FeaturePlot(Dog.combined.singlet.CD4.Treg, 
            features = c("FOXP3","IL2RA","CTLA4","CCR4"), ncol = 1,
            #split.by = "CCR4_Sep",
            order = T,
            cols = c("#E9F0FF","#0A5096"))
ggsave("Data_Healthy/CD4/CCR4_Gene_expression.png", width = 6, height = 5)

#
View(Dog.combined.singlet.B[[]])
DimPlot(Dog.combined.singlet.B, label = T, reduction = "umap")
FeaturePlot(Dog.combined.singlet.B, label = F, #pt.size = 0.8, 
            features = c("B_class_switched"), 
            order = F,  cols = viridis(100)# cols = c("#ADDAE6","#E63222")
) #+NoLegend() 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Feature_B_T1_IFN_legend.png", width = 5, height = 5) 

DimPlot(Dog.combined.singlet.B, label = T, reduction = "umap")
FeaturePlot(Dog.combined.singlet.B, label = F, #pt.size = 0.8, 
            features = c("OSBPL10"), 
            order = T, cols = c("#ADDAE6","#E63222")
) #+NoLegend() 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Feature_B_T1_IFN_legend.png", width = 5, height = 5) 

plot_density(Dog.combined.singlet.B, features = c("Tumor_Immune_Escape"), 
             reduction = "umap", pal = "inferno") + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/Density_CD8_PDCD1.png", width = 5, height = 5) 
View(Dog.combined.singlet.B[[]])
#
DimPlot(Dog.combined.singlet.CD8NKT, label = T, reduction = "umap")
FeaturePlot(Dog.combined.singlet.CD8NKT, label = F, #pt.size = 0.8, 
            features = c("OSBPL10"), 
            order = T,  cols = c("#ADDAE6","#E63222")
            # cols = viridis(100)
) +NoLegend() 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/Feature_ZEB2.png", width = 5, height = 5) 
library(Nebulosa)
plot_density(Dog.combined.singlet.CD8NKT, features = c("Tumor_Immune_Escape"), 
             reduction = "umap", pal = "inferno") + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/Density_CD8_PDCD1.png", width = 5, height = 5) 

VlnPlot(Dog.combined.singlet.CD8NKT, features = c("CD8_memory","CD8_effector"), sort = T)

FeaturePlot(Dog.combined.singlet.CD4, label = F, #pt.size = 0.8, 
            features = c("Response_To_ICB"), 
            order = T,  cols = c("#ADDAE6","#E63222")
            # cols = viridis(100)
) +NoLegend() 

FeaturePlot(Dog.combined.singlet.myeloid, label = F, #pt.size = 0.8, 
            features = c("CD16A"), 
            order = T,  cols = c("#ADDAE6","#E63222")
            # cols = viridis(100)
) +NoLegend() 

# prop
library(speckle)
library(ggplot2)
library(limma)

table(Idents(Dog.combined.singlet))
table(Dog.combined.singlet$orig.ident)
prop.table(table(Idents(Dog.combined.singlet)))
proptable <- prop.table(table(Idents(Dog.combined.singlet), Dog.combined.singlet$orig.ident), margin = 2)
write.table(proptable, "Data_Healthy/Prop_Table.csv", append = F, sep = ",")

color_pal <- scales::hue_pal()(51)
dittoBarPlot(
  object = Dog.combined.singlet,
  var = "seurat_clusters",
  group.by = "orig.ident") #+ NoLegend()
ggsave("Data_Healthy/Prop_Plot.png", width = 5, height = 6)
#
FeaturePlot(Dog.combined.singlet, 
            features = c("DPYD"),pt.size = 0.8,
            reduction = "tsne", order = T,
            cols = c("#ADDAE6","#E63222")
) +NoLegend()
ggsave("Data_Healthy/CD4_CD274.png", width=3, height = 3)

VlnPlot(Dog.combined.singlet.CD4,  features = c("Exhaustion","T_Cell_Terminal_Differentiation"), 
        flip = T, stack = F, sort = T) + NoLegend()
ggsave("Data_Healthy/CD4/Exhaustion.png", width=5, height = 3)

FeaturePlot(Dog.combined.singlet.CD4, 
            features = c("Exhaustion","T_Cell_Terminal_Differentiation"),
            blend = T, pt.size = 1, cols = c("#ADDAE6", "#7E6097", "#FF4B20"),
            reduction = "umap", order = T) +NoLegend()

FeaturePlot(Dog.combined.singlet.CD4, 
            features = c("IL2RA"),pt.size = 0.8,
            reduction = "umap", order = T,
            cols = c("#ADDAE6","#E63222")
            #cols = viridis(100)
) +NoLegend()
ggsave("Data_Healthy/CD4/CD4_IL2RA.png", width=5, height = 5)


View(Dog.combined.singlet.CD4[[]])
# CanFam3 and CanFam4
Dog.BIG.3 <- readRDS("PATH_TO_FILE")
Dog.BIG.4 <- readRDS("PATH_TO_FILE")

VlnPlot(Dog.BIG.3, features = c("nFeature_RNA"), 
        cols = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6","#7E6097"), pt.size = 0) + NoLegend()
ggsave("PATH_TO_FILE", width = 3, height = 3)

VlnPlot(Dog.BIG.4, features = c("percent.mt"), 
        cols = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6","#7E6097"), pt.size = 0.1) + NoLegend()
ggsave("PATH_TO_FILE", width = 3, height = 3)

Dog.BIG.3.1 <- subset(Dog.BIG.3, Subset8 %in% c("Gyeomdong")) 
Dog.BIG.3.2 <- subset(Dog.BIG.3, Subset8 %in% c("Makdoong"))
Dog.BIG.3.3 <- subset(Dog.BIG.3, Subset8 %in% c("Batt"))
Dog.BIG.3.4 <- subset(Dog.BIG.3, Subset8 %in% c("Koko"))
Dog.BIG.3.5 <- subset(Dog.BIG.3, Subset8 %in% c("Sojin"))
Dog.BIG.3.6 <- subset(Dog.BIG.3, Subset8 %in% c("Bori"))

Dog.BIG.4.1 <- subset(Dog.BIG.4, Subset8 %in% c("Gyeomdong"))
Dog.BIG.4.2 <- subset(Dog.BIG.4, Subset8 %in% c("Makdoong"))
Dog.BIG.4.3 <- subset(Dog.BIG.4, Subset8 %in% c("Batt"))
Dog.BIG.4.4 <- subset(Dog.BIG.4, Subset8 %in% c("Koko"))
Dog.BIG.4.5 <- subset(Dog.BIG.4, Subset8 %in% c("Sojin"))
Dog.BIG.4.6 <- subset(Dog.BIG.4, Subset8 %in% c("Bori"))

mean(Dog.BIG.3.1$nFeature_RNA) #1140.002
mean(Dog.BIG.3.2$nFeature_RNA) #1289.355
mean(Dog.BIG.3.3$nFeature_RNA) #1343.48
mean(Dog.BIG.3.4$nFeature_RNA) #1195.503
mean(Dog.BIG.3.5$nFeature_RNA) #1029.389
mean(Dog.BIG.3.6$nFeature_RNA) #1296.635
mean(Dog.BIG.3$nFeature_RNA) #1190.398

mean(Dog.BIG.4.1$nFeature_RNA) #1224.956
mean(Dog.BIG.4.2$nFeature_RNA) #1388.427
mean(Dog.BIG.4.3$nFeature_RNA) #1358.856
mean(Dog.BIG.4.4$nFeature_RNA) #1147.215
mean(Dog.BIG.4.5$nFeature_RNA) #1091.038
mean(Dog.BIG.4.6$nFeature_RNA) #1368.249
mean(Dog.BIG.4$nFeature_RNA) #1226.69

mean(Dog.BIG.3.1$percent.ens) #13.49681
mean(Dog.BIG.3.2$percent.ens) #15.37725
mean(Dog.BIG.3.3$percent.ens) #13.14333
mean(Dog.BIG.3.4$percent.ens) #11.92306
mean(Dog.BIG.3.5$percent.ens) #12.97885
mean(Dog.BIG.3.6$percent.ens) #13.21343
mean(Dog.BIG.3$percent.ens) #13.11932

mean(Dog.BIG.4.1$percent.loc) #8.039804
mean(Dog.BIG.4.2$percent.loc) #14.19685
mean(Dog.BIG.4.3$percent.loc) #10.03465
mean(Dog.BIG.4.4$percent.loc) #10.62799
mean(Dog.BIG.4.5$percent.loc) #9.949519
mean(Dog.BIG.4.6$percent.loc) #13.18889
mean(Dog.BIG.4$percent.loc) #10.4681

mean(Dog.BIG.3.1$percent.mt) #0.5230792
mean(Dog.BIG.3.2$percent.mt) #0.7024712
mean(Dog.BIG.3.3$percent.mt) #0.5569845
mean(Dog.BIG.3.4$percent.mt) #0.5094427
mean(Dog.BIG.3.5$percent.mt) #0.6532966
mean(Dog.BIG.3.6$percent.mt) #0.4080037
mean(Dog.BIG.3$percent.mt) #0.5598473

mean(Dog.BIG.4.1$percent.mt) #0.000275252
mean(Dog.BIG.4.2$percent.mt) #0.0002061768
mean(Dog.BIG.4.3$percent.mt) #0.0001526362
mean(Dog.BIG.4.4$percent.mt) #1.724343e
mean(Dog.BIG.4.5$percent.mt) #0.0001223809
mean(Dog.BIG.4.6$percent.mt) #0.0001812653
mean(Dog.BIG.4$percent.mt) #0.0001422935

meta.3 <- Dog.BIG.3@meta.data
meta.4 <- Dog.BIG.4@meta.data
write.table(meta.3,"meta.4.csv", append = "F", sep = ",")

VlnPlot(Dog.BIG.3.1, features = "nFeature_RNA", cols = "#FF4B20") + NoLegend()
ggsave("PATH_TO_FILE", width = 3, height = 4)
VlnPlot(Dog.BIG.4.1, features = "nFeature_RNA", cols = "#FF4B20") + NoLegend()
ggsave("PATH_TO_FILE", width = 3, height = 4)

#
ES.immune <-readRDS("ES.immune.singlet.automatic_annotation.rds")

monaco.CD4 <- readRDS("Data_Healthy/Annotation/monaco.Dog.combined.singlet.CD4.rds")
hpca.CD4 <- readRDS("Data_Healthy/Annotation/hpca.Dog.combined.singlet.CD4.rds")
nover.CD4 <- readRDS("Data_Healthy/Annotation/nover.Dog.combined.singlet.CD4.rds")
blue.CD4 <-readRDS("Data_Healthy/Annotation/blue.Dog.combined.singlet.CD4.rds")
immune.CD4 <-readRDS("Data_Healthy/Annotation/immune.Dog.combined.singlet.CD4.rds")
immgen.CD4 <-readRDS("Data_Healthy/Annotation/immgen.Dog.combined.singlet.CD4.rds")

monaco.B <- readRDS("Data_Healthy/Annotation/monaco.Dog.combined.singlet.B.rds")
hpca.B <- readRDS("Data_Healthy/Annotation/hpca.Dog.combined.singlet.B.rds")
nover.B <- readRDS("Data_Healthy/Annotation/nover.Dog.combined.singlet.B.rds")
blue.B <-readRDS("Data_Healthy/Annotation/blue.Dog.combined.singlet.B.rds")
immune.B <-readRDS("Data_Healthy/Annotation/immune.Dog.combined.singlet.B.rds")
immgen.B <- readRDS("Data_Healthy/Annotation/immgen.Dog.combined.singlet.B.rds")

monaco.myeloid <- readRDS("Data_Healthy/Annotation/monaco.Dog.combined.singlet.myeloid.rds")
hpca.myeloid <- readRDS("Data_Healthy/Annotation/hpca.Dog.comidined.singlet.myeloid.rds")
nover.myeloid <- readRDS("Data_Healthy/Annotation/nover.Dog.combined.singlet.myeloid.rds")
blue.myeloid <-readRDS("Data_Healthy/Annotation/blue.Dog.combined.singlet.myeloid.rds")
immune.myeloid <-readRDS("Data_Healthy/Annotation/immune.Dog.combined.singlet.myeloid.rds")
immgen.myeloid <-readRDS("Data_Healthy/Annotation/immgen.Dog.combined.singlet.myeloid.rds")

monaco.CD8NKT <- readRDS("Data_Healthy/Annotation/monaco.Dog.combined.singlet.CD8NKT.rds")
hpca.CD8NKT <- readRDS("Data_Healthy/Annotation/hpca.Dog.combined.singlet.CD8NKT.rds")
nover.CD8NKT <- readRDS("Data_Healthy/Annotation/nover.Dog.combined.singlet.CD8NKT.rds")
blue.CD8NKT <-readRDS("Data_Healthy/Annotation/blue.Dog.combined.singlet.CD8NKT.rds")
immune.CD8NKT <-readRDS("Data_Healthy/Annotation/immune.Dog.combined.singlet.CD8NKT.rds")
immgen.CD8NKT <-readRDS("Data_Healthy/Annotation/immgen.Dog.combined.singlet.CD8NKT.rds")

#
FeaturePlot(Dog.combined.singlet, label = T,
            features = c("PIK3CA"), reduction = "tsne")

VlnPlot(Dog.combined, features = c("percent.ens"), pt.size = 0, 
        group.by = "orig.ident") + 
  NoLegend() # + stat_summary(fun.y=median, geom="point", size=2, color="black")
ggsave("percent.ens.png", width = 3, height = 3)

#
Dog.combined <- JackStraw(Dog.combined, num.replicate = 100)
Dog.combined <- ScoreJackStraw(Dog.combined, dims = 1:20)
JackStrawPlot(Dog.combined, reduction = "pca", dims = 1:20, ymax = 1, xmax = 1) +NoLegend()
ggsave("PATH_TO_FILE", width = 3, height = 4)

ElbowPlot(Dog.combined, ndims = 30)
#
Dog.combined <- FindVariableFeatures(Dog.combined, selection.method = "vst", nfeatures = 3000)
write.table(VariableFeatures(Dog.combined), "Data_Healthy/variable.features.csv", append = F, sep = ",")
Dog.combined.3 <- FindVariableFeatures(Dog.combined.3, selection.method = "vst", nfeatures = 3000)
write.table(VariableFeatures(Dog.combined.3), "Data_Healthy/variable.features.canfam3.1.csv", append = F, sep = ",")

canfam3.1 <- VariableFeatures(Dog.combined.3)
canfam4 <- VariableFeatures(Dog.combined)

canfam4 <- as.data.frame(canfam4)    
canfam3.1 <- as.data.frame(canfam3.1)    

union(canfam4$canfam4, canfam3.1$canfam3.1)
intersect(canfam4$canfam4, canfam3.1$canfam4)
deg.canfam4 <- setdiff(canfam4$canfam4, canfam3.1$canfam3.1)
deg.canfam3 <- setdiff(canfam3.1$canfam3.1, canfam4$canfam4)
deg.union <- union(canfam4$canfam4, canfam3.1$canfam3.1)

write.table(deg.canfam4, "Data_Healthy/variable.features.canfam4specific.csv", append = F, sep = ",")
write.table(deg.canfam3, "Data_Healthy/variable.features.canfam3specific.csv", append = F, sep = ",")
write.table(deg.union, "Data_Healthy/variable.features.canfam4_3_common.csv", append = F, sep = ",")

genes <- read.csv("Data_Healthy/variable.features.canfam4specific.csv")
genes <- read.csv("Data_Healthy/variable.features.canfam4specific_noloc.csv")
genes <- read.csv("Data_Healthy/CanFam4_specific_loc.csv")
genes <- read.csv("Data_Healthy/Tumorescape.csv")

genes <- read.csv("Data_Healthy/variable.features.canfam4.csv")

counts <- GetAssayData(Dog.combined.singlet, assay = "RNA", layer = "data")
counts <- t(as.matrix(counts[rownames(counts) %in% genes$Genes,])) 
meta <- Dog.combined.singlet[[]][,c("Subset", "seurat_clusters")]

counts <- data.frame(meta, counts)
ncol(counts)
View(counts)
heatmap <- counts %>%
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:875, mean)) # 1490, 
heatmap <- as.data.frame(heatmap)

write.table(heatmap, "Data_Healthy/CanFam4_specific_hvg_615.csv", append = F, sep = ",")
#need to make rownames to match the heatmap to the annotation
nrow(heatmap)
View(heatmap)
rownames(heatmap) <- paste0("C", 0:50) # rownames will be generated as X1, X2, ..., X35
headers <- heatmap[,1:2]
heatmap <- heatmap[,3:877] # All T

#Defining the color scheme
SubsetColors <-  c("#FF4B20","blue")
names(SubsetColors) <- c("Healthy","Healthy")
clusterColors <- c(scales::hue_pal()(51)) 
names(clusterColors) <- 0:50
colors <- list(Subset = SubsetColors, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:50)

#Plotting non-scaled version
pdf("Data_Healthy/CanFam4_Specific_Variable_1490_row_height.pdf", width = 12, height = 10)
pdf("Data_Healthy/CanFam4_Specific_Variable_NoLoc_615_row_height.pdf", width = 12, height = 80)
pdf("Data_Healthy/CanFam4_Specific_Variable_Loc_875_row_height.pdf", width = 12, height = 99)
pheatmap::pheatmap(t(heatmap2),show_colnames = T, scale = "row",
                   annotation_col = headers, cluster_cols= T, cluster_rows = T,  
                   color = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
                   annotation_colors = colors, fontsize = 12, 
                   #  cellwidth = 10, 
                   # cellheight = 12,
                   legend = T, legend_labels = T, annotation_legend = T)
dev.off()

library(RColorBrewer)

#

top10 <- head(VariableFeatures(Dog.combined), 10)
VariableFeaturePlot(Dog.combined) ->plot1
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0) #+ NoLegend()
ggsave("PATH_TO_FILE", width = 4, height = 4)

#
# setwd("PATH_TO_PROJECT")
DimPlot(Dog.combined, group.by = "Subset6", shuffle = T, reduction = "tsne",
        cols = brewer.pal(3, "Set1")
) #+ NoLegend() #
DimPlot(Dog.combined, group.by = "Subset6", shuffle = T, reduction = "tsne",
        cols = c("#FF4B20", "#0348A6", "#7AC5FF")) + NoLegend()
ggsave("Tsne_Breed.png", width = 4, height = 4)

DimPlot(Dog.combined, split.by = "orig.ident", shuffle = T, reduction = "tsne",
        label = F,repel = T, #label.size = 3,
        # cols = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6","#7E6097")
) + NoLegend()
ggsave("PATH_TO_FILE", width = 7, height = 3)

DimPlot(Dog.combined, reduction = "tsne", label = F) +NoLegend()
ggsave("Tsne_Before_Singlet.png", width = 5, height = 4)

# change the resolution 
DefaultAssay(Dog.combined) <- "integrated"
Dog.combined <- FindClusters(Dog.combined, resolution = 3.7) #2.5, 3.0, 3.1, 3.7
DefaultAssay(Dog.combined) <- "RNA"
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T)
DimPlot(Dog.combined.singlet, reduction = "umap", label = T)
FeaturePlot(Dog.combined, features = c("CD79B","CD68","CD3E"))

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) + NoLegend()
ggsave("tSNE_singlet.png", width = 5, height = 5)

sample(1:30)
#
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters
num.clusters <- length(unique(Idents(Dog.combined.singlet)))
color_pal <- scales::hue_pal()(num.clusters)
color_pal <- as.vector(color_pal)

#
FeaturePlot(Dog.combined.singlet, pt.size = 2,
            features = c("Dog_MMVD"), 
            cols = viridis(100), 
            order = T, reduction = "tsne") + NoLegend()
ggsave("Data_Healthy/tSNE_Dog_MMVD_Feature.png", width = 5, height = 5)

#canfam4
DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # Plain
        cols =c("grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey")) + NoLegend()

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T, # Plain
        cols =c("#F8766D", 
                "#F37B59", "#EE8042", "#E7851E", "#E08A00", "#D98F00", "#D09400", "#C79800", "#BD9D00", "#B2A100","#A5A500",  #1-10
                "#98A800", "#89AC00", "#77AF00", "#62B200", "#45B500", "#00B709", "#00BA38", "#00BC51", "#00BD65","#00BF77",  #11-20
                "#00C087", "#00C096", "#00C1A4", "#00C0B2", "#00C0BE", "#00BECA", "#00BCD6", "#00BAE0", "#00B7E9","#00B3F2",  #21-30
                "#00AEF9", "#00A9FF", "#28A3FF", "#619CFF", "#8295FF", "#9C8DFF", "#B186FF", "#C37EFF", "#D277FF","#DF70F9",  #31-40
                "#E96AF1", "#F166E8", "#F863DE", "#FC61D3", "#FF61C7", "#FF62BA", "#FF65AD", "#FF689E", "#FF6D8F","#FC717F")) #41-50

DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # Miscellaneous
        cols =c("grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "#00C0BE", "grey", "#00BCD6", "grey", "grey","grey",
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "#D277FF","grey", 
                "grey", "grey", "#F863DE", "grey", "grey", "grey", "grey", "#FF689E", "grey","grey")) + NoLegend()
ggsave("PATH_TO_FILE", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # B
        cols =c("grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "#77AF00", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "#00BECA", "grey", "grey", "grey","grey",
                "grey", "grey", "grey", "grey", "grey", "grey", "#B186FF", "grey", "grey","grey", 
                "#E96AF1", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","#FC717F")) + NoLegend()
ggsave("tSNE_B.png", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # CD8
        cols =c("grey", 
                "#F37B59", "grey", "grey", "grey", "grey", "#D09400", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "#00B709", "grey", "grey", "grey","grey", 
                "grey", "#00C096", "grey", "#00C0B2", "grey", "grey", "grey", "grey", "grey","grey",
                "grey", "#00A9FF", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey")) + NoLegend()
ggsave("tSNE_CD8.png", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # Myeloid
        cols =c("grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "#C79800", "grey", "#B2A100","#A5A500", 
                "#98A800", "grey", "grey", "#62B200", "grey", "grey", "#00BA38", "#00BC51", "#00BD65","#00BF77", 
                "grey", "grey", "#00C1A4", "grey", "grey", "grey", "grey", "grey", "grey","#00B3F2",
                "grey", "grey", "#28A3FF", "#619CFF", "#8295FF", "#9C8DFF", "grey", "#C37EFF", "grey","grey", 
                "grey", "#F166E8", "grey", "#FC61D3", "#FF61C7", "#FF62BA", "#FF65AD", "grey", "#FF6D8F","grey")) + NoLegend()
ggsave("tSNE_Myeloid.png", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = F, # CD4
        cols =c("#F8766D", 
                "grey", "#EE8042", "#E7851E", "#E08A00", "#D98F00", "grey", "grey", "#BD9D00", "grey","grey", 
                "grey", "#89AC00", "grey", "grey", "#45B500", "grey", "grey", "grey", "grey","grey", 
                "#00C087", "grey", "grey", "grey", "grey", "grey", "grey", "#00BAE0", "#00B7E9","grey",
                "#00AEF9", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", 
                "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey")) + NoLegend()
ggsave("tSNE_CD4.png", width = 4, height = 4)



#Canfam3.1
DimPlot(Dog.combined.singlet, reduction = "tsne",label = T,
        cols =
          c("#F8766D", 
            "#F37C58", "#ED813E", "#E68613", "#DE8C00", "#D69100", "#CD9600", "#C29A00", "#B79F00", "#ABA300", "#9DA700", 
            "#8EAB00", "#7CAE00", "#66B200", "#49B500", "#0BB702", "#00BA38", "#00BC52", "#00BE67", "#00BF7A", "#00C08B", 
            "#00C19A", "#00C1A9", "#00C0B7", "#00BFC4", "#00BDD1", "#00BBDC", "#00B8E7", "#00B4F0", "#00AFF8", "#00A9FF", 
            "#22A3FF", "#619CFF", "#8494FF", "#9F8CFF", "#B584FF", "#C77CFF", "#D674FD", "#E36EF6", "#ED68ED", "#F564E3",
            "#FB61D8", "#FF61CC", "#FF62BF", "#FF64B0", "#FF68A1", "#FF6C91", "#FC7180"))

# CD4 
DimPlot(Dog.combined.singlet.3, reduction = "tsne",label = F,
        cols =
          c("#F8766D", 
            "grey", "#ED813E", "#E68613", "grey", "grey", "#CD9600", "grey", "grey", "grey","#9DA700",
            "#8EAB00", "#7CAE00", "#66B200", "grey", "grey", "grey", "grey", "grey", "#00BF7A","grey",
            "grey", "grey", "grey", "#00BFC4", "#00BDD1", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "#D674FD", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey")) + NoLegend()
ggsave("Data_Healthy/tsne_CD4_frame_corrected.png", width = 3, height = 3)
dev.off()

# CD8
DimPlot(Dog.combined.singlet.3, reduction = "tsne",label = F,
        cols =
          c("grey", 
            "#F37C58", "grey", "grey", "grey", "#D69100", "grey", "#C29A00", "grey", "#B79F00","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "#00C19A", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey")) + NoLegend()
ggsave("Data_Healthy/tsne_CD8_frame_corrected.png", width = 3, height = 3)

# Myeloid
DimPlot(Dog.combined.singlet, reduction = "tsne",label = F,
        cols =
          c("grey", 
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "#B79F00", "grey","grey",
            "grey", "grey", "grey", "#49B500", "#0BB702", "#00BA38", "#00BC52", "#00BE67", "grey","#00C08B",
            "grey", "grey", "#00C0B7", "grey", "grey", "#00BBDC", "#00B8E7", "grey", "#00AFF8","#00A9FF",
            "grey", "#619CFF", "grey", "#9F8CFF", "grey", "grey", "grey", "grey", "grey","grey",
            "#FB61D8", "#FF61CC", "#FF62BF", "grey", "grey", "#FF6C91", "grey")) + NoLegend()

DimPlot(Dog.combined.singlet, reduction = "tsne",label = T,
        cols =
          c("#F8766D", 
            "#F37C58", "#ED813E", "#E68613", "#DE8C00", "#D69100", "#CD9600", "#C29A00", "#B79F00", "#ABA300", "#9DA700", 
            "#8EAB00", "#7CAE00", "#66B200", "#49B500", "#0BB702", "#00BA38", "#00BC52", "#00BE67", "#00BF7A", "#00C08B", 
            "#00C19A", "#00C1A9", "#00C0B7", "#00BFC4", "#00BDD1", "#00BBDC", "#00B8E7", "#00B4F0", "#00AFF8", "#00A9FF", 
            "#22A3FF", "#619CFF", "#8494FF", "#9F8CFF", "#B584FF", "#C77CFF", "#D674FD", "#E36EF6", "#ED68ED", "#F564E3",
            "#FB61D8", "#FF61CC", "#FF62BF", "#FF64B0", "#FF68A1", "#FF6C91", "#FC7180"))

ggsave("tsne_Myeloid_frame.png", width = 3, height = 3)

# B and plasma 
DimPlot(Dog.combined.singlet, reduction = "tsne",label = F,
        cols =
          c("grey", 
            "grey", "grey", "grey", "#DE8C00", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "#8494FF", "grey", "#B584FF", "#C77CFF", "grey", "grey", "grey","#F564E3",
            "grey", "grey", "grey", "grey", "grey", "grey", "#FC7180")) + NoLegend()

DimPlot(Dog.combined.singlet, reduction = "tsne",label = T,
        cols =
          c("#F8766D", 
            "#F37C58", "#ED813E", "#E68613", "#DE8C00", "#D69100", "#CD9600", "#C29A00", "#B79F00", "#ABA300", "#9DA700", 
            "#8EAB00", "#7CAE00", "#66B200", "#49B500", "#0BB702", "#00BA38", "#00BC52", "#00BE67", "#00BF7A", "#00C08B", 
            "#00C19A", "#00C1A9", "#00C0B7", "#00BFC4", "#00BDD1", "#00BBDC", "#00B8E7", "#00B4F0", "#00AFF8", "#00A9FF", 
            "#22A3FF", "#619CFF", "#8494FF", "#9F8CFF", "#B584FF", "#C77CFF", "#D674FD", "#E36EF6", "#ED68ED", "#F564E3",
            "#FB61D8", "#FF61CC", "#FF62BF", "#FF64B0", "#FF68A1", "#FF6C91", "#FC7180"))
ggsave("tsne_B_frame.png", width = 3, height = 3)

#ETC
DimPlot(Dog.combined.singlet.3, reduction = "tsne",label = F,
        cols =
          c("grey", 
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey",
            "grey", "#00C1A9", "grey", "grey", "grey", "grey", "grey", "#00B4F0", "grey","grey",
            "#22A3FF", "grey", "grey", "grey", "grey", "grey", "#D674FD", "grey", "#ED68ED","grey",
            "grey", "grey", "grey", "#FF64B0", "#FF68A1", "grey", "grey")) + NoLegend()
ggsave("Data_Healthy/tsne_ETC_frame_corrected.png", width = 3, height = 3)

DimPlot(Dog.combined.singlet, reduction = "tsne",label = T,
        cols =
          c("#F8766D", 
            "#F37C58", "#ED813E", "#E68613", "#DE8C00", "#D69100", "#CD9600", "#C29A00", "#B79F00", "#ABA300", "#9DA700", 
            "#8EAB00", "#7CAE00", "#66B200", "#49B500", "#0BB702", "#00BA38", "#00BC52", "#00BE67", "#00BF7A", "#00C08B", 
            "#00C08B", "#00C1A9", "#00C0B7", "#00BFC4", "#00BDD1", "#00BBDC", "#00B8E7", "#00B4F0", "#00AFF8", "#00A9FF",
            "#22A3FF", "#619CFF", "#8494FF", "#9F8CFF", "#B584FF", "#C77CFF", "#D674FD", "#E36EF6", "#ED68ED", "#F564E3",
            "#FB61D8", "#FF61CC", "#FF62BF", "#FF64B0", "#FF68A1", "#FF6C91", "#FC7180"))
colorcode <- show_col(hue_pal()(51))
colorcode <- show_col(hue_pal()(length(unique(Idents(Dog.combined.singlet)))))

colorcode <- c("#F8766D", 
               "#F37B59", "#EE8042", "#E7851E", "#E08A00", "#D98F00", "#D09400", "#C79800", "#BD9D00", "#B2A100","#A5A500",  #1-10
               "#98A800", "#89AC00", "#77AF00", "#62B200", "#45B500", "#00B709", "#00BA38", "#00BC51", "#00BD65","#00BF77",  #11-20
               "#00C087", "#00C096", "#00C1A4", "#00C0B2", "#00C0BE", "#00BECA", "#00BCD6", "#00BAE0", "#00B7E9","#00B3F2",  #21-30
               "#00AEF9", "#00A9FF", "#28A3FF", "#619CFF", "#8295FF", "#9C8DFF", "#B186FF", "#C37EFF", "#D277FF","#DF70F9",  #31-40
               "#E96AF1", "#F166E8", "#F863DE", "#FC61D3", "#FF61C7", "#FF62BA", "#FF65AD", "#FF689E", "#FF6D8F","#FC717F")
#
# setwd("PATH_TO_PROJECT")

#
Dog.combined@active.ident <- Dog.combined$seurat_clusters
DimPlot(Dog.combined, reduction = "tsne", label = T, label.size = 3) +NoLegend()
ggsave("Data_Healthy/Dimplot_tsne.png", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) + NoLegend()
ggsave("Data_Healthy/Dimplot_tsne.png", width = 5, height = 5)

# Annotation cluster performance
FeaturePlot(Dog.combined.singlet, features = c("S100A4"), 
            reduction = "tsne", label = T, 
            #pt.size = 0.8, 
            order = F,repel = T,
            cols = c("#ADDAE6","#E63222") #E63222
            #c("#ADDAE6","#E63222"), #cols = c("lightgrey", "darkred"), #cols = c("lightgrey", "darkgreen")
) + NoLegend()
ggsave("Featureplot_Cluster_Marker/CD8A.png",
       width = 4, height = 4)

cluster3948.markers <- FindMarkers(Dog.combined.singlet, ident.1 = c(39,48), logfc.threshold = 1,only.pos = TRUE, min.pct = 0.5)
head(cluster3948.markers, 10)

dotPlot(Dog.combined.singlet, 
        features = c("PDCD1","CTLA4","HAVCR2","TIGIT","LAG3","TOX","CD274","TCF7"),
        #angle.y = 180, angle.x = 90,
        idents = c(9,2,12,10,21,38,0,11,3,19,25,13,6,24,37,39,31,44,45,38))
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Dotplot_T_exhaustion.png", 
       width = 5.75, height = 2.5)

plot_density(Dog.combined.singlet, 
             features = c("CAECAM1"), size = 0.5, pal = "plasma",# "cividis", ""plasma", "magma"
             reduction = "tsne") +NoLegend()
ggsave("Featureplot_Cluster_Marker/Myeloid_PLD4.png",
       width = 4, height = 4)
DotPlot(Dog.combined.singlet, features = c("LOC607937"))
VlnPlot(Dog.combined.singlet, features = c("CD247"))
# Do the heatmap to check the global expression of genes defining each cluster
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE, features = rownames(Dog.combined.singlet)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet <- JoinLayers(Dog.combined.singlet)
Dog.markers <- FindAllMarkers(Dog.combined.singlet, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top10, "Cluster_Marker/Top10_0.36_0.5.csv",append = "F",sep = ",")
ggsave("Data_Healthy/Tsne_Labeled_Screening_3.5.png", width=7, height=5)

DoHeatmap(subset(Dog.combined.singlet, downsample = 200),
          features = top10$gene, 
          #disp.max = 2,
          disp.min = -2,
          size = 3, angle = 50, cells = 1:200)+
  theme(text = element_text(size=8)) + 
  scale_fill_viridis(option = "E") + 
  labs(x = "Viridis E", y = NULL) #+ NoLegend()

#+ scale_fill_viridis(100)
#+ NoLegend()
ggsave("Data_Healthy/Heatmap_Labeled_Singlet_legend.png", width=8, height=7)

#markers <- subset(Dog.markers, cluster == 24)
#class(genes)
#genes <- markers$gene

#
Dog.combined <-ScaleData(Dog.combined, verbose = FALSE, features = rownames(Dog.combined)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined <- JoinLayers(Dog.combined, layers = "RNA")

DimPlot(Dog.combined, reduction = "tsne")
Dog.combined.singlet@active.ident <- factor(Dog.combined.singlet@active.ident, 
                                            levels = c("T_CD4", "T_CD8","B","Myeloid","NK","T_GD","Un"))
Dog.combined@active.ident <- factor(Dog.combined@active.ident, 
                                    levels = c("T_CD4", "T_CD8","T_DN1","T_DN2","T_DN3","Tgd",
                                               "Myeloid_1","Myeloid_2","Myeloid_3","Myeloid_4",
                                               "DC","PMN", "B","B_Plasma"))


colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

# DC
dotPlot(Dog.combined.singlet, features = c("FCER1A","CD1C","CLEC1B","GCSAM","OSCAR","TCF4",
                                           "IGFBP7","DLA-DRA",
                                           "IL3RA","FLT3","SCAMP5"))
FeaturePlot(Dog.combined, features = c("CD1C","CLEC10A","FCER1A","TCF4"), 
            reduction = "tsne", label = F, label.size = 3)
DimPlot(Dog.combined, reduction = "tsne", label = T) + NoLegend()
#
Dog.combined <- RenameIdents(Dog.combined,
                             '0' = "T", '1' = "T", '2' = "T",'3' = "T", '4' = "T", 
                             '5' = "Myeloid", '6' = "T", '7' = "T", '8' = "T", '9'= "B", '10' = "T",
                             '11'= "Myeloid", '12' = "Myeloid", '13' = "Myeloid", '14' = "T", '15' = "Myeloid",
                             '16' = "Myeloid", '17' = "T", 
                             '18' = "Myeloid", '19' = "Myeloid", '20' = "T", '21' = "T",
                             '22' = "T", '23' = "T", '24' = "T", '25' = "Myeloid", '26' = "T", 
                             '27' = "T", '28' = "T",
                             '29' = "B", '30' = "Plasma", '31' = "DC", '32' = "T", '33' = "Plasma", 
                             '33' = "Plasma", '34' = "Myeloid",
                             '35' = "Myeloid", '36' = "DC", '37' = "Tgd", '38' = "pDC", '39' = "B", '40' = "DC")

Dog.combined <- RenameIdents(Dog.combined, # more specific
                             '0' = "T_CD4", '1' = "T_CD4", '2' = "T_CD8",'3' = "T_CD4", '4' = "T_CD4", 
                             '5' = "Myeloid", '6' = "T_CD4", '7' = "T_CD4", '8' = "T_CD4", '9'= "B", '10' = "T_CD8",
                             '11'= "Myeloid_1", '12' = "Myeloid", '13' = "Myeloid", '14' = "T_CD8", '15' = "Myeloid",
                             '16' = "Myeloid", '17' = "T_CD4", 
                             '18' = "Myeloid", '19' = "Myeloid", '20' = "Myeloid", '21' = "T_DN",
                             '22' = "T_CD4", '23' = "Myeloid", '24' = "T_DN", '25' = "Myeloid", '26' = "T_DN", 
                             '27' = "Myeloid", '28' = "T_DN",
                             '29' = "B", '30' = "B", '31' = "Myeloid", '32' = "T_DN", '33' = "B_Plasma", 
                             '33' = "B_Plasma", '34' = "Myeloid", '35' = "Myeloid", 
                             '36' = "Myeloid", '37' = "Tgd", '38' = "Myeloid", '39' = "B", '40' = "Myeloid")
Dog.combined@active.ident <- Dog.combined$seurat_clusters

#
Dog.combined <- RenameIdents(Dog.combined, # more specific
                             '0' = "T_CD4", '1' = "T_CD4", '2' = "T_CD8",'3' = "T_CD4", '4' = "T_CD4", 
                             '5' = "Myeloid_1", '6' = "T_CD4", '7' = "T_CD4", '8' = "T_CD4", '9'= "B", '10' = "T_CD8",
                             '11'= "Myeloid_4", '12' = "Myeloid_1", '13' = "Myeloid_2", '14' = "T_CD8", '15' = "Myeloid_1",
                             '16' = "Myeloid_2", '17' = "T_CD4", 
                             '18' = "Myeloid_2", '19' = "Myeloid_2", '20' = "Myeloid_1", '21' = "T_DN1",
                             '22' = "T_CD4", '23' = "Myeloid_1", '24' = "T_DN2", '25' = "Myeloid_2", '26' = "T_CD4", 
                             '27' = "Myeloid_1", '28' = "T_DN3",
                             '29' = "B", '30' = "B_Plasma", '31' = "DC", '32' = "T_CD4",
                             '33' = "B_Plasma", '34' = "PMN", '35' = "Myeloid_3", 
                             '36' = "DC", '37' = "Tgd", '38' = "DC", '39' = "B", '40' = "Myeloid_3")
Dog.combined@active.ident <- Dog.combined$seurat_clusters

#
Dog.combined <- RenameIdents(Dog.combined, # more specific
                             '0' = "T_CD4", '1' = "T_CD8", '2' = "T_CD4",'3' = "T_CD4", '4' = "T_CD8", 
                             '5' = "B", '6' = "T_CD4", '7' = "T_CD8", '8' = "T_CD4", '9'= "T_CD4", '10' = "T_CD4",
                             '11'= "Myeloid_1", '12' = "T_CD4", '13' = "Myeloid_2", '14' = "T_CD4", 
                             '15' = "Myeloid_1",
                             '16' = "T_CD4", '17' = "Myeloid_1", 
                             '18' = "Myeloid_4", '19' = "Myeloid_2", '20' = "Myeloid_2", '21' = "Myeloid_1",
                             '22' = "Myeloid_2", '23' = "T_CD4", '24' = "Myeloid_3", '25' = "T_CD4", '26' = "T_CD4", 
                             '27' = "Myeloid_2", '28' = "T_CD4",
                             '29' = "Myeloid_4", '30' = "Mix_Mito", '31' = "Myeloid_1", '32' = "Myeloid_2", '33' = "T_Un1", 
                             '34' = "Myeloid_3", '35' = "Myeloid_1", 
                             '36' = "T_Replicate", '37' = "B", '38' = "B_Plasma", '39' = "T_CD8", '40' = "T_DN",
                             '41' = "B_Plasma", '42' = "DC", '43' = "T_Replicate", '44' = "T_Un2", '45' = "B",
                             '46' = "T_Replicate", '47' = "Eos", '48' = "DC",
                             '49' = "T_gd", '50' = "Myeloid_4", '51' = "PMN",
                             '52' = "DC_Metabolism", '53' = "T_IFN", '54' = "B_Plasma", '55' = "T_Replicate", '56'= "Un_6")

Dog.combined.singlet.B <- RenameIdents(Dog.combined.singlet.B, # more specific
                                       '0' = "C0", '1' = "C1", '2' = "C2",'3' = "C3", '4' = "C4", 
                                       '5' = "C5", '6' = "C6", '7' = "C7", '8' = "C8", '9'= "C9")
DimPlot(Dog.combined.singlet.B, reduction = "umap",label = T)
Dog.combined.singlet.B$cell_type <- Dog.combined.singlet.B@active.ident

unique(Dog.combined.singlet.B$cell_type)
DimPlot(Dog.combined, reduction = "umap",label = T)
DimPlot(Dog.combined, reduction = "tsne",label = T)
FeaturePlot(Dog.combined,features = c("CD96"), reduction = "umap")
VlnPlot(Dog.combined, features = c("CD96"))
length(unique(Dog.combined@active.ident))
ggsave("Data_Healthy/resolution_3.7_tsne.png", width = 10, height = 7)
#
library(celldex)
library(SingleR)

#
getOption('timeout')
options(timeout=99999)

#
Hpca <- celldex::HumanPrimaryCellAtlasData() # HumanPrimaryCellAtlasData() 
Monaco <- celldex::MonacoImmuneData() # MonacoImmuneData()
Dice <- celldex::DatabaseImmuneCellExpressionData() #DatabaseImmuneCellExpressionData()
Immgen <- celldex::ImmGenData() # ImmGenData()
Blue <- celldex::BlueprintEncodeData() #BlueprintEncodeData()
Nover <- celldex::NovershternHematopoieticData()

library(SingleR)

#
hpca <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Hpca, labels = Hpca$label.fine, assay.type.test=1)
monaco <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Monaco, labels = Monaco$label.fine, assay.type.test=1)
immune <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Dice, labels = Dice$label.fine, assay.type.test=1)
immgen <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Immgen, labels = Immgen$label.fine, assay.type.test=1)
blue <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Blue, labels = Blue$label.fine, assay.type.test=1)
nover <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet), ref = Nover, labels = Nover$label.fine, assay.type.test=1)

saveRDS(immgen, "PATH_TO_FILE")
saveRDS(hpca, "PATH_TO_FILE")
saveRDS(monaco, "PATH_TO_FILE")
saveRDS(immune, "PATH_TO_FILE")
saveRDS(blue, "PATH_TO_FILE")
saveRDS(nover, "PATH_TO_FILE")

monaco <- readRDS("PATH_TO_FILE")
hpca <- readRDS("PATH_TO_FILE")
nover <- readRDS("PATH_TO_FILE")
blue <-readRDS("PATH_TO_FILE")
immune <-readRDS("PATH_TO_FILE")

Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = hpca$labels, col.name = "hpca")
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = monaco$labels, col.name = "monaco")
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = immune$labels, col.name = "immune")
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = immgen$labels, col.name = "immgen")
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = blue$labels, col.name = "blue")
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, metadata = nover$labels, col.name = "nover")

DimPlot(Dog.combined.singlet,
        reduction = "tsne", group.by  = "monaco",
        label = T, shuffle = T, repel = T, label.box = T, label.size = 6, pt.size = 1,
        cols = colorRampPalette(brewer.pal(12, "Paired"))(29))  + NoLegend()
ggsave("PATH_TO_FILE", width = 12, height = 12)


meta <- Dog.combined.singlet[[]]
meta <- meta[,c(5,16,27:55)]
ncol(meta)
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:29, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 0:47) #rownames(heatmap) <- paste0("X", 1:22)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:31)]


clusterColors <- c(scales::hue_pal()(48))
names(clusterColors) <- 0:47
SubsetColors2 <- c("#FF4B20") 
names(SubsetColors2) <- c("PBMC")


colors <- list(
  Subset2 = SubsetColors2, 
  seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 1:48) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("GSEA_Healthy.pdf", width = 15, height =20) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "none", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   #  cellwidth = 10, 
                   #   cellheight = 12,
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   ))
dev.off()




#
DimPlot(subset(Dog.combined, seurat_clusters ==35), group.by = "monaco", reduction = "tsne", label = T,repel = T)
ggsave("Data_Healthy/Cluster28_blue.png", width = 15, height = 15)

#
FeaturePlot(subset(Dog.combined, seurat_clusters ==24), order = T,
            features = c("CD4","CD8A"), 
            reduction = "tsne", label = T,repel = T)
ggsave("Data_Healthy/Cluster24_immune.png", width = 15, height = 15)

#
cluster <- Dog.combined[[]]
cluster24 <- subset(cluster, seurat_clusters == 24)
unique(cluster24$nover)
#
DimPlot(Dog.combined, reduction = "tsne", label = T) + NoLegend()
FeaturePlot(Dog.combined, features = c("PTPRC"), 
            reduction = "tsne", cols = viridis(100), 
            label = F, label.size = 3, repel = T, order = T)
ggsave("PTPRC.png", width = 4, height = 4)
#
Dog.combined@active.ident <- Dog.combined$seurat_clusters

#
FeaturePlot(Dog.combined.singlet.CD4, features = c("NKG7","PI3"), reduction = "umap", label = T)
DimPlot(Dog.combined, reduction = "tsne", label = T)
#
dotPlot(Dog.combined.singlet, features = c("PTPRC","CD3D","CD3E","CD5","CD4","CD8A","SCART1",
                                           "CSF1R","CSF3R","ENSCAFG00000005430","ENSCAFG00000029470","CD177",
                                           "CD79B","CD19","MZB1","JCHAIN","GZMB","GZMA","PRF1",
                                           "SELL","TCF7","CCR7","CD44",
                                           "CD27","CXCR4","SDC1",
                                           "PASK","ITGA4","ITGB1",
                                           "IL2RA","FOXP3","IL7R","FUT7",
                                           "CD38",
                                           "ITGAM","TNFRSF13B",
                                           "ADGRE5","SLC9A9","TCF4",
                                           "CD1C","PLD4","DLA-DRA","KLF12","IL3RA","RORA","VSIR")) 
ggsave("Data_Healthy/resolution_3.7_marker.png", width = 10, height = 9)
dotPlot(Dog.combined, features = c("S100A12","FCGR1A","CD27","CXCR4","TNFRSF4","KLRB1",
                                   "TCF4","CD209","IL3RA","ITGAM","DLA-DRA","FCER1A"))
ggsave("Data_Healthy/resolution_3.7_marker2.png", width = 10, height = 9)

FeaturePlot(Dog.combined, features = "CD79B", reduction = "tsne",label = T )
Dog.combined@active.assay <- Dog.combined$orig.ident
dotPlot(Dog.combined, features = c("IL3RA","DLA-DRA","FLT3","SCAMP5","GCSAM"))

#
sign <- read.delim("PATH_TO_FILE")

full <- as.list(sign)
unique <- names(full) 
list <- list()
for (i in seq_along(unique)) {
  tmp <- full[[i]]
  tmp <- tmp[tmp != ""]
  tmp <- unique(toupper(tmp))
  tmp <- GSEABase::GeneSet(tmp, setName=paste(unique[i]))
  list[[i]] <- tmp}
list <- GSEABase::GeneSetCollection(list)

# 
gene.sets <- list(Bcells = c("MS4A1","CD79B","CD79A"), # make this first 
                  Myeloid = c("SPI1","FCER1G","CSF1R"),
                  Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"),
                  CD8nv = c("ITGA1", "LEF1", "PTGDR", "IL2RB", "ADGRG1", "NBEA"),
                  Mo = c("MT1E","MT2A","HMOX1","C1QC","C1QA","GSTP1","CTSS")
                  # looks not specific for CD8 T cells actually 
)
Dog.combined <- JoinLayers(Dog.combined) # do this second 
library(escape) #install this and run this package

ES.immune <- escape.matrix(Dog.combined.singlet,  # run this code 
                           gene.sets = list,
                           min.size = 3,
                           method = "ssGSEA")
saveRDS(ES.immune, "ES.immune.singlet.rds")

Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, ES.immune) # add metadata

FeaturePlot(Dog.combined.singlet, features = c("Mo"), 
            reduction = "tsne", cols = viridis(100) # plasma(100)
)
DimPlot(Dog.combined.singlet, reduction = "tsne", label=T) +NoLegend()
library(viridis)
addmeta
#
saveRDS(ES.immune, "PATH_TO_FILE")
Dog.combined <- JoinLayers(Dog.combined)

ES.immune <- readRDS("PATH_TO_FILE")


#DimPlot(Dog.combined, group.by = "Bcells", reduction = "tsne")

#
Dog.combined.singlet@meta.data <- Dog.combined.singlet@meta.data[,c(1:7,16)]

Dog.combined.singlet$percent.loc <- NULL


#
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, ES.immune) 
meta <- Dog.combined.singlet[[]]
meta$seurat_clusters <- Dog.combined.singlet@active.ident
ncol(meta)
View(meta)
View(heatmap)
meta <- meta[,c(5,16,21:118)] 
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:98, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 1:51) #rownames(heatmap) <- paste0("X", 1:22)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:100)]

unique(Dog.combined.singlet$Subset2)
subset <- unique(Dog.combined$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet$Subset2) # PBMC TIL
subset3 <- unique(Dog.combined$Subset3) # Healthy PBMC AL NAC MGT CRC LAC
subset4 <- unique(Dog.combined$Subset4) # 
subset5 <- unique(Dog.combined$Subset5) # 
subset6 <- unique(Dog.combined$Subset6) # 


SubsetColors <- c("#FCE540", "#3F114E") 
names(SubsetColors) <- subset2 #c("Healthy", "Tumor") 
SubsetColors2 <- c("darkred", "#0348A6") 
names(SubsetColors2) <- c("PBMC", "TIL") 
SubsetColors3 <- c("#DDDDDD", "#FFB433", "#C6FDEC", "#7AC5FF", "#7E6097","#FF4B20","darkgreen")  
names(SubsetColors3) <- subset3 #c("NA","AL","NAC","MGT","CRC","LAC") 
SubsetColors4 <- Glasbey
names(SubsetColors4) <- subset4 #c("NA","AL","NAC","MGT","CRC","LAC") 
clusterColors <- c(scales::hue_pal()(51))
names(clusterColors) <- 0:50

clusterColors <- c(scales::hue_pal()(57))
names(clusterColors) <- 0:56
SubsetColors6 <- c("#FF4B20", "#FFB433", "#7AC5FF") 
names(SubsetColors6) <- subset6 #c("Healthy", "Tumor") 

install.packages("Polychrome")
library(Polychrome)
library(ggplot2)
library(pals)
library(wesanderson)

Glasbey<-glasbey.colors(13)
Glasbey<-alphabet.colors(13)

swatch(Glasbey)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)
wes_palette("Moonrise2", 13, type = "continuous")
colors <- list(#Subset = SubsetColors, 
  Subset2 = SubsetColors, 
  # Subset3 = SubsetColors3, 
  # Subset4 = SubsetColors4, 
  # Subset5 = SubsetColors5, 
  #Subset6 = SubsetColors6, 
  
  seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 1:51) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Seurat_V5/GSEA_All_1.pdf", width = 15, height =10) #5.5 x 3.8, 9x4
pdf("GSEA_Healthy.pdf", width = 15, height =20) #5.5 x 3.8, 9x4
pdf("GSEA_Annotation_Singlet.pdf", width = 15, height =20) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   #   cellwidth = 10, 
                   #   cellheight = 10,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()
DotPlot(Dog.combined.singlet, features = "CD177")
FeaturePlot(Dog.combined, features = c("Migration_Inflammatory", "Adhesion_Arterial_Endothelium"), 
            reduction = "tsne", blend = T, order = F)
ggsave("M2_Anti_Combined.png", width = 10, height = 3)

FeaturePlot(Dog.combined, features = c("Migration_Inflammatory"), 
            reduction = "tsne")

DimPlot(Dog.combined, reduction = "tsne", label = T) +NoLegend()
DimPlot(Dog.combined, reduction = "umap", label = T) +NoLegend()

# How to find singlet and differentiated doublets 
library(Seurat)
Dog.combined <- JoinLayers(Dog.combined)
sce <- as.SingleCellExperiment(Dog.combined)  
sce <- scDblFinder(sce, clusters = , samples = ) # FIND THE DOUBLET
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined <- AddMetaData(Dog.combined, doublets)

DimPlot(Dog.combined, reduction = "tsne", group.by = "db.class")
unique(Dog.combined$db.class)
#

Dog.combined.singlet <- subset(Dog.combined, db.class %in% c("singlet"))

Dog.combined.singlet <- subset(Dog.combined, subset = db.class == c("singlet")) #Dog.combined$db.class = 
DefaultAssay(Dog.combined.singlet) <- "integrated"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <-RunPCA(Dog.combined.singlet, npcs = 30, verbose = FALSE)
Dog.combined.singlet <-RunUMAP(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- RunTSNE(Dog.combined.singlet, dims = 1:30, reduction = "pca")
Dog.combined.singlet <- FindNeighbors(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- FindClusters(Dog.combined.singlet, resolution = 3.7)
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <- JoinLayers(Dog.combined.singlet)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T)+ NoLegend()
DimPlot(Dog.combined, reduction = "tsne", label = T)+ NoLegend()

saveRDS(Dog.combined.singlet, "Dog.combined.singlet_good2.rds")

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T, 
        cols = colorRampPalette(brewer.pal(9, "Set1"))(48)) + NoLegend()
#
cols = brewer.pal(11, "Spectral")
DimPlot(Dog.combined.singlet, reduction = "umap", label = T) + NoLegend()

# doublet 
Dog.combined.doublet <- subset(Dog.combined, subset = db.class == c("doublet"))
DefaultAssay(Dog.combined.doublet) <- "integrated"
Dog.combined.doublet <-ScaleData(Dog.combined.doublet, verbose = FALSE)
Dog.combined.doublet <-RunPCA(Dog.combined.doublet, npcs = 30, verbose = FALSE)
Dog.combined.doublet <-RunUMAP(Dog.combined.doublet, reduction = "pca", dims = 1:30)
Dog.combined.doublet <- RunTSNE(Dog.combined.doublet, dims = 1:30, reduction = "pca")
Dog.combined.doublet <- FindNeighbors(Dog.combined.doublet, reduction = "pca", dims = 1:30)
Dog.combined.doublet <- FindClusters(Dog.combined.doublet, resolution = 0.6)
DefaultAssay(Dog.combined.doublet) <- "RNA"
Dog.combined.doublet <-ScaleData(Dog.combined.doublet, verbose = FALSE)
Dog.combined.doublet <- JoinLayers(Dog.combined.doublet)

saveRDS(Dog.combined.doublet, "Dog.combined.doublet.rds")

DimPlot(Dog.combined.doublet, reduction = "tsne", label = F) + NoLegend()
ggsave("Data_Healthy/Doublet_Subclustering.png", width = 5, height = 5)
DimPlot(Dog.combined.doublet, reduction = "umap", label = T) + NoLegend()

DimPlot(Dog.combined, group.by = "db.class", reduction = "tsne")
FeaturePlot(Dog.combined.doublet, order = T, reduction = "tsne",
            features = c("CD3D","CD3E","CD5","CD4","CD8A","LYZ","ITGAM",
                         "CSF3R","DPYD","S100A9","CD79B"))
FeaturePlot(Dog.combined.doublet, blend = T, features = c("CD3E","DPYD"))


#
doublet.percentage <- table(Idents(Dog.combined), Dog.combined$db.class)
write.table(doublet.percentage, "PATH_TO_FILE", sep = ",", append=F)

Cellnumbers.subset <- table(Idents(Dog.combined.singlet.3), Dog.combined.singlet.3$orig.ident)
write.table(Cellnumbers.subset,
            "Cellnumber_canfam3.orig.ident.csv", sep = ",", append=F)

Cellnumbers.subset <- table(Idents(Dog.combined.singlet.myeloid), Dog.combined.singlet.myeloid$orig.ident)
write.table(Cellnumbers.subset,
            "Cellnumber_canfam3.orig.ident.csv", sep = ",", append=F)

DimPlot(Dog.combined, reduction = "tsne", label = T) + NoLegend()
DimPlot(Dog.combined, reduction = "tsne", label = F, cols = c("grey","red"),
        group.by = "db.class") + NoLegend()
ggsave("Single_Doublet.png", width = 5, height = 5) + NoLegend()


FeaturePlot(Dog.combined.singlet.myeloid, 
            features = c("RHEX"), 
            order = T, pt.size =1.5,
            label = T, 
            label.size = 3,
            reduction = "tsne") 

DimPlot(Dog.combined.singlet, reduction = "umap", label = T) 
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) 
DimPlot(Dog.combined, reduction = "tsne", label = F, 
        group.by = "db.class", 
        cols = c("#ADDAE6","#E63222")) +NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)

plot_density(Dog.combined.singlet,# size = 0.1, label = T,
             reduction = "tsne", 
             features = c("AK8"))

#
Dog.combined.singlet$seurat_clusters
DimPlot(Dog.combined.singlet, label = T, reduction = "tsne")
DimPlot(Dog.combined.singlet, label = T, reduction = "umap")
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet, 
                                     '0' = "T_CD4", '1' = "T_CD8", '2' = "T_CD4",'3' = "T_CD4", '4' = "B", 
                                     '5' = "T_CD8", '6' = "T_CD4", '7' = "T_CD8", '8' = "Neu_1", '9'= "T_CD4", 
                                     '10' = "T_CD4", '11'= "T_CD4", '12' = "T_CD4", '13' = "T_CD4", '14' = "Neu_2", 
                                     '15' = "Mo_2", '16' = "Mo_1", '17' = "Neu_2", 
                                     '18' = "Mo_2", '19' = "T_CD4", '20' = "Mo_2", '21' = "T_CD4",
                                     '22' = "T_CD4", '23' = "Neu_1", '24' = "T_CD4", '25' = "T_CD4", '26' = "Mo_2", 
                                     '27' = "Mo_1", '28' = "Mix_Mt",
                                     '29' = "Neu_2", '30' = "Mo_2", '31' = "T_Prolif", '32' = "Neu_1", '33' = "Plasma", 
                                     '34' = "DC", '35' = "Plasma", 
                                     '36' = "B", '37' = "T_CD4", '38' = "T_CD4", '39' = "T_CD4", '40' = "B",
                                     '41' = "Eos_Baso", '42' = "DC_P", '43' = "Neu_1", '44' = "T_Prolif", '45' = "T_Gd",
                                     '46' = "Baso", '47' = "Plasma")

Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters 
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet, 
                                     '0' = "T_CD4", '1' = "T_CD8", '2' = "T_CD4",'3' = "T_CD4", '4' = "B", 
                                     '5' = "T_CD8", '6' = "T_CD4", '7' = "T_CD8", '8' = "Myeloid", '9'= "T_CD4", 
                                     '10' = "T_CD4", '11'= "T_CD4", '12' = "T_CD4", '13' = "T_CD4", '14' = "Myeloid", 
                                     '15' = "Myeloid", '16' = "Myeloid", '17' = "Myeloid", 
                                     '18' = "Myeloid", '19' = "T_CD4", '20' = "Myeloid", '21' = "T_CD4",
                                     '22' = "T_CD4", '23' = "Myeloid", '24' = "T_CD4", '25' = "T_CD4", '26' = "Myeloid", 
                                     '27' = "Myeloid", '28' = "Mix_Mt",
                                     '29' = "Myeloid", '30' = "Myeloid", '31' = "T_Prolif", '32' = "Myeloid", '33' = "B", 
                                     '34' = "Myeloid", '35' = "B", 
                                     '36' = "B", '37' = "T_CD4", '38' = "T_CD4", '39' = "T_CD4", '40' = "B",
                                     '41' = "Myeloid", '42' = "Myeloid", '43' = "Myeloid", '44' = "T_Prolif", '45' = "T_Gd",
                                     '46' = "Myeloid", '47' = "B")

#
FeaturePlot(Dog.combined.singlet, order = T,
            features = c("NK_PolyIC_Injection"),cols = viridis(100),
            reduction = "tsne", label = T)

#
Dog.combined.singlet.T <- subset(Dog.combined.singlet, idents = c(1,24,16,32,22,6,27,25,39,48,43,
                                                                  25,13,6,19,22,24,3,0,11,10,21,38,9,2,12))

DefaultAssay(Dog.combined.singlet.T) <- "integrated"
Dog.combined.singlet.T <-ScaleData(Dog.combined.singlet.T, verbose = FALSE)
Dog.combined.singlet.T <-RunPCA(Dog.combined.singlet.T, npcs = 30, verbose = FALSE)
Dog.combined.singlet.T <-RunUMAP(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <-RunTSNE(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <- FindNeighbors(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <- FindClusters(Dog.combined.singlet.T, resolution = 2)
DefaultAssay(Dog.combined.singlet.T) <- "RNA"
Dog.combined.singlet.T <-ScaleData(Dog.combined.singlet.T, verbose = FALSE)
Dog.combined.singlet.T <- JoinLayers(Dog.combined.singlet.T)
DimPlot(Dog.combined.singlet.T, reduction = "tsne", label = T)
FeaturePlot(Dog.combined.singlet.T, features = c("FOXP3","IL2RA","PDCD1"), order = T, reduction = "tsne")
saveRDS(Dog.combined.singlet.T, "PATH_TO_FILE")


#
library(Seurat)
FeaturePlot(Dog.combined.singlet, features = c("CD8A","CD8B"), reduction = "tsne", label = T)
Dog.combined.singlet.CD8NKT <- subset(Dog.combined.singlet, idents = c(1,24,16,32,22,6,27,25,39,48,43))

DefaultAssay(Dog.combined.singlet.CD8NKT) <- "integrated"
Dog.combined.singlet.CD8NKT <-ScaleData(Dog.combined.singlet.CD8NKT, verbose = FALSE)
Dog.combined.singlet.CD8NKT <-RunPCA(Dog.combined.singlet.CD8NKT, npcs = 30, verbose = FALSE)
Dog.combined.singlet.CD8NKT <-RunUMAP(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <-RunTSNE(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <- FindNeighbors(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <- FindClusters(Dog.combined.singlet.CD8NKT, resolution = 1.4)
DefaultAssay(Dog.combined.singlet.CD8NKT) <- "RNA"
Dog.combined.singlet.CD8NKT <-ScaleData(Dog.combined.singlet.CD8NKT, verbose = FALSE)
Dog.combined.singlet.CD8NKT <- JoinLayers(Dog.combined.singlet.CD8NKT)
saveRDS(Dog.combined.singlet.CD8NKT, "Dog.combined.singlet.CD8NKT.rds")

DimPlot(Dog.combined.singlet.CD8NKT, reduction = "umap", label = T) + NoLegend()
ggsave("CD8NKT_UMAP.png", width = 5, height = 5)
DimPlot(Dog.combined.singlet.CD8NKT, reduction = "tsne", label = T) + NoLegend()
ggsave("CD8NKT_TSNE.png", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.CD8NKT, 
            features = c("CD3E","CD4","CD8A","RHEX","ENSCAFG00000028587"), 
            order = F, 
            reduction = "umap", label = T)
#
Dog.combined.singlet.CD8NKT <- JoinLayers(Dog.combined.singlet.CD8NKT)
sce <- as.SingleCellExperiment(Dog.combined.singlet.CD8NKT)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, doublets) 

Dog.combined.CD8NKT.singlet <- subset(Dog.combined.singlet.CD8NKT, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.CD8NKT.singlet) <- "integrated"
Dog.combined.CD8NKT.singlet <-ScaleData(Dog.combined.CD8NKT.singlet, verbose = FALSE)
Dog.combined.CD8NKT.singlet <-RunPCA(Dog.combined.CD8NKT.singlet, npcs = 30, verbose = FALSE)
Dog.combined.CD8NKT.singlet <-RunUMAP(Dog.combined.CD8NKT.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD8NKT.singlet <- RunTSNE(Dog.combined.CD8NKT.singlet, dims = 1:30, reduction = "pca")
Dog.combined.CD8NKT.singlet <- FindNeighbors(Dog.combined.CD8NKT.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD8NKT.singlet <- FindClusters(Dog.combined.CD8NKT.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.CD8NKT.singlet) <- "RNA"
Dog.combined.CD8NKT.singlet <-ScaleData(Dog.combined.CD8NKT.singlet, verbose = FALSE)
Dog.combined.CD8NKT.singlet <- JoinLayers(Dog.combined.CD8NKT.singlet)

DimPlot(Dog.combined.singlet.CD8NKT, reduction = "umap", pt.size = 0.1, 
        label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/CD8/CD8NKT_Doublet_tsne.png", width = 4, height = 4)
CD8NKT.doublet.percentage <- table(Idents(Dog.combined.singlet.CD8NKT), Dog.combined.singlet.CD8NKT$db.class)
write.table(CD8NKT.doublet.percentage, "CD8NKT_db_%.csv", sep = ",", append=F)

#
DimPlot(Dog.combined.doublet, reduction = "tsne", label = T) #+ NoLegend()
ggsave("Doublet_TSNE_legend.png", width = 5, height = 5)

FeaturePlot(Dog.combined.CD8NKT.singlet, 
            features = c("PDCD1","LAG3","TIGIT","CXCR5","ICOS","BCL6","MKI67"), order = T)

#
# Subclustering for CD4 helper T cells
Dog.combined.singlet.CD4 <- subset(Dog.combined.singlet, idents = c(25,13,6,19,22,24,3,0,11,10,21,38,9,2,12)) 
DefaultAssay(Dog.combined.singlet.CD4) <- "integrated"
Dog.combined.singlet.CD4 <-ScaleData(Dog.combined.singlet.CD4, verbose = FALSE)
Dog.combined.singlet.CD4 <-RunPCA(Dog.combined.singlet.CD4, npcs = 30, verbose = FALSE)
Dog.combined.singlet.CD4 <-RunUMAP(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <-RunTSNE(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <- FindNeighbors(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <- FindClusters(Dog.combined.singlet.CD4, resolution = 1.2)
DefaultAssay(Dog.combined.singlet.CD4) <- "RNA"
Dog.combined.singlet.CD4 <-ScaleData(Dog.combined.singlet.CD4, verbose = FALSE)
Dog.combined.singlet.CD4 <-JoinLayers(Dog.combined.singlet.CD4)

DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", label = T) #+ NoLegend()
ggsave("CD4_tsne_legend.png", width = 4, height = 4)


FeaturePlot(Dog.combined.singlet.CD4, features = c("TOX","PDCD1"), reduction = "umap", label = T, order = T)
saveRDS(Dog.combined.singlet.CD4, "PATH_TO_FILE")

#
Dog.combined.singlet.CD4 <- JoinLayers(Dog.combined.singlet.CD4)
sce <- as.SingleCellExperiment(Dog.combined.singlet.CD4)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, doublets) 

Dog.combined.CD4.singlet <- subset(Dog.combined.singlet.CD4, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.CD4.singlet) <- "integrated"
Dog.combined.CD4.singlet <-ScaleData(Dog.combined.CD4.singlet, verbose = FALSE)
Dog.combined.CD4.singlet <-RunPCA(Dog.combined.CD4.singlet, npcs = 30, verbose = FALSE)
Dog.combined.CD4.singlet <-RunUMAP(Dog.combined.CD4.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD4.singlet <- RunTSNE(Dog.combined.CD4.singlet, dims = 1:30, reduction = "pca")
Dog.combined.CD4.singlet <- FindNeighbors(Dog.combined.CD4.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD4.singlet <- FindClusters(Dog.combined.CD4.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.CD4.singlet) <- "RNA"
Dog.combined.CD4.singlet <-ScaleData(Dog.combined.CD4.singlet, verbose = FALSE)
Dog.combined.CD4.singlet <- JoinLayers(Dog.combined.CD4.singlet)

DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("CD4_Doublet_tsne.png", width = 4, height = 4)

CD4.doublet.percentage <- table(Idents(Dog.combined.singlet.CD4), Dog.combined.singlet.CD4$db.class)
write.table(CD4.doublet.percentage, "CD4_db_%_.csv", sep = ",", append=F)

#
FeaturePlot(Dog.combined.singlet.CD4, features = c("SELL","CD4"),
            order = F, reduction = "umap" , label = F)

#
dotPlot(Dog.combined.singlet.CD4, features = c("ISG15"), rotation = T,angle.x = 45,
        colormap = "Blues", #idents = c(3,5,6),
        color.direction = 1) #+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/Dot_CD4_356.png", width = 5, height = 4) #width = 2, height = 2.7

dotPlot(Dog.combined.singlet.CD4, features = c("RPL8","LOC106558744","LOC119876465"), rotation = T,angle.x = 45,
        colormap = "Blues", #idents = c(3,5,6),
        color.direction = 1) #+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/Dot_CD4_356.png", width = 5, height = 4) #width = 2, height = 2.7

VlnPlot(Dog.combined.singlet.CD4, features = c("CCDC3"), pt.size = 0) +NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/Vlnplot_CCDC3.png", width = 5, height = 4) #width = 2, height = 2.7

#
FeaturePlot(Dog.combined.singlet.CD4, features = c("CCDC3"),
            order = F, reduction = "umap" ,label = T)
FeaturePlot(Dog.combined, features = c("CD3E","CD4","CSF3R","CD8A"), reduction = "tsne", order = T, label = T)
DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", label = T)
ggsave("Data_Healthy/Dog.combined.singlet.CD4.png", width = 5, height = 5)
CD4.doublet.percentage <- table(Idents(Dog.combined.CD4), Dog.combined.CD4$db.class)
write.table(CD4.doublet.percentage, "CD4_db_%_.csv", sep = ",", append=F)

Cellnumbers.subset <- table(Idents(Dog.combined.CD4.singlet), Dog.combined.CD4.singlet$orig.ident)
write.table(Cellnumbers.subset,
            "CD4_Cellnumber_orig.ident.csv", sep = ",", append=F)

DimPlot(Dog.combined.CD4, group.by = "db.class", reduction = "tsne", cols = c("grey","red")) +NoLegend()
ggsave("CD4.db.tsne.png", width = 5, height = 5)

DimPlot(Dog.combined.singlet.CD4, reduction = "umap", label = T) + NoLegend()
ggsave("Data_Healthy/CD4_Singlet_Umap.png", width = 5, height = 5)

Idents(Dog.combined.singlet.CD4) <- Dog.combined.singlet.CD4$seurat_clusters
DefaultAssay(Dog.combined.singlet.CD4) <- "RNA"
Dog.combined.singlet.CD4 <-ScaleData(Dog.combined.singlet.CD4, verbose = FALSE, features = rownames(Dog.combined.singlet.CD4)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet.CD4 <- JoinLayers(Dog.combined.singlet.CD4)
Dog.markers <- FindAllMarkers(Dog.combined.singlet.CD4, only.pos = TRUE, min.pct = 0.36, logfc.threshold = 0.36) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top50, "Data_Healthy/Featureplot_Cluster_Marker/CD4/Top20_0.5.csv", append = "F",sep = ",")
write.table(top50, "PATH_TO_FILE", append = "F",sep = ",")

grep("FOXP3",top20)
dotPlot(Dog.combined.singlet.CD4, features = c("LOC119876465"))
Dog.combined.singlet.CD4@active.ident <- factor(Dog.combined.singlet.CD4@active.ident, 
                                                levels = c("T_CD4", "T_CD8","B","Monocytes","DC","PMN","T_Gd","T","Un"))
DoHeatmap(subset(Dog.combined.singlet.CD4, downsample = 100),
          features = top5$gene, 
          #disp.max = 2,
          #  disp.min = -2,
          size = 3, angle = 50, cells = 1:50)+ NoLegend() #100 would be the best
theme(text = element_text(size=8)) #+NoLegend()
# scale_fill_viridis(option = "E")#+ 
# labs(x = "Viridis E", y = NULL) #+ NoLegend() #+ scale_fill_viridis(100) #+ NoLegend()
ggsave("PATH_TO_FILE", width=8, height=7)
dotPlot(Dog.combined.singlet.CD4, features = c("CD3E","CD5","CD4","CD8A","CD8B","LOC608848","LOC490269",
                                               "CSF3R","DPYD","LOC100686073","GRM7","IL2RA","FOXP3","IL7R",
                                               "CTLA4","IKZF2","JUND",
                                               "CCR7","SELL","TCF7","LEF1","CD44","CD27","S100A4","S100A6","CD74",
                                               "ITGB1","ITGA4","CD28","TBX21","GATA3","RORC",
                                               "STAT1","STAT4","IL12RB2","STAT6","STAT3","IL2",
                                               "TNFRSF4","TNFRSF18","IL5","IL10",
                                               "IL13","IL17F","TNF","ICOS",
                                               "GZMK","GZMB","PRF1","IFNG",
                                               "PDCD1","rna_TOX","LAG3","HAVCR2"))
dotPlot(Dog.combined.singlet.CD4, 
        features = c("CSF3R","DPYD","CD8A","CD8B","CD5","CD3E","CD4",
                     "IL2RA","FOXP3","CTLA4","IKZF2","BATF","IL32",
                     "CCR7","SELL","TCF7","LEF1",
                     "CD44","CD27","IL7R","S100A4",
                     "ITGB1","ITGA4","CD28","TBX21","GATA3","RORC",
                     "STAT1","STAT4","IL12RB2",#"rna_IL2",
                     "TNFRSF4","TNFRSF18","rna_IL10","TNF","ICOS",
                     "GZMK","GZMB","PDCD1","rna_TOX","LAG3","HAVCR2",
                     "MX1","ISG15","CCDC3","STAT3","LOC119876465"),
        cluster.idents = T,
        colormap = "Blues"
) 
ggsave("Data_Healthy/CD4/Dot_Marker.png", width = 7, height = 10)

Dog.combined.singlet.CD4 <- RenameIdents(Dog.combined.singlet.CD4,
                                         '0' = "CD4", '1' = "CD4", '2' = "CD4",'3' = "CD4", '4' = "CD4", 
                                         '5' = "CD4", '6' = "CD4", '7' = "CD4", '8' = "CD4", '9'= "CD4", 
                                         '10' = "CD4", '11'= "CD4", '12' = "C12", '13' = "CD4", '14' = "C14", 
                                         '15' = "CD4", '16' = "CD4", '17' = "CD4")
Dog.combined.singlet.CD4@active.ident <- Dog.combined.singlet.CD4$seurat_clusters 


VlnPlot(Dog.combined.singlet.CD4, idents = c(7,17,5,2,12,6,3,9,16,14),
        #cols = c("#0348A6","#0348A6","#0348A6","#0348A6","#0348A6",
        #         "#0348A6","#0348A6","#0348A6","#0348A6","#E63222",
        #         "#E63222"), 
        features = c("MX1","IFI44","OAS1","CMPK2","LGALS9",
                     "ISG15","MX2","ISG20","XAF1"
                     #   "PDCD1","CTLA4"
        ), 
        stack = T, flip = T, sort = T)  + NoLegend() 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/C14_IFN_genes.png", width=3.8, height=3)

VlnPlot(Dog.combined.singlet.CD4, idents = c(7,17,5,2,12,6,3,9,16,14),
        cols = c("#0348A6","#0348A6","#0348A6"), 
        features = c("PDCD1","CTLA4","CCR4"), 
        stack = T, flip = T, sort = T)  + NoLegend() 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/C14_PDCD1_CTLA4_CCR4_genes.png", width=3.7, height=1.5)

VlnPlot(Dog.combined.singlet.CD4.Treg, features = "Exhaustion", group.by = "CCR4_Sep", cols = c("#E9F0FF","#0A5096"))
p <-VlnPlot(Dog.combined.singlet.CD4.Treg, 
            features = c("Tumor_Immune_Escape"),  # Tumor_Immune_Escape, Response_To_ICB
            group.by = "CCR4_Sep", 
            cols = c("#E9F0FF","#0A5096")) 
p + stat_compare_means(method = "t.test")
ggsave("Data_Healthy/CD4/CCR4_Treg_Response_To_ICB.png", width = 4, height = 4)


Dog.combined.singlet.CD4@active.ident <- Dog.combined.singlet.CD4$seurat_clusters
FeatureScatter(subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(14)), 
               feature1 = "T1_Interferon", feature2 = "Dog_MMVD",
               jitter = T,span = T, 
               cols = "#CD78F5",
               pt.size = 3,plot.cor = T, log = F)
ggsave("Data_Healthy/CD4/C14_scatterfeature_T1_IFN_MGT_Complex.png", width=6, height=5)
C14<-subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(14))
cor.test(C14$T1_Interferon, C14$Dog_MGT_Complex) # P= 5.843e-10
model <- lm(C14$T1_Interferon ~ C14$Dog_MMVD)
model <- lm(C14$T1_Interferon ~ C14$Dog_MGT_Complex)

summary(model)

FeatureScatter(subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(7,17)), 
               feature1 = "CTLA4", feature2 = "BATF",
               jitter = T,span = T, 
               pt.size = 3,plot.cor = T, log = F)


cells.to.include <- WhichCells(
  object = Dog.combined.singlet.CD4,
  expression = BATF > 0 & FOXP3 > 0,
  idents = c(7, 17)
)

FeatureScatter(
  object = subset(Dog.combined.singlet.CD4, cells = cells.to.include),
  feature1 = "BATF",
  feature2 = "FOXP3",
  jitter = TRUE,
  span = TRUE,
  pt.size = 3,
  plot.cor = TRUE,
  log = FALSE
)

FeatureScatter(subset(Dog.combined.singlet.myeloid, seurat_clusters %in% c(20)), 
               feature1 = "Dog_MMVD", feature2 = "Dog_IMHA",
               jitter = T,span = T, 
               cols = "#CD78F5",
               pt.size = 3,plot.cor = T, log = F)
ggsave("Data_Healthy/Myeloid//C20_scatterfeature_MMVD_IMHA.png", width=6, height=5)
C20 <-subset(Dog.combined.singlet.myeloid, seurat_clusters %in% c(20))
model <- lm(C20$Dog_IMHA ~ C20$Dog_MMVD)
summary(model)


DimPlot(Dog.combined.singlet.CD4, split.by = "PD1")
FeatureScatter(Dog.combined.singlet.CD4, 
               feature1 = "PDCD1", feature2 = "Proinflammatory",
               jitter = T,span = T, 
               #cols = "#CD78F5",
               pt.size = 3,plot.cor = T, log = F)
ggsave("Data_Healthy/CD4/C14_scatterfeature_T1_IFN_MGT_Complex.png", width=6, height=5)

C18<-subset(Dog.combined.singlet.myeloid, seurat_clusters %in% c(18))
cor.test(C18$Dog_MMVD, C18$T1_Interferon)
model <- lm(C18$Dog_MMVD ~ C18$T1_Interferon)
summary(model)

FeatureScatter(Dog.combined.singlet.myeloid, 
               feature1 = "T1_Interferon", feature2 = "Dog_MMVD",
               jitter = T,span = T, 
               #  cols = "#E86AB5",
               pt.size = 3, plot.cor = T, log = T) + NoLegend() 
cor.test(Dog.combined.singlet.myeloid$Dog_MMVD, Dog.combined.singlet.myeloid$T2_Interferon) # p-value < 2.2e-16

FeatureScatter(subset(Dog.combined.singlet.B, seurat_clusters %in% c(9)), 
               feature1 = "T1_Interferon", feature2 = "Dog_MMVD",
               jitter = T,span = T, 
               cols = "#E86AB5",
               pt.size = 3,plot.cor = T, log = F) + NoLegend() 
ggsave("Data_Healthy/B/C9_scatterfeature_T1_IFN_MMVD.png", width=6, height=5)

B9 <-subset(Dog.combined.singlet.B, seurat_clusters %in% c(9))
cor.test(B9$T1_Interferon, B9$Dog_MMVD) # P= 0.1723
model <- lm(B9$T1_Interferon ~ B9$Dog_MMVD)
summary(model)


FeatureScatter(subset(Dog.combined.singlet.myeloid, seurat_clusters %in% c(18)), 
               feature1 = "T2_Interferon", feature2 = "Dog_MMVD",
               jitter = T,span = T, 
               cols = "#CF77F2",
               pt.size = 3, plot.cor = T, log = F) + NoLegend() 
ggsave("Data_Healthy/Myeloid/C18_scatterfeature_T1_IFN_MMVD.png", width=6, height=5)
C18<-subset(Dog.combined.singlet.myeloid, seurat_clusters %in% c(18))
cor.test(C18$T1_Interferon, C18$Dog_MMVD) # P= 1.133e-08

FeaturePlot(Dog.combined.singlet.CD4, reduction = "umap", 
            features = c("Exhaustion"),
            pt.size = 1,
            label = F,order = T, #group.by = "db.class",
            cols = viridis(100), #c("#ADDAE6","#E63222")
) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/CD4_ISG15.png", width = 4, height = 4)

FeatureScatter(subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos")), 
               feature1 = "Adaptive_Immunity", feature2 = "Exhaustion",
               jitter = T, span = T, 
               pt.size = 3, plot.cor = T, log = F)
ggsave("Data_Healthy/CD4/PD1+CD4T_scatterfeature_Adaptive_Exhaustion.png", width=6, height=5)

VlnPlot(Dog.combined.singlet.CD4, group.by = "PD1", features = c("Proinflammatory"), 
        cols = c("#440154FF","#FDE725FF"))
#

#
unique(Dog.combined.singlet.CD4@active.ident)

CD4.cluster14 <- FindMarkers(Dog.combined.singlet.CD4, ident.1 = 14, min.pct = 0.25) 
CD8.cluster7 <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = 7, min.pct = 0.25) 

write.table(CD4.cluster14, "Data_Healthy/Featureplot_Cluster_Marker/CD4/cluster14.csv", append = F, sep = ",")
write.table(CD8.cluster7, "Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/cluster7.csv", append = F, sep = ",")

FeaturePlot(Dog.combined.singlet.CD4, 
            features = c("Th17"), 
            cols = viridis(100), pt.size = 1,
            order = F) #+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/Density_Th2_Lang.png", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.CD4, 
            features = c("Th17"), 
            cols = c("gray90", "red"),
            #   cols = viridis(100),pt.size = 1,
            min.cutoff = "q1",
            max.cutoff = "q99",
            order =T)# + NoLegend()
min(Dog.combined.singlet.CD4[["RNA"]]@misc["Th17", ])

ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/FOXP3.png", width = 9, height = 3)

FeaturePlot(Dog.combined.singlet.CD4, 
            features = c("PDCD1"), 
            #  cols = c("#ADDAE6","#E63222"), 
            cols = plasma(100), 
            order = T)
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD4/.png", width = 12, height = 12)

# monaco, nover, blue, hpca, immune, immgen
DimPlot(subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(8)),
        reduction = "umap", 
        group.by = "hpca", 
        # label = T, 
        shuffle = T, 
        #repel = T, 
        label.box = T, 
        #label.size = 6,
        #pt.size = 1,
        cols = colorRampPalette(brewer.pal(12, "Paired"))(11))#+ NoLegend()
ggsave("PATH_TO_FILE", width = 12, height = 12)
dev.off()
library(Nebulosa)
library(RColorBrewer)
library(viridis)
Dog.combined.singlet.CD4$CD4
FeaturePlot(Dog.combined.singlet.CD4, reduction = "umap", features = c("TIGIT"),
            pt.size = 1,
            label = F,order = T, #group.by = "db.class",
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/CD4/CD4_TIGIT.png", width = 4, height = 4)
CD4.doublet.percentage <- table(Idents(Dog.combined.singlet.CD4), Dog.combined.singlet.CD4$db.class)
write.table(CD4.doublet.percentage, "Data_Healthy/CD4_db_%_.csv", sep = ",", append=F)

FeaturePlot(Dog.combined.singlet.CD4, reduction = "umap", features = c("TIGIT"),
            pt.size = 1,
            label = F,order = T, #group.by = "db.class",
            cols = viridis(100)
) + NoLegend()
ggsave("Data_Healthy/CD4/Density_CD4_IFN.png", width = 5, height = 5)

plot_density(Dog.combined.singlet.CD4, reduction = "umap",
             features = c("CTLA4"), 
             size = 1, pal = "plasma",# "cividis", "inferno", "plasma"
) #+NoLegend()
ggsave("Data_Healthy/CD4/Density_TOX.png",
       width = 5, height = 5)

Dog.combined.singlet.CD4$pdcd1 <- Dog.combined.singlet.CD4@assays$RNA$data["PDCD1",]
Dog.combined.singlet.CD4$PDL1_Sep <- ifelse(Dog.combined.singlet.CD4$cd274>0, "Pos","Neg")

FeaturePlot(Dog.combined.singlet.CD4, features = c("IL32"), #pt.size = 0.1, 
            order = F, cols =  c("#ADDAE6","#E63222"),
            # split.by = "PDL1_Sep", 
            label = T)
ggsave("Data_Healthy/CD4/Feature_IL32.png", width = 5, height = 5)
Dog.combined.singlet.CD4@active.ident <- Dog.combined.singlet.CD4$seurat_clusters

DimPlot(Dog.combined.singlet.CD4, group.by = "PDL1_Sep", pt.size = 2,
        cols =  c("#ADDAE6","#E63222"), order = T) + NoLegend()
ggsave("Data_Healthy/Umap_CD4_group.by_PDL1.png", width = 5, height = 5)

cor.test(subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(14),Dog.combined.singlet.CD4$Dog_MMVD,Dog.combined.singlet.CD4$T1_Interferon ))
cor.test.subset <- subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(14))

cor.test(cor.test.subset$Dog_MMVD,cor.test.subset$T1_Interferon)
cor.test(cor.test.subset$Dog_MMVD,cor.test.subset$Dog_MGT_Complex)

#
DefaultAssay(Dog.combined.singlet.CD8NKT) <- "RNA"
Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$seurat_clusters
Dog.combined.singlet.CD8NKT <-ScaleData(Dog.combined.singlet.CD8NKT, verbose = FALSE, features = rownames(Dog.combined.singlet.CD8NKT)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet.CD8NKT <- JoinLayers(Dog.combined.singlet.CD8NKT)
Dog.markers <- FindAllMarkers(Dog.combined.singlet.CD8NKT, only.pos = TRUE, min.pct = 0.36, logfc.threshold = 0.36) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top20, "Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/Top20_FC.csv", append = "F",sep = ",")
write.table(top50, "PATH_TO_FILE", append = "F",sep = ",")

VlnPlot(Dog.combined.singlet.CD8NKT, features = c("PDCD1","IL10","FOXP3","TIGIT","HAVCR2","LAG3"))
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD413/Vln.png", width = 10, height = 10)
Dog.combined.singlet.CD8NKT@active.ident <- factor(Dog.combined.singlet.CD8NKT@active.ident, 
                                                   levels = c("T_CD8NKT", "T_CD8","B","Monocytes","DC","PMN","T_Gd","T","Un"))
DoHeatmap(subset(Dog.combined.singlet.CD8NKT, downsample = 100),
          features = top10$gene,
          size = 3, angle = 50, cells = 1:50)+ NoLegend() #100 would be the best
theme(text = element_text(size=8)) #+NoLegend()
# scale_fill_viridis(option = "E")#+ 
# labs(x = "Viridis E", y = NULL) #+ NoLegend() #+ scale_fill_viridis(100) #+ NoLegend()
ggsave("PATH_TO_FILE", width=8, height=7)
library(CellChat)
dotPlot(Dog.combined.singlet.CD8NKT, features = c("CD3D","CD247","CD5","CD4","CD8A","CD8B",
                                                  "LOC491694","LOC611565","RHEX",
                                                  "CCR7","SELL","LEF1","TCF7",
                                                  "IL7R","CD7","CD27","ITGB1","CD69",
                                                  "FCER1G","PI3","ITGAE","CXCR6",
                                                  "CD44","CD28","ZNF683","PRDM1",
                                                  "CD160","GZMA","GZMK","GZMB","PRF1",
                                                  "PDCD1","CTLA4","ID3","TOX","LAG3","HAVCR2","FOXP3","IL2RA","GATA3",
                                                  "LOC102155278","LOC111090118",
                                                  "MKI67","PCLAF","CSF3R","DPYD"),
        cluster.idents = T,
        colormap = "Blues") 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/Dot_Marker.png", width = 7, height = 8.5)
gene <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = c(10,11), min.pct = 0.5, only.pos = T)
head(gene, 5)
FeaturePlot(Dog.combined.singlet.CD8NKT, 
            features = c("T_Cell_Proliferaiton"), order = F, #pt.size = 1,
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8_CTLA4.png", width = 5, height = 5)

mean(Dog.combined.singlet.CD8NKT$T_Cell_Terminal_Differentiation)
mean(Dog.combined.singlet.CD8NKT$Exhaustion)

FeaturePlot(Dog.combined.singlet.CD8NKT, 
            features = c("T_Cell_Terminal_Differentiation","Exhaustion"), blend = T,
            order = T,
) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8_T_Diff_Exh.png", width = 8, height = 2.5)

FeaturePlot(subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(9,10,11,12,13)), 
            features = c("T_Cycle"), order = T,label = T,label.size = 7,
            #   pt.size = 0.1,
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/C910111213_CD247.png", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.CD8NKT,
            features = c("T_Cycle"), order = T,#label = T,label.size = 7,
            pt.size = 1,
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/T_Cycle.png", width = 5, height = 5)

plot_density(Dog.combined.singlet.CD8NKT, reduction = "umap",
             features = c("CCR7","LEF1","PRF1","GZMB","GZMA","GZMK"), pal = "inferno")# + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/Density_Niave_Effector.png", width = 8, height = 4)

plot_density(Dog.combined.singlet.CD8NKT, reduction = "umap",
             features = c("RHEX"), pal = "inferno") + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/Density_STAT3.png", width =5, height = 5)

FeatureScatter(Dog.combined.singlet.myeloid, feature1 = c("Dog_MMVD"), feature2 = c("T1_Interferon"))

Dog.combined.singlet[[]]
DimPlot(subset(Dog.combined.singlet.CD8NKT, seurat_clusters%in% c(10,11)),
        reduction = "umap", group.by  = "blue",  # no immgen, blue, hpca, nover, monaco
        label = T, shuffle = T, repel = T, label.box = T, 
        label.size = 2,
        cols = colorRampPalette(brewer.pal(12, "Paired"))(22))+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/SingleR_monaco.png", width = 9, height = 9)

meta <- Dog.combined.singlet.CD8NKT[[]] 
View(meta)
#meta <- meta[,c(4,16,22:96,110,117:124,145:182,189)] 
#meta <- meta[,c(4,16,86:96,117:118,189)] 
#meta <- meta[,c(4,16,27,37,44,38,146,55,180)] 
meta <- meta[,c(4,16,66,75,76,82,93,124,125, 129:132)] 

meta <- meta[,c(4,16,27,26,37,38,44,55,131:132)] 


View(meta)
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:8, mean)) #122
View(heatmap)
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C",0:18) 
headers <- heatmap[,c(1:2)] 
#heatmap <- heatmap[c(1:9,14,15),c(3:124)] # 124
#heatmap <- heatmap[c(1:16,18,19),c(3:9)] # 124
heatmap <- heatmap[,c(3:10)] # 124

SubsetColors <- c("darkred") 
names(SubsetColors) <- c("Healthy") 
clusterColors <- c(scales::hue_pal()(19))
names(clusterColors) <- 0:18

colors <- list(Subset = SubsetColors, seurat_clusters = clusterColors)
normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:18) #rownames(heatmap2) <- paste0("X", 1:22)
#rownames(heatmap2) <- paste0("C",c(1:16,18,19)) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/GSEA_CD8_selected_small3.pdf", width = 8, height =7) 
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols=T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   cellheight = 12,
                   #    color = inferno(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)) 
)
dev.off()


#

meta <- Dog.combined.singlet.CD4[[]] 
View(meta)
meta <- meta[,c(16,58:63,66:68,70:75,191)] 
View(meta)
heatmap <- meta %>% 
  group_by(seurat_clusters, CTLA4_Sep) %>% 
  summarise(across(1:15, mean)) #122
View(heatmap)
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C",0:35) 
headers <- heatmap[,c(1:2)] 
heatmap <- heatmap[,c(3:17)] # 124

SubsetColors <- c("darkred") 
names(SubsetColors) <- c("Healthy") 
Ctla4Colors <- c("#FCE540", "#3F114E") 
names(Ctla4Colors) <- c("Pos","Neg") 
clusterColors <- c(scales::hue_pal()(19))
names(clusterColors) <- 0:18

colors <- list(Subset = SubsetColors, 
               CTLA4_Sep = Ctla4Colors,
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:35) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/GSEA_CD4_CTLA4.pdf", width = 12, height =7) 
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols=T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   cellheight = 12,
                   #    color = inferno(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)) 
)
dev.off()
FeatureScatter(Dog.combined.singlet.CD4, 
               feature1 = c("T_Cell_Terminal_Differentiation"),
               feature2 = c("CTLA4"))

#
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <- read.csv("PATH_TO_FILE")
genes <- read.csv("PATH_TO_FILE")
genes <- read.csv("PATH_TO_FILE")

counts <- GetAssayData(Dog.combined.singlet.CD8NKT, assay = "RNA", slot = "data")
counts <- t(as.matrix(counts[rownames(counts) %in% genes$Genes,])) 
meta <- Dog.combined.singlet.CD8NKT[[]][,c("Subset", "seurat_clusters")]

counts <- data.frame(meta, counts)
ncol(counts)
View(counts)
heatmap <- counts %>%
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:13, mean)) #
heatmap <- as.data.frame(heatmap)

#need to make rownames to match the heatmap to the annotation
nrow(heatmap)
View(heatmap)
rownames(heatmap) <- paste0("C", 0:18) # rownames will be generated as X1, X2, ..., X35
headers <- heatmap[,1:2]
heatmap <- heatmap[c(1,2,3,4,5,6,7,8,9,15),3:7] # All T
heatmap <- heatmap[,c(3:15)] # All T

#Defining the color scheme
SubsetColors <-  c("#FF4B20","blue")
names(SubsetColors) <- c("Healthy","Healthy")
clusterColors <- c(scales::hue_pal()(19)) 
names(clusterColors) <- 0:18
colors <- list(Subset = SubsetColors, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", c(1,2,3,4,5,6,7,8,9,15))
rownames(heatmap2) <- paste0("C", 0:18)

#Plotting non-scaled version
pdf("Data_Healthy/CD8NKT_Terminal_Genes.pdf", width = 10, height = 8)
pdf("Data_Healthy/CD8NKT_TumorEscapeGene.pdf", width = 10, height = 8)
pdf("Data_Healthy/CD8NKT_ResposnetoICIgene.pdf", width = 10, height = 8)

pheatmap::pheatmap(t(heatmap2),show_colnames = T, scale = "row",
                   annotation_col = headers, cluster_cols= T, cluster_rows = T,  
                   color = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
                   annotation_colors = colors, fontsize = 12, 
                   cellwidth = 10, 
                   cellheight = 12,
                   legend = T, legend_labels = T, annotation_legend = T)
dev.off()
VlnPlot(Dog.combined.singlet.CD8NKT, 
        features = c("CXCR6"), sort= T, pt.size = 0
        # stack = T,flip = T, sort = F
) +NoLegend() #cols = c("#FF4B20", "#FFB433", "#7AC5FF", "#C6FDEC", "#0348A6","#3F114E")
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/CXCR6_Tgd_vlnplot.png", width = 4.5, height = 2.5)
dotPlot(Dog.combined.singlet.CD8NKT, cluster.idents = T,
        features = c("TAPBP","CXCR6")) #"DLA-DQA1","HLA-DRB1",
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8/CXCR6_Tgd_DOTplot_LEGEND.png", width = 6, height = 4)

genes <- read.csv("PATH_TO_FILE")
counts <- GetAssayData(Dog.combined.singlet.CD8NKT, assay = "RNA", slot = "data")
counts <- t(as.matrix(counts[rownames(counts) %in% genes$Genes,])) 
meta <- Dog.combined.singlet.CD8NKT[[]][,c("PD1_LAG3", "Subset")]

counts <- data.frame(meta, counts)
ncol(counts)
View(counts)
heatmap <- counts %>%
  group_by(Subset, PD1_LAG3) %>% 
  summarise(across(1:5, mean)) #
heatmap <- as.data.frame(heatmap)

#need to make rownames to match the heatmap to the annotation
nrow(heatmap)
View(heatmap)
rownames(heatmap) <- paste0("X", 1:2) # rownames will be generated as X1, X2, ..., X35
headers <- heatmap[,1:2]
heatmap <- heatmap[,3:7] # All T

#Defining the color scheme
SubsetColors <-  c("#FF4B20","blue")
names(SubsetColors) <- c("Healthy","Healthy")
PD1LAG3Colors <-  c("#FF4B20","blue")
names(PD1LAG3Colors) <- c("Pos","Neg")
clusterColors <- c(scales::hue_pal()(19)) 
names(clusterColors) <- 0:18
colors <- list(Subset = SubsetColors, 
               PD1_LAG3 = PD1LAG3Colors, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", c(1:2))

#Plotting non-scaled version
pdf("Data_Healthy/CD8NKT_Terminal_Genes.pdf", width = 10, height = 8)
pheatmap::pheatmap(t(heatmap2),show_colnames = F, scale = "none",
                   annotation_col = headers, cluster_cols= T, cluster_rows = T,  
                   color = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
                   annotation_colors = colors, fontsize = 12, 
                   cellwidth = 10, 
                   cellheight = 12,
                   legend = T, legend_labels = T, annotation_legend = T)
dev.off()

#
Dog.combined.singlet.myeloid <- subset(Dog.combined.singlet, 
                                       idents = c(29,14,17,30,20,15,18,26,
                                                  43,32,8,23,
                                                  34,42,41,46,
                                                  16,27))

DefaultAssay(Dog.combined.singlet.myeloid) <- "integrated"
Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE)
Dog.combined.singlet.myeloid <-RunPCA(Dog.combined.singlet.myeloid, npcs = 30, verbose = FALSE)
Dog.combined.singlet.myeloid <-RunUMAP(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <-RunTSNE(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <- FindNeighbors(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <- FindClusters(Dog.combined.singlet.myeloid, resolution = 1.2)
DefaultAssay(Dog.combined.singlet.myeloid) <- "RNA"
Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE)
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
saveRDS(Dog.combined.singlet.myeloid, "Dog.combined.singlet.myeloid.rds")

DimPlot(Dog.combined.singlet.myeloid, reduction = "umap", label = T) + NoLegend()
ggsave("myeloid_UMA.png", width = 5, height = 5)
DimPlot(Dog.combined.singlet.myeloid, reduction = "tsne", label = T) + NoLegend()
ggsave("myeloid_TSNE.png", width = 5, height = 5)
FeaturePlot(Dog.combined.singlet.myeloid, features = c("CD8A","DPYD"))

Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
sce <- as.SingleCellExperiment(Dog.combined.singlet.myeloid)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, doublets) 

Dog.combined.myeloid.singlet <- subset(Dog.combined.singlet.myeloid, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.myeloid.singlet) <- "integrated"
Dog.combined.myeloid.singlet <-ScaleData(Dog.combined.myeloid.singlet, verbose = FALSE)
Dog.combined.myeloid.singlet <-RunPCA(Dog.combined.myeloid.singlet, npcs = 30, verbose = FALSE)
Dog.combined.myeloid.singlet <-RunUMAP(Dog.combined.myeloid.singlet, reduction = "pca", dims = 1:30)
Dog.combined.myeloid.singlet <- RunTSNE(Dog.combined.myeloid.singlet, dims = 1:30, reduction = "pca")
Dog.combined.myeloid.singlet <- FindNeighbors(Dog.combined.myeloid.singlet, reduction = "pca", dims = 1:30)
Dog.combined.myeloid.singlet <- FindClusters(Dog.combined.myeloid.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.myeloid.singlet) <- "RNA"
Dog.combined.myeloid.singlet <-ScaleData(Dog.combined.myeloid.singlet, verbose = FALSE)
Dog.combined.myeloid.singlet <- JoinLayers(Dog.combined.myeloid.singlet)

DimPlot(Dog.combined.singlet.myeloid, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Myeloid_Doublet_tsne.png", width = 4, height = 4)
Myeloid.doublet.percentage <- table(Idents(Dog.combined.singlet.myeloid), Dog.combined.singlet.myeloid$db.class)
write.table(Myeloid.doublet.percentage, "Myeloid_db_%_.csv", sep = ",", append=F)

#
DimPlot(Dog.combined.singlet.myeloid, reduction = "umap", label = T)
ggsave("Data_Healthy/myeloid_umap.png", width = 5, height = 5)
FeaturePlot(Dog.combined.singlet.myeloid, features = c("CD24"), 
            #pt.size = 1, 
            reduction = "umap", order = T,
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/Myeloid/CD24.png",
       width = 4, height = 4)

plot_density(Dog.combined.singlet.myeloid, reduction = "umap",
             features = c("ITGAX"), size = 0.5, pal = "plasma",# "cividis", "inferno", "plasma"
) +NoLegend()
ggsave("Data_Healthy/Myeloid/CD4_Density.png",
       width = 4, height = 4)

Dog.markers <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = c(17), only.pos = FALSE, #min.pct = 0.3, #logfc.threshold = 0.5
)
write.table(Dog.markers, "Data_Healthy/Featureplot_Cluster_Marker/CD8/Tgd_C17_DEG.csv", append = "F",sep = ",")

VlnPlot(Dog.combined.singlet.CD8NKT, features = c("TAPBP"))
head(Dog.markers, 10)
DefaultAssay(Dog.combined.singlet.myeloid) <- "RNA"
Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE, features = rownames(Dog.combined.singlet.myeloid)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
Dog.markers <- FindAllMarkers(Dog.combined.singlet.myeloid, only.pos = TRUE, min.pct = 0.36, logfc.threshold = 0.36) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top20, "PATH_TO_FILE", append = "F",sep = ",")

DotPlot(Dog.combined.singlet.myeloid, f
        
)
Dog.combined.singlet.myeloid@active.ident <- factor(Dog.combined.singlet.myeloid@active.ident, 
                                                    levels = c("T_CD4", "T_CD8","B","Monocytes","DC","PMN","T_Gd","T","Un"))
DoHeatmap(subset(Dog.combined.singlet, downsample = 200),
          features = top5$gene, 
          #disp.max = 2,
          #  disp.min = -2,
          size = 3, angle = 50, cells = 1:200)+ NoLegend() #100 would be the best
theme(text = element_text(size=8)) #+NoLegend()
# scale_fill_viridis(option = "E")#+ 
# labs(x = "Viridis E", y = NULL) #+ NoLegend() #+ scale_fill_viridis(100) #+ NoLegend()
ggsave("PATH_TO_FILE", width=8, height=7)


DimPlot(Dog.combined.singlet.myeloid, label = T)

Dog.combined.singlet.myeloid <- RenameIdents(Dog.combined.singlet.myeloid,
                                             '0' = "0", '1' = "0", '2' = "0",'3' = "1", '4' = "1", 
                                             '5' = "1", '6' = "1", '7' = "1", '8' = "1", '9'= "2", 
                                             '10' = "2", '11'= "2", '12' = "2", '13' = "2", '14' = "6", 
                                             '15' = "3", '16' = "7", '17' = "6", 
                                             '18' = "3", '19' = "4", '20' = "5", '21' = "5",
                                             '22' = "4")
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters 
Dog.combined.singlet.myeloid$seurat_clusters <- Dog.combined.singlet.myeloid@active.ident
Dog.combined.singlet.myeloid$seurat_clusters2 <- Dog.combined.singlet.myeloid@active.ident
DimPlot(Dog.combined.singlet.myeloid)

library(escape)
ES.immune.T <- escape.matrix(Dog.combined.singlet.T,  # run this code 
                             gene.sets = list,
                             min.size = 3,
                             method = "ssGSEA")
saveRDS(ES.immune.T, "Data_Healthy/ES.immune.singlet.T.rds")

Dog.combined.singlet.T <- runEscape(Dog.combined.singlet.T, 
                                    method = "ssGSEA",
                                    normalize = T,
                                    gene.sets = list, 
                                    groups = 5000, 
                                    min.size = 3,
                                    new.assay.name = "escape.ssGSEA")
Dog.combined.singlet.T <- AddMetaData(Dog.combined.singlet.T, ES.immune.T)

saveRDS(Dog.combined.singlet.T, "PATH_TO_FILE")

FeatureScatter(Dog.combined.singlet.CD8NKT, shuffle = T, plot.cor = T, jitter = T, log = T,
               span = T, cols = viridis(38),
               feature1 = "T_Cell_Terminal_Differentiation",
               feature2 = "CTLA4") +NoLegend()
ggsave("Data_Healthy/Exh_adaptiveimmunity_CD4.png", width = 3, height = 3)
geyserEnrichment(Dog.combined.singlet.CD4, 
                 assay = "escape.ssGSEA",order.by = "mean",
                 scale = T,
                 gene.set = "Exhaustion", 
) #+ stat_compare_means(method = "anova")

#
ES.immune4 <- escape.matrix(Dog.combined.singlet.CD4,  # run this code 
                            gene.sets = list,
                            min.size = 3,
                            method = "ssGSEA")
saveRDS(ES.immune4, "Data_Healthy/ES.immune.singlet.CD8NKT.rds")

Dog.combined.singlet.CD4 <- runEscape(Dog.combined.singlet.CD4, 
                                      method = "ssGSEA",
                                      normalize = T,
                                      gene.sets = list, 
                                      groups = 5000, 
                                      min.size = 3,
                                      new.assay.name = "escape.ssGSEA")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, ES.immune4)

#
ES.immune4 <- escape.matrix(Dog.combined.singlet.myeloid,  # run this code 
                            gene.sets = list,
                            min.size = 3,
                            method = "ssGSEA")
saveRDS(ES.immune4, "Data_Healthy/ES.immune.singlet.myeloid.rds")

Dog.combined.singlet.myeloid <- runEscape(Dog.combined.singlet.myeloid, 
                                          method = "ssGSEA",
                                          normalize = T,
                                          gene.sets = list, 
                                          groups = 5000, 
                                          min.size = 3,
                                          new.assay.name = "escape.ssGSEA")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, ES.immune4)

#
ES.immune5 <- escape.matrix(Dog.combined.singlet.B,  # run this code 
                            gene.sets = list,
                            min.size = 3,
                            method = "ssGSEA")
saveRDS(ES.immune5, "Data_Healthy/ES.immune.singlet.B.rds")
Dog.combined.singlet.B[[]]
Dog.combined.singlet.B <- runEscape(Dog.combined.singlet.B, 
                                    method = "ssGSEA",
                                    normalize = T,
                                    gene.sets = list, 
                                    groups = 5000, 
                                    min.size = 3,
                                    new.assay.name = "escape.ssGSEA")

ES.immune6 <- escape.matrix(Dog.combined.singlet.B,  # run this code 
                            gene.sets =gene.list,
                            min.size = 3,
                            method = "ssGSEA")
saveRDS(ES.immune, "Data_Healthy/ES.immune.singlet.B.rds")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, ES.immune6)
Dog.combined.singlet.B <- runEscape(Dog.combined.singlet.B, 
                                    method = "ssGSEA",
                                    normalize = T,
                                    gene.sets = list, 
                                    groups = 5000, 
                                    min.size = 3,
                                    new.assay.name = "escape.ssGSEA")
#
saveRDS(ES.immune, "Data_Healthy/ES.immune.singlet.CD4.rds")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, ES.immune2) # add metadata
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, ES.immune4) # add metadata
Dog.combined.doublet <- AddMetaData(Dog.combined.doublet, ES.immune) # add metadata

#
VlnPlot(subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(7,17)), 
        features = c("Response_To_ICB")) + 
  stat_compare_means("t.test")

Idents(Dog.combined.singlet.CD4) <- Dog.combined.singlet.CD4$CCR4_Sep
Idents(Dog.combined.singlet.CD4.Treg) <- Dog.combined.singlet.CD4.Treg$CCR4_Sep
t.test(Dog.combined.singlet.CD4.Treg$Treg, alternative = "two.sided", var.equal = T)

Idents(Dog.combined.singlet.CD4.Treg) <- Dog.combined.singlet.CD4.Treg@active.ident
DimPlot(Dog.combined.singlet.CD4.Treg,cols = c("#095096","#E8F0FE"),
        shuffle = T, split.by = "CCR4_Sep") + NoLegend()
ggsave("Data_Healthy/CD4/CCR4+Treg_Umap.png", width = 4, height = 2)
VlnPlot(Dog.combined.singlet.CD4.Treg, features = c("Treg"))

Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, ES.immune) # add metadata
meta <- Dog.combined.singlet.CD4[[]] 
View(meta)
meta <- meta[,c(4,16,21:57)] 
meta <- meta[,c(4,16,58,60,61,63,67,68,73,74,87,132,154,155,176:181)] 
meta <- meta[,c(4,16,62,63,67,68,73,74,122:124,155,189,190, 192)] 

View(meta)
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset, CCR4_Sep) %>% 
  summarise(across(1:12, mean))
View(heatmap)
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C",0:17) 
rownames(heatmap) <- paste0("C",1:33) 

headers <- heatmap[,c(1:2)] 
headers <- heatmap[,c(1:3)] 

heatmap <- heatmap[c(1:10,12:18),c(3:14)]
heatmap <- heatmap[,c(4:15)]

SubsetColors <- c("darkred") 
names(SubsetColors) <- c("Healthy") 
CCR4Colors <- c("#095096","#E8F0FE") 
names(CCR4Colors) <- c("Pos","Neg") 
clusterColors <- c(scales::hue_pal()(18))
names(clusterColors) <- 0:17

colors <- list(Subset = SubsetColors, seurat_clusters = clusterColors, CCR4_Sep = CCR4Colors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:17) #rownames(heatmap2) <- paste0("X", 1:22)
rownames(heatmap2) <- paste0("C", c(0:9,11:17)) #rownames(heatmap2) <- paste0("X", 1:22)
rownames(heatmap2) <- paste0("C", c(1:33)) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/CD4/GSEA_CD4_CCR4.pdf", width = 10, height =7) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols=T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   cellheight = 12,
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()
DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", label = T)
FeaturePlot(Dog.combined.singlet.CD4, features = c("CD5"), cols = viridis(100), order = T)
library(RColorBrewer)

#
ES.immune <- readRDS("Data_Healthy/ES.immune.singlet.myeloid.rds")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, ES.immune) # add metadata

meta <- Dog.combined.singlet.myeloid[[]] 
View(meta)
meta <- meta[,c(4,16,24:30,32,33,81,97,120)] 
meta <- meta[,c(4,16,27:30,32,33,81,97,120)] 

meta <- meta[,c(4,16,24:30,32,33,35,43,66,74,72,75,81,95:97,114,124,125)] 
meta <- meta[,c(4,16,87:97,124:125)] 

View(meta)
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:9, mean))
View(heatmap)
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C",0:22) 
headers <- heatmap[,c(1:2)] 
heatmap <- heatmap[,c(3:11)]

SubsetColors <- c("darkred") 
names(SubsetColors) <- c("Healthy") 
clusterColors <- c(scales::hue_pal()(23))
names(clusterColors) <- 0:22

colors <- list(Subset = SubsetColors, seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:22) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/Myeloid/GSEA_Myeloid_2.pdf", width = 10, height =10) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols=T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()
library(RColorBrewer)



#
Dog.combined.singlet.B <- subset(Dog.combined.singlet, 
                                 idents = c(4,40,36,47))
DefaultAssay(Dog.combined.singlet.B) <- "integrated"
Dog.combined.singlet.B <-ScaleData(Dog.combined.singlet.B, verbose = FALSE)
Dog.combined.singlet.B <-RunPCA(Dog.combined.singlet.B, npcs = 30, verbose = FALSE)
Dog.combined.singlet.B <-RunUMAP(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <-RunTSNE(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <- FindNeighbors(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <- FindClusters(Dog.combined.singlet.B, resolution = 1)
DefaultAssay(Dog.combined.singlet.B) <- "RNA"
Dog.combined.singlet.B <-ScaleData(Dog.combined.singlet.B, verbose = FALSE)
Dog.combined.singlet.B <- JoinLayers(Dog.combined.singlet.B)
saveRDS(Dog.combined.singlet.B, "Dog.combined.singlet.B.rds")

DimPlot(Dog.combined.singlet.B, reduction = "umap", label = T)# + NoLegend()
ggsave("B_UMA_legend.png", width = 5, height = 5)
DimPlot(Dog.combined.singlet.B, reduction = "tsne", label = T) + NoLegend()
ggsave("B_TSNE.png", width = 5, height = 5)
FeaturePlot(Dog.combined.singlet.B, features = c("CD8A","DPYD"))

Dog.combined.singlet.B <- JoinLayers(Dog.combined.singlet.B)
sce <- as.SingleCellExperiment(Dog.combined.singlet.B)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, doublets) 

Dog.combined.B.singlet <- subset(Dog.combined.singlet.B, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.B.singlet) <- "integrated"
Dog.combined.B.singlet <-ScaleData(Dog.combined.B.singlet, verbose = FALSE)
Dog.combined.B.singlet <-RunPCA(Dog.combined.B.singlet, npcs = 30, verbose = FALSE)
Dog.combined.B.singlet <-RunUMAP(Dog.combined.B.singlet, reduction = "pca", dims = 1:30)
Dog.combined.B.singlet <- RunTSNE(Dog.combined.B.singlet, dims = 1:30, reduction = "pca")
Dog.combined.B.singlet <- FindNeighbors(Dog.combined.B.singlet, reduction = "pca", dims = 1:30)
Dog.combined.B.singlet <- FindClusters(Dog.combined.B.singlet, resolution = 1)
DefaultAssay(Dog.combined.B.singlet) <- "RNA"
Dog.combined.B.singlet <-ScaleData(Dog.combined.B.singlet, verbose = FALSE)
Dog.combined.B.singlet <- JoinLayers(Dog.combined.B.singlet)

DimPlot(Dog.combined.singlet.B, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("B_Doublet_tsne.png", width = 4, height = 4)
B.doublet.percentage <- table(Idents(Dog.combined.singlet.B), Dog.combined.singlet.B$db.class)
write.table(B.doublet.percentage, "B_db_%_.csv", sep = ",", append=F)

ES.immune <- escape.matrix(Dog.combined.singlet.B,  # run this code 
                           gene.sets = list, #list
                           min.size = 3,
                           method = "ssGSEA")
saveRDS(ES.immune, "Data_Healthy/ES.immune.singlet.B.rds")
Dog.combined.singlet.B[[]]
Dog.combined.singlet.B <- runEscape(Dog.combined.singlet.B, 
                                    method = "ssGSEA",
                                    normalize = T,
                                    gene.sets = list, #list
                                    groups = 5000, 
                                    min.size = 3,
                                    new.assay.name = "escape.ssGSEA")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, ES.immune)

meta <- Dog.combined.singlet.B[[]]
meta$seurat_clusters <- Dog.combined.singlet.B@active.ident
ncol(meta)
View(meta)
meta <- meta[,c(5,16,86:96,118,119,126,167,168)]  # dog's disease
meta <- meta[,c(5,16,94:96,118,119,126)]  # dog's disease
meta <- meta[,c(5,16,22:27,31:36,97:108)]  # dog's disease
meta <- meta[,c(5,16,68, 25:27,31:32,33,51,101,105,107)]  #GSEA
meta <- meta[,c(5,16,97:109, 132, 133, 157:161)]  #GSEA
meta <- meta[,c(5,16,26,27,31,32,105, 118,119, 94,96, 126,95)]  #GSEA

View(heatmap)
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:11, mean)) #1:13
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C", 0:9)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[c(1:7,9,10),c(3:13)]

subset <- unique(Dog.combined.singlet.B$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet.B$Subset2) # PBMC TIL

SubsetColors2 <- c("darkred", "#0348A6") 
names(SubsetColors2) <- c("PBMC", "TIL") 
clusterColors <- c(scales::hue_pal()(10))
names(clusterColors) <- 0:9

colors <- list(Subset2 = SubsetColors2, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:9) #rownames(heatmap2) <- paste0("X", 1:22)
rownames(heatmap2) <- paste0("C", c(0:6,8,9)) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/Featureplot_Cluster_Marker/B/GSEA_B_Dog_Select2.pdf", width = 7, height =5) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 12, 
                   cellheight = 12,
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   ))
dev.off()

#
DimPlot(Dog.combined.singlet.T, reduction = "tsne", label = T)

meta <- Dog.combined.singlet.T[[]]
meta$seurat_clusters <- Dog.combined.singlet.T@active.ident
ncol(meta)
View(meta)
meta <- meta[,c(5,16,21:38,131,132)]

heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:20, mean)) 
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C", 0:32)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:22)]

subset <- unique(Dog.combined.singlet.T$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet.T$Subset2) # PBMC TIL

unique(Dog.combined.singlet.T$seurat_clusters)
SubsetColors2 <- c("darkred", "darkred") 
names(SubsetColors2) <- c("PBMC", "PBMC") 
clusterColors <- c(scales::hue_pal()(33))
names(clusterColors) <- 0:32

colors <- list(Subset2 = SubsetColors2, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:32) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/Featureplot_Cluster_Marker/GSEA_T_Dog.pdf", width = 12, height =7) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 12, 
                   cellheight = 12,
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   ))
dev.off()
VlnPlot(Dog.combined.singlet.CD4, features = c("CTLA4"))

library(RColorBrewer)
FeaturePlot(Dog.combined.singlet.T, features = c("PDCD1","CD8A","CD4"), label = T, reduction = "tsne", order = T)

FeaturePlot(Dog.combined.singlet.B, features = c("Plasma_Cell_Differentiation"),
            order = T,
            cols = viridis(100)) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Feature_Plasma_Cell_Differentiation.png", width = 5, height = 5)
#
DefaultAssay(Dog.combined.singlet.B) <- "RNA"
Dog.combined.singlet.B <-ScaleData(Dog.combined.singlet.B, verbose = FALSE, features = rownames(Dog.combined.singlet.B)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet.B <- JoinLayers(Dog.combined.singlet.B)
Dog.markers <- FindAllMarkers(Dog.combined.singlet.B, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top10, "Data_Healthy/Featureplot_Cluster_Marker/B/B_Top10_0.5_0.5.csv",append = "F",sep = ",")
ggsave("Data_Healthy/Tsne_Labeled_Screening_3.5.png", width=7, height=5)

DoHeatmap(subset(Dog.combined.singlet.B, downsample = 100),
          features = top10$gene, 
          disp.max = 2,
          disp.min = -1,
          size = 3, angle = 50, cells = 1:50)+
  theme(text = element_text(size=8)) + 
  scale_fill_viridis(option = "E") + 
  labs(x = "Viridis E", y = NULL) #+ NoLegend()

#+ scale_fill_viridis(100)
#+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Heatmap_top10.png", width=8, height=7)

DimPlot(Dog.combined.singlet.B, reduction = "umap", label = T, label.size = 6) #+ NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/B_UMAP_legend.png", width = 5, height = 5)

Dog.combined.singlet.B@active.ident <- Dog.combined.singlet.B$seurat_clusters
dotPlot(Dog.combined.singlet.B, features = c("CD3E","DPYD","CD79B","MS4A1",
                                             "HLA-DRB1","DLA-DRA","BTG1","ID3",
                                             "LOC100685971","VPREB3","NIBAN3","MX1","XAF1","TCF7","GNG2",
                                             
                                             "JCHAIN",
                                             "MZB1","XBP1","TNFRSF17","PRDM1","IRF4",
                                             "CD38","EZH2","UBE2C","NCAPG","MKI67","OSBPL10"),
        cluster.idents = T, 
        colormap = "Blues") 
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Dot.png", width = 5, height = 6)

cluster9 <- FindMarkers(Dog.combined.singlet.B,
                        ident.1 = c(5), #only.pos = T,
                        #  ident.2 = c(16,5,9,8,12,18),
                        min.pct = 0.3)
write.table(cluster9, "Data_Healthy/Featureplot_Cluster_Marker/B/Cluster5_DEG.csv", append = F, sep = ",")

naivenmarkers <- FindMarkers(Dog.combined.singlet.B, ident.1 = c(5), min.pct = 0.36, only.pos = F, logfc.threshold = 0)
write.table(naivenmarkers, "Data_Healthy/Featureplot_Cluster_Marker/B/plasma5_DEG.csv", append = F, sep = ",")

head(naivenmarkers, 5)
FeaturePlot(Dog.combined.singlet.B, features = c("LOC102152513"), #  
            reduction = "umap", label = F, pt.size = 1.5, 
            order = T, repel = T,#pt.size = 0.8,
            cols = c("#ADDAE6","#E63222")
) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Feature_MZB1.png", width = 5, height = 5)

VlnPlot(Dog.combined.singlet.B, features = updwgene$updwgene, flip = T, stack = T, idents = c(4,5,8))

plot_density(Dog.combined.singlet, 
             features = c("LOC490151"),pal = "plasma",# "cividis", "inferno", "plasma"
             reduction = "tsne") +NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/Density_CD4.png",
       width = 4, height = 4)

#

unique(Dog.combined.doublet$seurat_clusters)
Dog.combined.doublet <- AddMetaData(Dog.combined.doublet, ES.immune) # add metadata
meta <- Dog.combined.doublet[[]] 
View(meta)
meta <- meta[,c(4,16,23:130)] 
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:108, mean))
View(heatmap)
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C",0:15) 
headers <- heatmap[,c(1:2)] 
heatmap <- heatmap[,c(3:110)]

SubsetColors <- c("darkred") 
names(SubsetColors) <- c("Healthy") 
clusterColors <- c(scales::hue_pal()(16))
names(clusterColors) <- 0:15

colors <- list(Subset = SubsetColors, seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:15) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/GSEA_CD4_Sorted.pdf", width = 10, height =7) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols=T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   #   cellheight = 12,
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()

#
cc.genes <- Seurat::cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Dog.combined.singlet.B <- CellCycleScoring(Dog.combined.singlet.B, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

freq_table <- Dog.combined.singlet.B[[]]
freq_table <- freq_table[,c("Subset", "seurat_clusters", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(0,1,2,3,4,5,6,8,9))

freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(0,1,2,3,4,5,6,8,9))
ggplot(freq_table, aes(x=seurat_clusters, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() +
  geom_text(aes(label = percent),
            position = position_stack(vjust = .5), size = 3)
ggsave("PATH_TO_FILE", height=5, width=15)
ggsave("PATH_TO_FILE", height=6, width=20)
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/Cluster_9_12_cellcycle.png", height=2, width=3)

FeatureScatter(Dog.combined.singlet.CD4,jitter = F, log = T,
               feature1 = "T1_Interferon" , feature2 = "Dog_MMVD")
DimPlot(Dog.combined.singlet.CD4)
meta <- Dog.combined.singlet.CD4[[]]
View(meta)

#
freq_table <- Dog.combined.singlet.B[[]]
freq_table <- freq_table[,c("Subset", "old.ident", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, old.ident, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(0,1,2,3,4,5,6,8,9))

freq_table <- freq_table %>%
  group_by(Subset, old.ident) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))

ggplot(freq_table, aes(x=old.ident, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() +
  geom_text(aes(label = percent),
            position = position_stack(vjust = .5), size = 3)
ggsave("Data_Healthy/Featureplot_Cluster_Marker/B/B_cellcycle.png", height=2.75, width=5)

#


# T DN
cc.genes <- Seurat::cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Dog.combined.singlet <- CellCycleScoring(Dog.combined.singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

freq_table <- Dog.combined.singlet[[]]
freq_table <- freq_table[,c("Subset", "seurat_clusters", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50))

freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50))
ggplot(freq_table, aes(x=seurat_clusters, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() +
  geom_text(aes(label = percent),
            position = position_stack(vjust = .5), size = 4.5)
ggsave("PATH_TO_FILE", height=5, width=15)
ggsave("PATH_TO_FILE", height=6, width=20)

#
freq_table <- Dog.combined[[]]
freq_table <- freq_table[,c("Subset", "old.ident", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, old.ident, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$old.ident <- factor(freq_table$old.ident, levels = c("T_CD4","T_CD8"))

freq_table <- freq_table %>%
  group_by(Subset, old.ident) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$old.ident <- factor(freq_table$old.ident, levels = c("T_CD4","T_CD8"))
ggplot(freq_table, aes(x=old.ident, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() +
  geom_text(aes(label = percent),
            position = position_stack(vjust = .5), size = 3)
ggsave("Data_Healthy/cellcycle.png", height=3, width=15)

#
hpca.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Hpca, labels = Hpca$label.fine, assay.type.test=1)
monaco.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Monaco, labels = Monaco$label.fine, assay.type.test=1)
immune.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Dice, labels = Dice$label.fine, assay.type.scale=1)
immgen.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Immgen, labels = Immgen$label.fine, assay.type.test=1)
blue.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Blue, labels = Blue$label.fine, assay.type.test=1)
nover.CD4 <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD4), ref = Nover, labels = Nover$label.fine, assay.type.test=1)

hpca.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Hpca, labels = Hpca$label.fine, assay.type.test=1)
monaco.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Monaco, labels = Monaco$label.fine, assay.type.test=1)
immune.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Dice, labels = Dice$label.fine, assay.type.test=1)
immgen.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Immgen, labels = Immgen$label.fine, assay.type.test=1)
blue.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Blue, labels = Blue$label.fine, assay.type.test=1)
nover.CD8NKT <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.CD8NKT), ref = Nover, labels = Nover$label.fine, assay.type.test=1)

hpca.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Hpca, labels = Hpca$label.fine, assay.type.test=1)
monaco.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Monaco, labels = Monaco$label.fine, assay.type.test=1)
immune.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Dice, labels = Dice$label.fine, assay.type.test=1)
immgen.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Immgen, labels = Immgen$label.fine, assay.type.test=1)
blue.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Blue, labels = Blue$label.fine, assay.type.test=1)
nover.B <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.B), ref = Nover, labels = Nover$label.fine, assay.type.test=1)

saveRDS(hpca.B, "hpca.Dog.combined.singlet.B.rds")
saveRDS(monaco.B, "monaco.Dog.combined.singlet.B.rds")
saveRDS(immune.B, "immune.Dog.combined.singlet.B.rds")
saveRDS(immgen.B, "immgen.Dog.combined.singlet.B.rds")
saveRDS(nover.B, "nover.Dog.combined.singlet.B.rds")
saveRDS(blue.B, "blue.Dog.combined.singlet.B.rds")

hpca.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Hpca, labels = Hpca$label.fine, assay.type.test=1)
monaco.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Monaco, labels = Monaco$label.fine, assay.type.test=1)
immune.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Dice, labels = Dice$label.fine, assay.type.test=1)
immgen.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Immgen, labels = Immgen$label.fine, assay.type.test=1)
blue.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Blue, labels = Blue$label.fine, assay.type.test=1)
nover.myeloid <- SingleR(test =as.SingleCellExperiment(Dog.combined.singlet.myeloid), ref = Nover, labels = Nover$label.fine, assay.type.test=1)

saveRDS(hpca.myeloid, "hpca.Dog.combined.singlet.myeloid.rds")
saveRDS(monaco.myeloid, "monaco.Dog.combined.singlet.myeloid.rds")
saveRDS(immune.myeloid, "immune.Dog.combined.singlet.myeloid.rds")
saveRDS(immgen.myeloid, "immgen.Dog.combined.singlet.myeloid.rds")
saveRDS(nover.myeloid, "nover.Dog.combined.singlet.myeloid.rds")
saveRDS(blue.myeloid, "blue.Dog.combined.singlet.myeloid.rds")
#
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = hpca.CD4$labels, col.name = "hpca")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = monaco.CD4$labels, col.name = "monaco")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = immune.CD4$labels, col.name = "immune")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = immgen.CD4$labels, col.name = "immgen")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = blue.CD4$labels, col.name = "blue")
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, metadata = nover.CD4$labels, col.name = "nover")

Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = hpca.B$labels, col.name = "hpca")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = monaco.B$labels, col.name = "monaco")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = immune.B$labels, col.name = "immune")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = immgen.B$labels, col.name = "immgen")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = blue.B$labels, col.name = "blue")
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, metadata = nover.B$labels, col.name = "nover")

Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = hpca.CD8NKT$labels, col.name = "hpca")
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = monaco.CD8NKT$labels, col.name = "monaco")
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = immune.CD8NKT$labels, col.name = "immune")
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = immgen.CD8NKT$labels, col.name = "immgen")
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = blue.CD8NKT$labels, col.name = "blue")
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, metadata = nover.CD8NKT$labels, col.name = "nover")

Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = hpca.myeloid$labels, col.name = "hpca")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = monaco.myeloid$labels, col.name = "monaco")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = immune.myeloid$labels, col.name = "immune")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = immgen.myeloid$labels, col.name = "immgen")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = blue.myeloid$labels, col.name = "blue")
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, metadata = nover.myeloid$labels, col.name = "nover")

#
library(RColorBrewer)
DimPlot(Dog.combined.singlet.myeloid,
        reduction = "umap", group.by  = "immgen", 
        label = T, shuffle = T, repel = T, label.box = T, 
        label.size = 6, #pt.size = 1,
        cols = colorRampPalette(brewer.pal(12, "Paired"))(36))+ NoLegend()
ggsave("SingleR_CD4_monaco_No_Legend.png", width = 9, height = 9)



# canfam6
# setwd("PATH_TO_PROJECT")
H1.data <- Read10X("Gyeonmdong/")
H2.data <- Read10X("Makdoong/")
H3.data <- Read10X("Batt/")
H4.data <- Read10X("Koko//")
H5.data <- Read10X("Sojin//")
H6.data <- Read10X("Bori//")
H1 <- CreateSeuratObject(counts = H1.data, assay = "RNA", project = "Healthy1")
H2 <- CreateSeuratObject(counts = H2.data, assay = "RNA", project = "Healthy2")
H3 <- CreateSeuratObject(counts = H3.data, assay = "RNA", project = "Healthy3")
H4 <- CreateSeuratObject(counts = H4.data, assay = "RNA", project = "Healthy4")
H5 <- CreateSeuratObject(counts = H5.data, assay = "RNA", project = "Healthy5")
H6 <- CreateSeuratObject(counts = H6.data, assay = "RNA", project = "Healthy6")
H1$Subset <- "Healthy"
H2$Subset <- "Healthy"
H3$Subset <- "Healthy"
H4$Subset <- "Healthy"
H5$Subset <- "Healthy"
H6$Subset <- "Healthy"
H1$Subset2 <- "PBMC" #PBMC TIL
H2$Subset2 <- "PBMC" 
H3$Subset2 <- "PBMC"
H4$Subset2 <- "PBMC"
H5$Subset2 <- "PBMC"
H6$Subset2 <- "PBMC"
H1$Subset3 <- "Healthy" # Tumor type
H2$Subset3 <- "Healthy" 
H3$Subset3 <- "Healthy"
H4$Subset3 <- "Healthy"
H5$Subset3 <- "Healthy"
H6$Subset3 <- "Healthy"
H1$Subset4 <- "Healthy_PBMC" # Tumor type
H2$Subset4 <- "Healthy_PBMC" 
H3$Subset4 <- "Healthy_PBMC"
H4$Subset4 <- "Healthy_PBMC"
H5$Subset4 <- "Healthy_PBMC"
H6$Subset4 <- "Healthy_PBMC"

H1$Subset5 <- "11y"
H2$Subset5 <- "8y"
H3$Subset5 <- "7y"
H4$Subset5 <- "10y"
H5$Subset5 <- "12y"
H6$Subset5 <- "8y"

H1.data <- Read10X("Gyeonmdong/")
H2.data <- Read10X("Makdoong/")
H3.data <- Read10X("Batt/")
H4.data <- Read10X("Koko//")
H5.data <- Read10X("Sojin//")
H6.data <- Read10X("Bori//")

H1$Subset6 <- "Mixed"
H2$Subset6 <- "Maltese"
H3$Subset6 <- "Mixed"
H4$Subset6 <- "Mixed"
H5$Subset6 <- "Poodle"
H6$Subset6 <- "Maltese"

H1$Subset7 <- "IM"
H2$Subset7 <- "CM"
H3$Subset7 <- "SF"
H4$Subset7 <- "CM"
H5$Subset7 <- "CM"
H6$Subset7 <- "CM"

H1$Subset8 <- "Gyeomdong"
H2$Subset8 <- "Makdoong"
H3$Subset8 <- "Batt"
H4$Subset8 <- "Koko"
H5$Subset8 <- "Sojin"
H6$Subset8 <- "Bori"

Dog.BIG <-merge(H1, y = c(H2,H3,H4,H5,H6), 
                add.cell.ids = c("H1","H2","H3","H4","H5","H6"), 
                project = "Healthy_PBMC")

Dog.BIG[["percent.mt"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^MT-")
Dog.BIG[["percent.loc"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^LOC")
Dog.BIG[["percent.ribo"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^RP[LS]")

VlnPlot(Dog.BIG, features = c("nFeature_RNA","percent.loc","percent.ribo","percent.mt"), 
        group.by = "orig.ident", ncol = 4, pt.size = 0, raster=F)
VlnPlot(Dog.BIG, features = c("percent.mt"), pt.size = 0, group.by = "Subset4")
VlnPlot(Dog.BIG, features = c("nFeature_RNA"), pt.size = 0, group.by = "Subset")
saveRDS(Dog.BIG, "Dog.BIG.CanFam6.rds")

# 
Dog.BIG <-subset(Dog.BIG, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 & percent.mt <10) # 4000 is the standard
Dog_list <- SplitObject(Dog.BIG, split.by = "orig.ident")
Dog_list <- lapply(X = Dog_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)})
features <- SelectIntegrationFeatures(object.list = Dog_list)
Dog.anchors <- FindIntegrationAnchors(object.list = Dog_list, anchor.features = features, )
Dog.combined <- IntegrateData(anchorset = Dog.anchors) 

#
DefaultAssay(Dog.combined) <- "integrated"
Dog.combined <-ScaleData(Dog.combined, verbose = FALSE)
Dog.combined <-RunPCA(Dog.combined, npcs = 30, verbose = FALSE)
Dog.combined <-RunUMAP(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <-RunTSNE(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <- FindNeighbors(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <- FindClusters(Dog.combined, resolution = 1.2)
DefaultAssay(Dog.combined) <- "RNA"
saveRDS(Dog.combined, "Dog.combined.rds")
DimPlot(Dog.combined, reduction = "tsne")

#
Dog.combined <- JoinLayers(Dog.combined)
sce <- as.SingleCellExperiment(Dog.combined)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined <- AddMetaData(Dog.combined, doublets) 

#
Dog.combined.singlet <- subset(Dog.combined, subset = db.class == c("singlet")) #Dog.combined$db.class = 
DefaultAssay(Dog.combined.singlet) <- "integrated"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <-RunPCA(Dog.combined.singlet, npcs = 30, verbose = FALSE)
Dog.combined.singlet <-RunUMAP(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- RunTSNE(Dog.combined.singlet, dims = 1:30, reduction = "pca")
Dog.combined.singlet <- FindNeighbors(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- FindClusters(Dog.combined.singlet, resolution = 3.7)
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <- JoinLayers(Dog.combined.singlet)
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) + NoLegend()

saveRDS(Dog.combined.singlet, "Dog.combined.singlet.rds")
FeaturePlot(Dog.combined.singlet, features = c("CD1A"), reduction = "tsne") 

write.table(top20, "top20.csv", append = F, sep = ",")

# canfam4
library(Seurat)
# setwd("PATH_TO_PROJECT")
H1.data <- Read10X("Gyeomdong_Healthy3/")
H2.data <- Read10X("Makdoong_Healthy4//")
H3.data <- Read10X("Batt_Healthy6//")
H4.data <- Read10X("Koko_Healthy1//")
H5.data <- Read10X("Sojin_Healthy5//")
H6.data <- Read10X("Bori_Healthy2//")
H1 <- CreateSeuratObject(counts = H1.data, assay = "RNA", project = "Healthy1")
H2 <- CreateSeuratObject(counts = H2.data, assay = "RNA", project = "Healthy2")
H3 <- CreateSeuratObject(counts = H3.data, assay = "RNA", project = "Healthy3")
H4 <- CreateSeuratObject(counts = H4.data, assay = "RNA", project = "Healthy4")
H5 <- CreateSeuratObject(counts = H5.data, assay = "RNA", project = "Healthy5")
H6 <- CreateSeuratObject(counts = H6.data, assay = "RNA", project = "Healthy6")
H1$Subset <- "Healthy"
H2$Subset <- "Healthy"
H3$Subset <- "Healthy"
H4$Subset <- "Healthy"
H5$Subset <- "Healthy"
H6$Subset <- "Healthy"
H1$Subset2 <- "PBMC" #PBMC TIL
H2$Subset2 <- "PBMC" 
H3$Subset2 <- "PBMC"
H4$Subset2 <- "PBMC"
H5$Subset2 <- "PBMC"
H6$Subset2 <- "PBMC"
H1$Subset3 <- "Healthy" # Tumor type
H2$Subset3 <- "Healthy" 
H3$Subset3 <- "Healthy"
H4$Subset3 <- "Healthy"
H5$Subset3 <- "Healthy"
H6$Subset3 <- "Healthy"
H1$Subset4 <- "Healthy_PBMC" # Tumor type
H2$Subset4 <- "Healthy_PBMC" 
H3$Subset4 <- "Healthy_PBMC"
H4$Subset4 <- "Healthy_PBMC"
H5$Subset4 <- "Healthy_PBMC"
H6$Subset4 <- "Healthy_PBMC"

H1$Subset5 <- "11y"
H2$Subset5 <- "8y"
H3$Subset5 <- "7y"
H4$Subset5 <- "10y"
H5$Subset5 <- "12y"
H6$Subset5 <- "8y"

H1$Subset6 <- "Mixed"
H2$Subset6 <- "Maltese"
H3$Subset6 <- "Mixed"
H4$Subset6 <- "Mixed"
H5$Subset6 <- "Poodle"
H6$Subset6 <- "Maltese"

H1$Subset7 <- "IM"
H2$Subset7 <- "CM"
H3$Subset7 <- "SF"
H4$Subset7 <- "CM"
H5$Subset7 <- "CM"
H6$Subset7 <- "CM"

H1$Subset8 <- "Gyeomdong"
H2$Subset8 <- "Makdoong"
H3$Subset8 <- "Batt"
H4$Subset8 <- "Koko"
H5$Subset8 <- "Sojin"
H6$Subset8 <- "Bori"

Dog.BIG <-merge(H1, y = c(H2,H3,H4,H5,H6), 
                add.cell.ids = c("H1","H2","H3","H4","H5","H6"), 
                project = "Healthy_PBMC")

Dog.BIG[["percent.mt"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^MT-")
Dog.BIG[["percent.loc"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^LOC")
Dog.BIG[["percent.ribo"]] <- PercentageFeatureSet(Dog.BIG, pattern = "^RP[LS]")

VlnPlot(Dog.BIG, features = c("nFeature_RNA","percent.loc","percent.ribo","percent.mt"), 
        group.by = "orig.ident", ncol = 4, pt.size = 0, raster=F)
VlnPlot(Dog.BIG, features = c("percent.mt"), pt.size = 0, group.by = "Subset4")
VlnPlot(Dog.BIG.4, features = c("percent.loc"), pt.size = 0) +NoLegend()
saveRDS(Dog.BIG, "Dog.BIG.CanFam4.rds")

mean(Dog.BIG$percent.mt*100)
# 
Dog.BIG <-subset(Dog.BIG, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 & percent.mt <10) # 4000 is the standard
Dog_list <- SplitObject(Dog.BIG, split.by = "orig.ident")
Dog_list <- lapply(X = Dog_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)})
features <- SelectIntegrationFeatures(object.list = Dog_list)
Dog.anchors <- FindIntegrationAnchors(object.list = Dog_list, anchor.features = features)
Dog.combined <- IntegrateData(anchorset = Dog.anchors) 

variablegenes <- FindVariableFeatures(Dog.combined, nfeatures = 3000)
#
DefaultAssay(Dog.combined) <- "integrated"
Dog.combined <-ScaleData(Dog.combined, verbose = FALSE)
Dog.combined <-RunPCA(Dog.combined, npcs = 30, verbose = FALSE)
Dog.combined <-RunUMAP(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <-RunTSNE(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <- FindNeighbors(Dog.combined, reduction = "pca", dims = 1:30)
Dog.combined <- FindClusters(Dog.combined, resolution = 1.2)
DefaultAssay(Dog.combined) <- "RNA"
#saveRDS(Dog.combined, "PATH_TO_FILE")
DimPlot(Dog.combined, reduction = "tsne")

#
Dog.combined <- JoinLayers(Dog.combined)
sce <- as.SingleCellExperiment(Dog.combined)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined <- AddMetaData(Dog.combined, doublets) 

#
Dog.combined.singlet <- subset(Dog.combined, subset = db.class == c("singlet")) #Dog.combined$db.class = 
DefaultAssay(Dog.combined.singlet) <- "integrated"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <-RunPCA(Dog.combined.singlet, npcs = 30, verbose = FALSE)
Dog.combined.singlet <-RunUMAP(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- RunTSNE(Dog.combined.singlet, dims = 1:30, reduction = "pca")
Dog.combined.singlet <- FindNeighbors(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- FindClusters(Dog.combined.singlet, resolution = 3.7)
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <- JoinLayers(Dog.combined.singlet)
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) 

#saveRDS(Dog.combined.singlet, "Dog.combined.singlet.rds")

FeaturePlot(Dog.combined.singlet, order=T,features = c("KLRF1","CD160","RHEX",
                                                       "GZMB"), reduction = "tsne")
VlnPlot(Dog.combined.singlet, features = c("nFeature_RNA"))

DefaultAssay(Dog.combined) <- "RNA"
DimPlot(Dog.combined, reduction = "tsne", group.by = "db.class", label = T) 
saveRDS(Dog.combined, "Dog.combined.rds")

# analysis of CanFam4
# setwd("PATH_TO_PROJECT")
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) 
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T, label.size = 3)+NoLegend()
ggsave("tSNE_nolabel.png", width = 4, height = 4)

DefaultAssay(Dog.combined) <- "RNA"
Dog.combined <- JoinLayers(Dog.combined)
FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("CTLA4"),
            reduction = "umap", label = T, #pt.size = 0.5, 
            order = T, repel = T, pt.size = 0.8,
            cols = c("#ADDAE6","#E63222")
            #cols = c("lightgrey", "darkred"), 
            #cols = c("lightgrey", "darkgreen")
) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/ITGAX.png",
       width = 4, height = 4)

ggsave("PATH_TO_FILE",
       width = 4, height = 4)

FeaturePlot(Dog.combined.singlet.3, features = c("BATF"), 
            reduction = "tsne", label = F, pt.size = 0.8, order = T, repel = T, #pt.size = 0.8,
            cols = c("#ADDAE6","#E63222")
            #c("#ADDAE6","#E63222"), #cols = c("lightgrey", "darkred"), #cols = c("lightgrey", "darkgreen")
) + NoLegend()
ggsave("PATH_TO_FILE",
       width = 4, height = 4)
ggsave("PATH_TO_FILE",
       width = 4, height = 4)

FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("FOXP3","IL2RA"), order = T)

Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters
dotPlot(Dog.combined.singlet, 
        features = c("PTPRC","CD3D","CD3E","LOC608395","BATF",
                     # "LOC119875827",
                     "LOC607937",#"LOC106559087","LOC106559086","LOC111097255","LOC119874697","LOC609004",
                     "CD79B","CD4","CD8A","CD8B","NKG7","CD14","CCL23","LOC100686511",
                     "CCR7","SELL","IL2RA","IL7R","S100A9","CAECAM1","DPYD",
                     "IFNG","IL4","IL5","IL10","IL13","IL17A","IL17F","TNF","ICOS",
                     "PDCD1","CTLA4","HAVCR2","TIGIT","LAG3","CD274","TOX","TRGC2",
                     "ITGB1","CD28","FOXP3","TBX21","GATA3","RORC","STAT1","STAT4","IL12RB2",
                     "STAT6","STAT3","IL2","TNFRSF4","TNFRSF18","GZMK","GZMB","PRF1","KLRF1","LOC486692",
                     "CD160","CSF3R","CD1A","CD1C","CD1E","TCF4","IL3RA","MZB1","CD19","RHEX","LOC491694"),
        #angle.y = 180, 
        cluster.idents = T, angle.x = 0
        # idents = c(9,2,12,10,21,38,0,11,3,19,25,13,6,24,37,39,31,44,45,38)
)
ggsave("Data_Healthy/Dotplot_Cluster.png", 
       width = 13, height = 12)
library(Nebulosa)
dotPlot(Dog.combined.singlet, features = c("CD1D"))
VlnPlot(Dog.combined.singlet, features = c("CD1D"), sort = T, pt.size = 0)
FeaturePlot(Dog.combined.singlet.myeloid, order = T,#size = 0.5,
            features = c("CD1E","CD1B","CD1C","CD1D","LOC608848"), #pt.size = 1.5,# order = T,
            #             size = 0.5,
            #pal = "plasma",# "cividis", "inferno", "plasma"
            reduction = "umap") +NoLegend()
ggsave("PATH_TO_FILE",
       width = 4, height = 4) 

plot_density(Dog.combined.singlet.3,  
             features = c("MMP12"), size = 0.5, #pal = "plasma",# "cividis", "inferno", "plasma"
             reduction = "tsne") +NoLegend()
ggsave("PATH_TO_FILE",
       width = 4, height = 4)
VlnPlot(Dog.combined.singlet, features = c("nFeature_RNA"))
genes <- read.csv("PATH_TO_FILE", header = T)
dotPlot(Dog.combined.singlet.3, features = genes)

#
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T)
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet,
                                     '0' = "T_CD4", '1' = "T_CD8", '2' = "T_CD4",'3' = "T_CD4", '4' = "T_CD4", 
                                     '5' = "T_CD4", '6' = "T_CD4", '7' = "Myeloid", '8' = "T_CD4", '9'= "Myeloid", 
                                     '10' = "Myeloid", '11'= "Myeloid", '12' = "T_CD4", '13' = "B", '14' = "Myeloid", 
                                     '15' = "T_CD4", '16' = "T_CD8", '17' = "Myeloid", 
                                     '18' = "Myeloid", '19' = "Myeloid", '20' = "Myeloid", '21' = "T_CD4",
                                     '22' = "T_CD8", '23' = "Myeloid", '24' = "T_CD8", '25' = "Un", '26' = "B", 
                                     '27' = "T", '28' = "T_CD4",
                                     '29' = "T_CD4", '30' = "Myeloid", '31' = "T_CD4", '32' = "T_CD8", '33' = "Myeloid", 
                                     '34' = "Myeloid", '35' = "Myeloid", 
                                     '36' = "Myeloid", '37' = "B", '38' = "Myeloid", '39' = "T", '40' = "T_CD4",
                                     '41' = "B", '42' = "Myeloid", '43' = "T_CD4", '44' = "Myeloid", '45' = "Myeloid",
                                     '46' = "Myeloid", '47' = "Myeloid", '48' = 'T', '49'='Myeloid', '50' = 'Un')
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters 

# Final
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) + NoLegend()
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet,
                                     '0' = "T_CD4", '1' = "T_CD8", '2' = "T_CD4",'3' = "T_CD4", '4' = "T_CD4", 
                                     '5' = "T_CD4", '6' = "T_CD8", '7' = "Myeloid", '8' = "T_CD4", '9'= "Myeloid", 
                                     '10' = "Myeloid", '11'= "Myeloid", '12' = "T_CD4", '13' = "B", '14' = "Myeloid", 
                                     '15' = "T_CD4", '16' = "T_CD8", '17' = "Myeloid", 
                                     '18' = "Myeloid", '19' = "Myeloid", '20' = "Myeloid", '21' = "T_CD4",
                                     '22' = "T_CD8", '23' = "Myeloid", '24' = "T_CD8", '25' = "Un", '26' = "B", 
                                     '27' = "T", '28' = "T_CD4",
                                     '29' = "T_CD4", '30' = "Myeloid", '31' = "T_CD4", '32' = "T_CD8", '33' = "Myeloid", 
                                     '34' = "Myeloid", '35' = "Myeloid", 
                                     '36' = "Myeloid", '37' = "B", '38' = "Myeloid", '39' = "T", '40' = "T_CD4",
                                     '41' = "B", '42' = "Myeloid", '43' = "T", '44' = "Myeloid", '45' = "Myeloid",
                                     '46' = "Myeloid", '47' = "Myeloid", '48' = 'T', '49'='Myeloid', '50' = 'Un')
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters 
dotPlot(Dog.combined.singlet.3, features = c("CD3D","CD3E","CD4","CD8A"))
#
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE, features = rownames(Dog.combined.singlet)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet <- JoinLayers(Dog.combined.singlet)
Dog.markers <- FindAllMarkers(Dog.combined.singlet, only.pos = TRUE, min.pct = 0.36, logfc.threshold = 0.36) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top50, "Data_Healthy/Top50_0.36_0.36.csv", append = "F",sep = ",")
ggsave("Data_Healthy/Tsne_Labeled_Screening_3.5.png", width=7, height=5)

Dog.combined.singlet@active.ident <- factor(Dog.combined.singlet@active.ident, 
                                            levels = c("T_CD4", "T_CD8","B","Monocytes","DC","PMN","T_Gd","T","Un"))
DoHeatmap(subset(Dog.combined.singlet, downsample = 200),
          features = top5$gene, 
          #disp.max = 2,
          #  disp.min = -2,
          size = 3, angle = 50, cells = 1:200)+ NoLegend() #100 would be the best
theme(text = element_text(size=8)) #+NoLegend()
# scale_fill_viridis(option = "E")#+ 
# labs(x = "Viridis E", y = NULL) #+ NoLegend() #+ scale_fill_viridis(100) #+ NoLegend()
ggsave("PATH_TO_FILE", width=8, height=7)

DoHeatmap(Dog.combined.singlet,
          features = top5$gene, 
          size = 3, angle = 50)+ NoLegend() 
theme(text = element_text(size=8))
ggsave("PATH_TO_FILE", width=15, height=15)

Cellnumbers.subset <- table(Idents(Dog.combined.singlet), Dog.combined.singlet$orig.ident)
write.table(Cellnumbers.subset,
            "PATH_TO_FILE", sep = ",", append=F)

Cellnumbers.subset <- table(Idents(Dog.combined), Dog.combined$orig.ident)
write.table(Cellnumbers.subset,
            "PATH_TO_FILE", sep = ",", append=F)
DimPlot(Dog.combined.singlet, reduction = "tsne")

#
ES.immune <- escape.matrix(Dog.combined.singlet,
                           gene.sets = gene.set,
                           min.size = 3,
                           method = "ssGSEA")
saveRDS(ES.immune, "ES.immune.singlet.automatic_annotation.rds")
gene.set <- list(CD8nv = c("ITGA1", "LEF1", "PTGDR", "IL2RB", "ADGRG1", "NBEA"),
                 CD8_memory = c("GZMK", "GZMB", "PI3", "BTBD11", "CTSW", "CCR5", "CCL4", "KLRG1", "FASLG"),
                 CD8_gd = c("PTHLH", "IGF2BP2", "ABTB2", "AKAP12", "SOX4", "CTSW", "SLC16A10", "PXT1", "ZNRF3", "SULT2B1"),
                 CD4_naive = c("LEF1", "CSTA", "RGS10", "ZNF536", "CCR7", "COL6A5", "LTB", "TNFSF8"),
                 CD4_memory = c("LEF1", "TSHZ2", "CD52", "CCR7", "IL7R", "CTPS1", "EFHC2", "CARMIL1"),
                 CD4_effector_memory_Th1 = c("IL7R", "PTPN13", "IL18R1", "CD28", "RCAN2", "CCR9", "CCR5", "IL12RB2", "CD52", "PRUNE2"),
                 CD4_effector_memory_Th2 = c("RNF220", "ITGA2", "GATA3", "CCDC3", "LGALS3", "PTPN13", "S100A2", "PPEF1", "CMA1"), 
                 CD4_effector_memory_Th17 = c("NTRK2", "PTPN13", "ADAM12", "NRG2", "RGS17", "DNAH8", "CCR6", "NPAS2", "RORA", "LTBP1"), 
                 CD4_Treg = c("IKZF2", "CTLA4", "RGS1", "ICOS", "IL2RA", "CD28", "ZNF831"), 
                 CD4_IFN  = c("CXCL10", "IFI44", "OAS1", "ISG15", "IFI44L", "IFGGB2", "CTLA4", "STAT1", "DDX58", "XAF1"), 
                 Monocyte = c("LYZ", "BPI", "LRMDA", "MT2A", "F13A1", "FN1", "NRG1", "CCDC88A", "CD83", "RETN"), 
                 Monocyte_CD4 = c("IL1B", "MAFB", "NFKBIA", "CXCL8", "FN1", "BLOC1S6", "CD83", "S100P", "BPI", "NRG1"),
                 M_MDSC = c("IL18", "IL1B", "LTF", "MEFV", "KCNJ2", "CPXM2", "S100A12", "STEAP4", "CSF3R", "IL31RA"), 
                 Monocyte_IFN = c("RSAD2", "OAS1", "OAS2", "DDX58", "HERC6", "OAS3", "RTP4", "EIF2AK2", "IFIT2"), 
                 Pre_DC = c("FGF12", "GPHA2", "MTUS2", "FCER1A", "PLCE1", "PTPRS", "IGF1", "NECTIN1", "IL3RA", "AK8"), 
                 cDC1 = c("ZNF366", "SDC2", "DISC1", "ECRG4", "TMEM163", "RIMS2", "KIT", "OTOF", "RTKN2", "RAB7B"), 
                 cDC2 = c("PKIB", "CD300H", "SDC2", "CD1C", "NCAM2", "CD86", "BATF3", "ZNF366", "PID1", "ECM1"), 
                 pDC = c("COBLL1", "RAB3C", "IGF1", "FCER1A", "RYR1", "PRKG1", "CCND1", "STYXL2", "ANK1", "OCIAD2"), 
                 Un_DC = c("PLCB4", "ZNF366", "KCNK13", "STRIP2", "SDC2", "OTOF", "HACD1", "C5", "SLC8A1", "CNTLN"), 
                 Neutrophil = c("S100A12", "CD4", "SERPINA1", "SGK1", "S100A8", "ALDH1A2", "FNDC3B", "GGH", "SRGN", "IL1R2"), 
                 PMN_MDSC = c("CAMP", "PGLYRP1", "CRISP2", "MMP9", "MMP8", "TCN1", "CD177", "LTF", "FADS1", "S100A12"), 
                 Eosinophil = c("C30H15orf48", "TGM2", "DACH1", "PADI3", "SMPD1", "CA8", "IL5RA"),  
                 Basophil = c("DACH1", "CA8", "IL5RA", "DAPK2", "TGFA", "ANKRD33B", "HK2", "PRR5L"), 
                 B_immature = c("SYT1", "PAX5", "VPREB3", "ERC2", "TMTC2", "KLHL14", "F8", "TEX9", "TDRP", "ADGRF1"), 
                 B_naive = c("TNFRSF13C", "BANK1","HTR1F", "PAX5", "EBF1", "BTLA", "NRIP1", "ADAM9"), 
                 B_class_switched = c("TNFRSF13C", "GOLM1", "BANK1", "BTLA", "EBF1", "DYNC1I1", "MTMR2", "PAX5"),
                 B_activated = c("IGKC", "CACNB2", "PAX5", "TNFRSF13C", "IGHM", "RASGRF2", "AOX2", "BCAR3", "ADAM32"), 
                 Plasma = c("JCHAIN", "MZB1", "TXNDC5", "LMAN1", "FKBP11", "LAP3", "DERL3", "CCR10", "MKI67", "TNFRSF13B"), 
                 T_DN = c("KIAA0825", "TMEM132D", "KANK1", "NMB", "CTLA4", "SYNJ2", "BICDL1", "SLF1", "ID3", "KIAA1549"), 
                 T_gd = c("PARD3B", "RHEX", "IL17RB", "CDH4", "GATA3", "FAT1", "TOX2", "ADARB1", "ZNF683", "TGFBR3"), 
                 NK = c("KLRF1", "STMN2", "PAX4", "NCR3", "F2RL3", "CD96", "IL2RB", "IGSF3", "FREM1", "FASLG"), 
                 T_Cycle = c("TOP2A", "MKI67", "RRM2", "H1-5", "DIAPH3", "TK1", "KIF11", "TPX2", "ASPM"), 
                 NKT = c("GPA33", "TGFBR3", "KLRK1", "CD96", "SYTL2", "MOV10L1", "SLA2", "DSTN", "RARRES1"), 
                 CD34_unclassified = c("TFPI", "ZNF521", "CD34", "NDST3", "GUCY1A1", "HPGD", "CLEC3B", "KIT", "CD109", "DNTT"),
                 CD8_effector = c("GZMA", "GZMB", "CCL5", "CD7", "CD2", "ITGB7", "CD96", "CD3E"), 
                 T_naive = c("VIM", "IL7R", "S100A8", "EEF1A1", "CXCR4", "LTB", "TMSB10"), 
                 CD8_tissue_resident_memory = c("CCL4", "RGS1", "IFNG", "NR4A2", "BCL2A1", "FASLG"))

FeaturePlot(Dog.combined.singlet, features = c("CAECAM1"), reduction = "tsne", cols = plasma(100))

#
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, ES.immune) 
meta <- Dog.combined.singlet[[]]
meta$seurat_clusters <- Dog.combined.singlet@active.ident
ncol(meta)
View(meta)
meta <- meta[,c(5,16,21:57)] 
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:37, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("C", 0:50)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:39)]

levels(Dog.combined$seurat_clusters)
subset <- unique(Dog.combined.singlet$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet$Subset2) # PBMC TIL

SubsetColors2 <- c("darkred", "#0348A6") 
names(SubsetColors2) <- c("PBMC", "TIL") 
clusterColors <- c(scales::hue_pal()(51))
names(clusterColors) <- 0:50

colors <- list(Subset2 = SubsetColors2, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("C", 0:50) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("GSEA_Dog_Automatic_Annotation.pdf", width = 15, height =15) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 12, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()

FeaturePlot(Dog.combined.singlet, features = c("B_naive","B_activated"), reduction = "tsne", cols = viridis(100)) + NoLegend()
FeaturePlot(Dog.combined.singlet, features = c("B_naive","B_class_switched"), reduction = "tsne",blend = T,order = T) + NoLegend()
FeaturePlot(Dog.combined.singlet, features = c("B_activated","B_class_switched"), reduction = "tsne",blend = T,order = T) + NoLegend()
FeaturePlot(Dog.combined.singlet, features = c("B_naive","B_activated"), reduction = "tsne",blend = T,order = T, interactive = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 8, height = 4)

ggsave("PATH_TO_FILE", width = 4, height = 4)

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) + NoLegend()
ggsave("tsne.canfam3.png", width = 4, height = 4)

DimPlot(Dog.combined.singlet.3, reduction = "tsne", label = T) + NoLegend()
ggsave("tsne.canfam4.png", width = 4, height = 4)

#
Dog.combined.singlet <- AddMetaData(Dog.combined.singlet, ES.immune) 
meta <- Dog.combined.singlet[[]]
meta$seurat_clusters <- Dog.combined.singlet@active.ident
ncol(meta)
View(meta)
View(heatmap)
meta <- meta[,c(5,16,58:145)] 
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:88, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 0:50) #rownames(heatmap) <- paste0("X", 1:22)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:90)]

levels(Dog.combined$seurat_clusters)
subset <- unique(Dog.combined.singlet$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet$Subset2) # PBMC TIL

SubsetColors2 <- c("darkred", "#0348A6") 
names(SubsetColors2) <- c("PBMC", "TIL") 
clusterColors <- c(scales::hue_pal()(51))
names(clusterColors) <- 0:50

colors <- list(Subset2 = SubsetColors2, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 0:50) 

pdf("GSEA_Canfam4.pdf", width = 15, height =20) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 12, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()


#
meta <- Dog.combined.singlet.CD8NKT[[]]
meta$seurat_clusters <- Dog.combined.singlet.CD8NKT@active.ident
ncol(meta)
View(meta)
View(heatmap)
meta <- meta[,c(5,16,58:145)] 
heatmap <- meta %>% 
  group_by(seurat_clusters, Subset2) %>% 
  summarise(across(1:88, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 0:50) #rownames(heatmap) <- paste0("X", 1:22)
headers <- heatmap[,c(1,2)] 
heatmap <- heatmap[,c(3:90)]

levels(Dog.combined$seurat_clusters)
subset <- unique(Dog.combined.singlet$Subset) # Healthy Tumor
subset2 <- unique(Dog.combined.singlet$Subset2) # PBMC TIL

SubsetColors2 <- c("darkred", "#0348A6") 
names(SubsetColors2) <- c("PBMC", "TIL") 
clusterColors <- c(scales::hue_pal()(51))
names(clusterColors) <- 0:50

colors <- list(Subset2 = SubsetColors2, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 0:50) 

pdf("GSEA_Canfam4.pdf", width = 15, height =20) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "row", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 12, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   )
)
dev.off()

# subset for Canfam4
Dog.combined.singlet.myeloid <- subset(Dog.combined.singlet, 
                                       idents = c(38,10,33,35,23,18,19,42,34,17,30,11,7,20,9,36,14,46,47,45,44,49))
DefaultAssay(Dog.combined.singlet.myeloid) <- "integrated"
Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE)
Dog.combined.singlet.myeloid <-RunPCA(Dog.combined.singlet.myeloid, npcs = 30, verbose = FALSE)
Dog.combined.singlet.myeloid <-RunUMAP(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <-RunTSNE(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <- FindNeighbors(Dog.combined.singlet.myeloid, reduction = "pca", dims = 1:30)
Dog.combined.singlet.myeloid <- FindClusters(Dog.combined.singlet.myeloid, resolution = 1.2)
DefaultAssay(Dog.combined.singlet.myeloid) <- "RNA"
Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE)
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
saveRDS(Dog.combined.singlet.myeloid, "Dog.combined.singlet.myeloid.rds")
saveRDS(Dog.combined.singlet.CD4, "Dog.combined.singlet.CD4.rds")
saveRDS(Dog.combined.singlet.CD8NKT, "Dog.combined.singlet.CD8NKT.rds")

DimPlot(Dog.combined.singlet.myeloid, reduction = "umap", label = T)# + NoLegend()
DimPlot(Dog.combined.singlet.myeloid, reduction = "umap", label = T) + NoLegend()
ggsave("myeloid_UMA_LEGEND_label.png", width = 5, height = 5)
DimPlot(Dog.combined.singlet.myeloid, reduction = "tsne", label = T) + NoLegend()
ggsave("myeloid_TSNE_label.png", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.myeloid, features = c("MIF"),
            reduction = "umap", label = T, #pt.size = 0.5, 
            order = T, repel = T, pt.size = 0.8,
            cols = c("#ADDAE6","#E63222"))
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Myeloid_CCR1.png",
       width = 4, height = 4)

plot_density(Dog.combined.singlet.myeloid,
             features = c("ITGAE"), size = 0.5, pal = "plasma",# "cividis", ""plasma", "magma"
             reduction = "umap") +NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Density_Myeloid_ITGAX.png",
       width = 4, height = 4)
VlnPlot(Dog.combined.singlet.myeloid, features = "ITGAE", sort = T, pt.size = 0) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Vln_Myeloid_LGALS9.png",
       width = 5, height = 2)
dotPlot(Dog.combined.singlet.myeloid, features = c("CD3E","DPYD","LOC490269",
                                                   "S100A9","CD4","CD24","CD177","CEACAM1","LRMDA","KLF4","IER5L",
                                                   "CD14","CCL23","CITED4","IL10RB",
                                                   "CCL3","CCL4",
                                                   "DLA-DRA","DLA-DOB",
                                                   "ITGAX","CD86","CD83","CD80","DLA-DOA","PLD4",
                                                   "LOC608848","CD1B","CD1C","CD1D","CD1E",
                                                   "FCER1A","GPR183",
                                                   "LOC476816","LOC611536","IL1B",
                                                   "ISG15","MX1","IFI44","PGLYRP1",
                                                   "IL3RA","TCF4",
                                                   "XCR1","CCR1","CCR4","CCR5"), rotation = T,angle.x = 45,cluster.idents = T,
        
        colormap = "Blues",
        #idents = c(15,17,21,22),
        color.direction = 1)
ggsave("Data_Healthy/Dot_myeloid.png", width = 3, height = 2)
ggsave("Data_Healthy/Dot_myeloid.png", width = 8, height = 8)
library(CellChat)
dotPlot(Dog.combined.singlet.myeloid, features = c("ITGAX","CD86","CD83","DLA-DOA","PLD4",
                                                   "FCER1A","LOC100686511","GPR183",
                                                   "LOC608848","CD1B","CD1C","CD1D","CD1E",
                                                   "IL3RA","TCF4",
                                                   "CXCR3","XCR1"
), rotation = T,angle.x = 45,
colormap = "Blues", #idents = c(15,17,21,22),
color.direction = 1) + NoLegend()
ggsave("Data_Healthy/Dot_Myeloid_DC.png", width = 2, height = 3.5) #width = 2, height = 2.7
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters

Dog.combined.singlet.myeloid <-ScaleData(Dog.combined.singlet.myeloid, verbose = FALSE, features = rownames(Dog.combined.singlet.myeloid)) #https://github.com/satijalab/seurat/issues/6542
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
Dog.markers <- FindAllMarkers(Dog.combined.singlet.myeloid, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5) # label each clster and re-run this function and then you could get labled cluster's genes
Dog.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Dog.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
Dog.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
Dog.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) -> top15
Dog.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
Dog.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC) -> top6

Dog.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
Dog.markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7
Dog.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC) -> top8

write.table(top20, "Data_Healthy/Featureplot_Cluster_Marker/Top20_0.3_0.36.csv",append = "F",sep = ",")
ggsave("Data_Healthy/Tsne_Labeled_Screening_3.5.png", width=7, height=5)

DoHeatmap(subset(Dog.combined.singlet.myeloid, downsample = 100),
          features = top8$gene, 
          # disp.max = 2,
          #  disp.min = -2,X
          size = 3, angle = 50, cells = 1:150)+
  theme(text = element_text(size=8)) + 
  scale_fill_viridis(option = "E") + 
  labs(x = "Viridis E", y = NULL) #+ NoLegend()
ggsave("Data_Healthy/Heatmap_Myeloid_top8.png", width = 4.5, height = 3.2)
Dog.combined.singlet.myeloid <- RenameIdents(Dog.combined.singlet.myeloid,
                                             '0' = "PMN", '1' = "PMN", '2' = "PMN",'3' = "PMN", '4' = "Mo1", 
                                             '5' = "Mo", '6' = "PMN", '7' = "PMN", '8' = "Mo", '9'= "Mo", '10' = "Mo1",
                                             '11'= "PMN", '12' = "Mo", '13' = "Mo1", '14' = "Mo/T", '15' = "DC",
                                             '16' = "Mo", '17' = "PMN", 
                                             '18' = "Mo", '19' = "PMN", '20' = "PMN", '21' = "DC",
                                             '22' = "DC")
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters

Dog.combined.singlet.myeloid <- RenameIdents(Dog.combined.singlet.myeloid,
                                             '0' = "PMN", '1' = "PMN", '2' = "PMN",'3' = "PMN", '4' = "Mo", 
                                             '5' = "Mo", '6' = "PMN", '7' = "PMN", '8' = "Mo", '9'= "Mo", '10' = "Mo",
                                             '11'= "PMN", '12' = "Mo", '13' = "Mo", '14' = "Mo/T", '15' = "DC",
                                             '16' = "Mo", '17' = "PMN", 
                                             '18' = "Mo", '19' = "PMN", '20' = "PMN", '21' = "DC",
                                             '22' = "DC")
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters
#+ scale_fill_viridis(100)
#+ NoLegend()
ggsave("Data_Healthy/Heatmap_Labeled_Singlet_legend.png", width=8, height=7)
DimPlot(Dog.combined.singlet.myeloid)
FeaturePlot(Dog.combined.singlet.myeloid, features = c("ISG15"), order = T, pt.size = 1,
            cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Myeoid_CD274.png", width=4, height=5)

VlnPlot(Dog.combined.singlet.myeloid, features = c(""), #pt.size = 0, 
        # idents = c(2,3,11),
        split.by = "PDL1_Sep", cols = c("#ADDAE6","#E63222")) + 
  NoLegend() + 
  stat_compare_means("t.test")
ggsave("Data_Healthy/Myeloid/PDL1+PMN_CD47.png", width=4, height=5)



VlnPlot(Dog.combined.singlet.myeloid, features = c("CD47"), idents = c(10,4,13,5,9,16,12,18,8),
        group.by = "PDL1_Sep", cols = c("#ADDAE6","#E63222")) + NoLegend() #+ stat_compare_means("t.test")


DimPlot(Dog.combined.singlet.myeloid, split.by = "PDL1_Sep",
        group.by = "PDL1_Sep", cols = c("#ADDAE6","#E63222")) + NoLegend() #+ stat_compare_means("t.test")
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Myeoid_CD274_Sep.png", width=4, height=3)

Dog.combined.singlet.myeloid <- RunUMAP(Dog.combined.singlet.myeloid, dims = 1:30)
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
sce <- as.SingleCellExperiment(Dog.combined.singlet.myeloid)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.myeloid <- AddMetaData(Dog.combined.singlet.myeloid, doublets) 

Dog.combined.myeloid.singlet <- subset(Dog.combined.singlet.myeloid, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.myeloid.singlet) <- "integrated"
Dog.combined.myeloid.singlet <-ScaleData(Dog.combined.myeloid.singlet, verbose = FALSE)
Dog.combined.myeloid.singlet <-RunPCA(Dog.combined.myeloid.singlet, npcs = 30, verbose = FALSE)
Dog.combined.myeloid.singlet <-RunUMAP(Dog.combined.myeloid.singlet, reduction = "pca", dims = 1:30)
Dog.combined.myeloid.singlet <- RunTSNE(Dog.combined.myeloid.singlet, dims = 1:30, reduction = "pca")
Dog.combined.myeloid.singlet <- FindNeighbors(Dog.combined.myeloid.singlet, reduction = "pca", dims = 1:30)
Dog.combined.myeloid.singlet <- FindClusters(Dog.combined.myeloid.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.myeloid.singlet) <- "RNA"
Dog.combined.myeloid.singlet <-ScaleData(Dog.combined.myeloid.singlet, verbose = FALSE)
Dog.combined.myeloid.singlet <- JoinLayers(Dog.combined.myeloid.singlet)

DimPlot(Dog.combined.singlet.myeloid, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Myeloid_Doublet_tsne.png", width = 4, height = 4)
Myeloid.doublet.percentage <- table(Idents(Dog.combined.singlet.myeloid), Dog.combined.singlet.myeloid$db.class)
write.table(Myeloid.doublet.percentage, "Myeloid_db_%_.csv", sep = ",", append=F)

#
Dog.combined.singlet.B <- subset(Dog.combined.singlet, 
                                 idents = c(13,26,41,37,50))
DefaultAssay(Dog.combined.singlet.B) <- "integrated"
Dog.combined.singlet.B <-ScaleData(Dog.combined.singlet.B, verbose = FALSE)
Dog.combined.singlet.B <-RunPCA(Dog.combined.singlet.B, npcs = 30, verbose = FALSE)
Dog.combined.singlet.B <-RunUMAP(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <-RunTSNE(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <- FindNeighbors(Dog.combined.singlet.B, reduction = "pca", dims = 1:30)
Dog.combined.singlet.B <- FindClusters(Dog.combined.singlet.B, resolution = 1)
DefaultAssay(Dog.combined.singlet.B) <- "RNA"
Dog.combined.singlet.B <-ScaleData(Dog.combined.singlet.B, verbose = FALSE)
Dog.combined.singlet.B <- JoinLayers(Dog.combined.singlet.B)
saveRDS(Dog.combined.singlet.B, "Dog.combined.singlet.B.rds")

DimPlot(Dog.combined.singlet.B, reduction = "umap", label = T, label.color = 5) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)

DimPlot(Dog.combined.singlet.B, reduction = "tsne", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.B, features = c("JCHAIN","CD19","CD79B","MZB1"))

Dog.combined.singlet.B <- JoinLayers(Dog.combined.singlet.B)
sce <- as.SingleCellExperiment(Dog.combined.singlet.B)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.B <- AddMetaData(Dog.combined.singlet.B, doublets) 

Dog.combined.B.singlet <- subset(Dog.combined.singlet.B, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.B.singlet) <- "integrated"
Dog.combined.B.singlet <-ScaleData(Dog.combined.B.singlet, verbose = FALSE)
Dog.combined.B.singlet <-RunPCA(Dog.combined.B.singlet, npcs = 30, verbose = FALSE)
Dog.combined.B.singlet <-RunUMAP(Dog.combined.B.singlet, reduction = "pca", dims = 1:30)
Dog.combined.B.singlet <- RunTSNE(Dog.combined.B.singlet, dims = 1:30, reduction = "pca")
Dog.combined.B.singlet <- FindNeighbors(Dog.combined.B.singlet, reduction = "pca", dims = 1:30)
Dog.combined.B.singlet <- FindClusters(Dog.combined.B.singlet, resolution = 1)
DefaultAssay(Dog.combined.B.singlet) <- "RNA"
Dog.combined.B.singlet <-ScaleData(Dog.combined.B.singlet, verbose = FALSE)
Dog.combined.B.singlet <- JoinLayers(Dog.combined.B.singlet)

DimPlot(Dog.combined.singlet.B, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("B_Doublet_tsne.png", width = 4, height = 4)
B.doublet.percentage <- table(Idents(Dog.combined.singlet.B), Dog.combined.singlet.B$db.class)
write.table(B.doublet.percentage, "B_db_%_.csv", sep = ",", append=F)

#
Dog.combined$PTPRC <- Dog.combined[["RNA"]]$data["PTPRC",]
Dog.combined$Sep <- ifelse(Dog.combined$PTPRC > 0 | Dog.combined$PTPRC < 0 , "Pos","Neg")
DimPlot(Dog.combined, group.by = "Sep", reduction = "tsne")
FeaturePlot(Dog.combined, features = c("CD3E"), reduction = "tsne")

Dog.combined.singlet$PTPRC <- Dog.combined.singlet[["RNA"]]$data["PTPRC",]
Dog.combined.singlet$Sep <- ifelse(Dog.combined.singlet$PTPRC > 0 | Dog.combined.singlet$PTPRC < 0 , "Pos","Neg")
DimPlot(Dog.combined.singlet, group.by = "Sep", reduction = "tsne")
FeaturePlot(Dog.combined, features = c("CD3E"), reduction = "tsne")

#
JSS.combined.B.singlet$IRF7 <- JSS.combined.B.singlet[["RNA"]]@data[c("IRF7"),]

#
DimPlot(Dog.combined, reduction = "tsne", label = T)# +NoLegend()
FeaturePlot(Dog.combined.singlet, features = c("CD8A","CD4"), reduction = "tsne", order = T)


Dog.combined.singlet.CD4 <- subset(Dog.combined.singlet, idents = c(15,5,0,40,28,
                                                                    21,2,4,29,31,
                                                                    8,3,12)) 

DefaultAssay(Dog.combined.singlet.CD4) <- "integrated"
Dog.combined.singlet.CD4 <-ScaleData(Dog.combined.singlet.CD4, verbose = FALSE)
Dog.combined.singlet.CD4 <-RunPCA(Dog.combined.singlet.CD4, npcs = 30, verbose = FALSE)
Dog.combined.singlet.CD4 <-RunUMAP(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <-RunTSNE(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <- FindNeighbors(Dog.combined.singlet.CD4, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD4 <- FindClusters(Dog.combined.singlet.CD4, resolution = 1)
DefaultAssay(Dog.combined.singlet.CD4) <- "RNA"
Dog.combined.singlet.CD4 <-ScaleData(Dog.combined.singlet.CD4, verbose = FALSE)
Dog.combined.singlet.CD4 <-JoinLayers(Dog.combined.singlet.CD4)

DimPlot(Dog.combined.singlet.CD4, reduction = "umap", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)
DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)


FeaturePlot(Dog.combined.singlet.CD4, features = c("FOXP3","CD8A","CD4"), reduction = "tsne")
saveRDS(Dog.combined.singlet.CD4, "Dog.combined.singlet.CD4.rds")

#
Dog.combined.singlet.CD4 <- JoinLayers(Dog.combined.singlet.CD4)
sce <- as.SingleCellExperiment(Dog.combined.singlet.CD4)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, doublets) 

Dog.combined.CD4.singlet <- subset(Dog.combined.singlet.CD4, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.CD4.singlet) <- "integrated"
Dog.combined.CD4.singlet <-ScaleData(Dog.combined.CD4.singlet, verbose = FALSE)
Dog.combined.CD4.singlet <-RunPCA(Dog.combined.CD4.singlet, npcs = 30, verbose = FALSE)
Dog.combined.CD4.singlet <-RunUMAP(Dog.combined.CD4.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD4.singlet <- RunTSNE(Dog.combined.CD4.singlet, dims = 1:30, reduction = "pca")
Dog.combined.CD4.singlet <- FindNeighbors(Dog.combined.CD4.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD4.singlet <- FindClusters(Dog.combined.CD4.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.CD4.singlet) <- "RNA"
Dog.combined.CD4.singlet <-ScaleData(Dog.combined.CD4.singlet, verbose = FALSE)
Dog.combined.CD4.singlet <- JoinLayers(Dog.combined.CD4.singlet)


DimPlot(Dog.combined.singlet.CD4, reduction = "tsne", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("CD4_Doublet_tsne.png", width = 4, height = 4)

CD4.doublet.percentage <- table(Idents(Dog.combined.singlet.CD4), Dog.combined.singlet.CD4$db.class)
write.table(CD4.doublet.percentage, "CD4_db_%_.csv", sep = ",", append=F)

# All T
DefaultAssay(Dog.combined.singlet) <- "RNA"
Dog.combined.singlet.T <- subset(Dog.combined.singlet, idents = c(22,32,16,1,24,6,43,39,48,27,25,
                                                                  15,5,0,40,28,21,2,4,29,31,8,3,12))
DefaultAssay(Dog.combined.singlet.T) <- "integrated"
Dog.combined.singlet.T <-ScaleData(Dog.combined.singlet.T, verbose = FALSE)
Dog.combined.singlet.T <-RunPCA(Dog.combined.singlet.T, npcs = 30, verbose = FALSE)
Dog.combined.singlet.T <-RunUMAP(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <-RunTSNE(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <- FindNeighbors(Dog.combined.singlet.T, reduction = "pca", dims = 1:30)
Dog.combined.singlet.T <- FindClusters(Dog.combined.singlet.T, resolution = 3.7)
DefaultAssay(Dog.combined.singlet.T) <- "RNA"
Dog.combined.singlet.T <-ScaleData(Dog.combined.singlet.T, verbose = FALSE)
Dog.combined.singlet.T <- JoinLayers(Dog.combined.singlet.T)
saveRDS(Dog.combined.singlet.T, "PATH_TO_FILE")

DimPlot(Dog.combined.singlet.T, reduction = "tsne", label = T) + NoLegend()

ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)

FeaturePlot(Dog.combined.singlet.T, 
            features = c("CD3E","CD4","CD8A","FOXP3","CCR7","CTLA4"), 
            order = T, 
            reduction = "tsne", label = T)
FeaturePlot(Dog.combined.singlet.T, blend = T, features = c("CD4","CD8A"))
dotPlot(Dog.combined.singlet.T, features = c("CD4","CD8A","CD247"))
DimPlot(Dog.combined.singlet.T, reduction = "tsne", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)

Dog.combined

# CD8, NK, gd, T, unclassified
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T)
Dog.combined.singlet.CD8NKT <- subset(Dog.combined.singlet, idents = c(22,32,16,1,24,6,43,39,48,27,25))
DefaultAssay(Dog.combined.singlet.CD8NKT) <- "integrated"
Dog.combined.singlet.CD8NKT <-ScaleData(Dog.combined.singlet.CD8NKT, verbose = FALSE)
Dog.combined.singlet.CD8NKT <-RunPCA(Dog.combined.singlet.CD8NKT, npcs = 30, verbose = FALSE)
Dog.combined.singlet.CD8NKT <-RunUMAP(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <-RunTSNE(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <- FindNeighbors(Dog.combined.singlet.CD8NKT, reduction = "pca", dims = 1:30)
Dog.combined.singlet.CD8NKT <- FindClusters(Dog.combined.singlet.CD8NKT, resolution = 1.6)
DefaultAssay(Dog.combined.singlet.CD8NKT) <- "RNA"
Dog.combined.singlet.CD8NKT <-ScaleData(Dog.combined.singlet.CD8NKT, verbose = FALSE)
Dog.combined.singlet.CD8NKT <- JoinLayers(Dog.combined.singlet.CD8NKT)
saveRDS(Dog.combined.singlet.CD8NKT, "PATH_TO_FILE")

DimPlot(Dog.combined.singlet.CD8NKT, reduction = "umap", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)
DimPlot(Dog.combined.singlet.CD8NKT, reduction = "tsne", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)


FeaturePlot(Dog.combined.singlet.CD8NKT, 
            features = c("CD3E","CD4","CD8A","RHEX","XCL"), 
            order = F, 
            reduction = "umap", label = T)

#
Dog.combined.singlet.CD8NKT <- JoinLayers(Dog.combined.singlet.CD8NKT)
sce <- as.SingleCellExperiment(Dog.combined.singlet.CD8NKT)  
sce <- scDblFinder(sce)
doublets <- data.frame(db.weight.score = sce$scDblFinder.weighted,
                       db.class = sce$scDblFinder.class,
                       db.score = sce$scDblFinder.score)
rownames(doublets) <- rownames(sce@colData)
Dog.combined.singlet.CD8NKT <- AddMetaData(Dog.combined.singlet.CD8NKT, doublets) 

Dog.combined.CD8NKT.singlet <- subset(Dog.combined.singlet.CD8NKT, subset = db.class == c("singlet"))
DefaultAssay(Dog.combined.CD8NKT.singlet) <- "integrated"
Dog.combined.CD8NKT.singlet <-ScaleData(Dog.combined.CD8NKT.singlet, verbose = FALSE)
Dog.combined.CD8NKT.singlet <-RunPCA(Dog.combined.CD8NKT.singlet, npcs = 30, verbose = FALSE)
Dog.combined.CD8NKT.singlet <-RunUMAP(Dog.combined.CD8NKT.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD8NKT.singlet <- RunTSNE(Dog.combined.CD8NKT.singlet, dims = 1:30, reduction = "pca")
Dog.combined.CD8NKT.singlet <- FindNeighbors(Dog.combined.CD8NKT.singlet, reduction = "pca", dims = 1:30)
Dog.combined.CD8NKT.singlet <- FindClusters(Dog.combined.CD8NKT.singlet, resolution = 1.2)
DefaultAssay(Dog.combined.CD8NKT.singlet) <- "RNA"
Dog.combined.CD8NKT.singlet <-ScaleData(Dog.combined.CD8NKT.singlet, verbose = FALSE)
Dog.combined.CD8NKT.singlet <- JoinLayers(Dog.combined.CD8NKT.singlet)

DimPlot(Dog.combined.singlet.CD8NKT, reduction = "umap", pt.size = 0.1, label = F,order = T, group.by = "db.class",
        cols = c("#ADDAE6","#E63222")) + NoLegend()
ggsave("Myeloid_Doublet_tsne.png", width = 4, height = 4)
DimPlot(Dog.combined.singlet.CD8NKT, label = T, label.size = 5) #+ NoLegend()
ggsave("Data_Healthy/CD8NKT_Umap.png", width = 5, height = 5)
ggsave("Data_Healthy/CD8NKT_Umap_legend.png", width = 5, height = 5)

CD8NKT.doublet.percentage <- table(Idents(Dog.combined.singlet.CD8NKT), Dog.combined.singlet.CD8NKT$db.class)
write.table(CD8NKT.doublet.percentage, "CD8NKT_db_%_.csv", sep = ",", append=F)

#
DimPlot(Dog.combined.doublet, reduction = "umap", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)

DimPlot(Dog.combined.doublet, reduction = "tsne", label = T) + NoLegend()
ggsave("PATH_TO_FILE", width = 5, height = 5)
ggsave("PATH_TO_FILE", width = 5, height = 5)


FeaturePlot(Dog.combined.doublet, features = c("CD3D","CD5","DPYD","CD79B","CD14","CAECAM1"), order = T,
            label = F, repel = T, ncol = 5, reduction = "tsne")
ggsave("Data_Healthy/Doublet_TSNE.png", width = 15, height = 3)

#
sign <- read.delim("PATH_TO_FILE")

genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")

genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")
genes <-read.csv("PATH_TO_FILE")

Dog.combined.singlet.myeloid@active.ident<-Dog.combined.singlet.myeloid$seurat_clusters


Dog.combined.singlet.CD4@assays$RNA

counts <- GetAssayData(Dog.combined.singlet.CD4, assay = "RNA", slot = "data")
counts <- t(as.matrix(counts[rownames(counts) %in% genes$Genes,])) 
meta <- Dog.combined.singlet.CD4[[]][,c("Subset", "seurat_clusters")]

counts <- data.frame(meta, counts)
ncol(counts)
View(counts)
heatmap <- counts %>%
  group_by(seurat_clusters, Subset) %>% 
  summarise(across(1:4, mean)) #
heatmap <- as.data.frame(heatmap)

#need to make rownames to match the heatmap to the annotation
nrow(heatmap)
rownames(heatmap) <- paste0("X", 1:23) # rownames will be generated as X1, X2, ..., X35
headers <- heatmap[,1:2]
heatmap <- heatmap[,3:6] # All T

#Defining the color scheme
SubsetColors <-  c("#FF4B20")
names(SubsetColors) <- c("Healthy")
clusterColors <- c(scales::hue_pal()(23)) 
names(clusterColors) <- 0:22
colors <- list(Subset = SubsetColors, 
               seurat_clusters = clusterColors)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 1:23)

#Plotting non-scaled version
pdf("Data_Healthy/GSEA_Gene_Dog_IMHC.pdf", width = 10, height = 80)
pheatmap::pheatmap(t(heatmap2),show_colnames = T, scale = "row",
                   annotation_col = headers, cluster_cols= T, cluster_rows = T,  
                   color = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
                   annotation_colors = colors, fontsize = 12, 
                   cellwidth = 10, 
                   cellheight = 10,
                   legend = T, legend_labels = T, annotation_legend = T)
dev.off()

FeaturePlot(Dog.combined.singlet.myeloid, features = c("ZEB2"),
            
            order = T) + NoLegend()
VlnPlot(Dog.combined.singlet.myeloid, 
        features = c("MX1","IFI44","OAS1","CMPK2",
                     "ISG15","RTP4","CGAS","MX2"), 
        idents = c(10,4,13,14,16,5,9,8,18,12),
        stack = T, flip = T, sort = T)
cluster_marker_all <- FindMarkers(Dog.combined.singlet.myeloid, ident.1 = c(20),
                                  ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,
                                              11,12,13,14,15,16,17,18,19,21,22), 
                                  min.pct = 0.3)
cluster_marker_mo <- FindMarkers(Dog.combined.singlet.myeloid,
                                 ident.1 = c(18),
                                 ident.2 = c(10,4,13,14,16,5,9,8,12),
                                 min.pct = 0.3)
cluster_marker_mo <- FindMarkers(Dog.combined.singlet.myeloid,
                                 ident.1 = c(1,2,6,7,11),
                                 ident.2 = c(0,3,19,17,20), only.pos = T,
                                 min.pct = 0.3)
cluster_marker_mo <- FindMarkers(Dog.combined.singlet.myeloid,
                                 ident.1 = c(12,16,13),
                                 ident.2 = c(10,4,14,5,9,8,18),
                                 min.pct = 0.3)
cluster_marker_mo <- FindMarkers(Dog.combined.singlet.myeloid,
                                 ident.1 = c(4,10),only.pos = T,
                                 ident.2 = c(16,5,9,8,12,18),
                                 min.pct = 0.3)
cluster_marker_revision <- FindMarkers(Dog.combined.singlet.myeloid,
                                       ident.1 = c(5,8,9,12,18,16),only.pos = T,
                                       min.pct = 0.3)
write.table(cluster_marker_revision, "DEG_C589121618.csv", append = F, sep = ",")



head(cluster_marker_mo, 10)
dotPlot(Dog.combined.singlet.myeloid,
        features = c("MX1"), cluster.idents = T)

Dog.combined.singlet.myeloid <- RenameIdents(Dog.combined.singlet.myeloid,
                                             '0' = "PMN", '1' = "PMN", '2' = "PMN",'3' = "PMN", '4' = "Mo", 
                                             '5' = "Mo", '6' = "PMN", '7' = "PMN", '8' = "Mo", '9'= "Mo", '10' = "Mo",
                                             '11'= "PMN", '12' = "Mo", '13' = "Mo", '14' = "Mo", '15' = "DC1",
                                             '16' = "Mo", '17' = "PMN", 
                                             '18' = "Mo", '19' = "PMN", '20' = "PMN", '21' = "DC2",
                                             '22' = "DC3")
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters

PATTERN <- FindMarkers(Dog.combined.singlet.myeloid, ident.1 = "Mo_IFN", ident.2 = "PMN", only.pos = T)
write.table(PATTERN, "Data_Healthy/DEG_Mo_IFN_vs_PMN.csv", append = F, sep = ",")
grep("MIF",rownames(PATTERN), value = F)

Dog.combined.singlet.myeloid@active.ident <- factor(Dog.combined.singlet.myeloid@active.ident, 
                                                    levels = c("Mo","Mo_IFN","DC","PMN","PMN_CD177"))
VlnPlot(Dog.combined.singlet.myeloid, features = "MIF",
        cols = c("#C6FDEC", "#7AC5FF", "#FFB433", "#FF4B20", "#0348A6"),
        pt.size = 0.1) + NoLegend()# + stat_compare_means(method = "anova") 
ggsave("Data_Healthy/Myeloid_Vln_MIF.png", width = 3.2, height = 2.5)

FeaturePlot(Dog.combined.singlet.myeloid, features = c("LOC490269"), order = T)
Dog.combined.singlet.myeloid <- RenameIdents(Dog.combined.singlet.myeloid,
                                             '0' = "PMN", '1' = "PMN", '2' = "PMN",'3' = "PMN", '4' = "Mo", 
                                             '5' = "Mo", '6' = "PMN", '7' = "PMN", '8' = "Mo", '9'= "Mo", '10' = "Mo",
                                             '11'= "PMN", '12' = "Mo", '13' = "Mo", '14' = "Mo", '15' = "DC",
                                             '16' = "Mo", '17' = "PMN", 
                                             '18' = "C18", '19' = "PMN", '20' = "PMN", '21' = "DC",
                                             '22' = "DC")
VlnPlot(Dog.combined.singlet.myeloid, features = "LGALS9",
        cols = c("#7AC5FF", "#FFB433", "#FF4B20", "#0348A6"),
        pt.size = 0) + NoLegend()# + stat_compare_means(method = "anova") 
ggsave("Data_Healthy/Myeloid_Vln_LGALS9_Nodot.png", width = 3.2, height = 2.5)

FeaturePlot(Dog.combined.singlet.myeloid, features = c("LGALS9"))
VlnPlot(Dog.combined.singlet.myeloid, cols = c("#FF4B20","#0348A6","#0348A6",
                                               "#0348A6","#0348A6","#0348A6","#0348A6",
                                               "#0348A6","#0348A6","#0348A6","#0348A6"), # brewer.pal(10,"Spectral")
        features = c("MX1","IFI44","OAS1","CMPK2","LGALS9",
                     "ISG15","CGAS","MX2","ISG20","XAF1","IL18"), 
        stack = T, flip = T, sort = T)  + NoLegend() + stat_compare_means(method = "anova")
ggsave("Data_Healthy/Mo_IFN_CD177_genes.png", width=2.7, height=3.7)

DimPlot(Dog.combined.singlet.myeloid, reduction = "umap", label = T)
ggsave("Data_Healthy/Annotation/Umap_Myeloid.png", width = 5, height = 5)
# CD47, SLAMF7, PDCD1LG2
ggsave("Data_Healthy/Featureplot_Cluster_Marker/Myeoid_CD177.png", width=5, height=5)
Dog.combined.singlet.myeloid <- JoinLayers(Dog.combined.singlet.myeloid)
Dog.combined.singlet.myeloid <- runEscape(Dog.combined.singlet.myeloid, 
                                          method = "ssGSEA",
                                          gene.sets = list, 
                                          #  groups = 5000, 
                                          # min.size = 0,
                                          new.assay.name = "escape.ssGSEA")

Dog.combined.singlet.myeloid <- performNormalization(Dog.combined.singlet.myeloid, 
                                                     assay = "escape.ssGSEA", 
                                                     gene.sets = list)
rownames(Dog.combined.singlet.myeloid@assays$escape.ssGSEA@data)[1:100]
heatmapEnrichment(Dog.combined.singlet.myeloid, 
                  group.by = "seurat_clusters",
                  scale = T,
                  cluster.rows = T,
                  cluster.columns = T,
                  gene.set.use = rownames(Dog.combined.singlet.myeloid@
                                            assays$escape.ssGSEA@
                                            data)[c(1:89)], #c(2:6,10:11,44,52,63,76,59)
                  assay = "escape.ssGSEA")

geyserEnrichment(Dog.combined.singlet.myeloid, 
                 assay = "escape.ssGSEA",scale = T,#f
                 #facet.by = "orig.ident",
                 color.by = "T1-Interferon", #palette = "Spectral",
                 gene.set = "T1-Interferon", order.by = "mean") + NoLegend() 

ridgeEnrichment(Dog.combined.singlet.myeloid, 
                assay = "escape.ssGSEA",
                gene.set = "T1-Interferon",
                add.rug = TRUE,
                scale = TRUE)

densityEnrichment(Dog.combined.singlet.myeloid, 
                  gene.set.use = "T1-Interferon", 
                  group.by = "orig.ident",
                  gene.sets = list)

scatterEnrichment(Dog.combined.singlet.myeloid, 
                  assay = "escape.ssGSEA",
                  x.axis = "T1-Interferon",scale = T,
                  y.axis = "T2-Interferon", facet.by = "Subset5",
                  style = "hex")

Dog.combined.singlet.myeloid <- performPCA(Dog.combined.singlet.myeloid, 
                                           assay = "escape.ssGSEA",
                                           n.dim = 1:30)
pcaEnrichment(Dog.combined.singlet.myeloid, 
              dimRed = "escape.PCA",
              x.axis = "PC1",
              y.axis = "PC2")

pcaEnrichment(Dog.combined.singlet.myeloid, 
              dimRed = "escape.PCA",
              x.axis = "PC1",
              y.axis = "PC2",
              add.percent.contribution = TRUE,
              display.factors = TRUE,
              number.of.factors = 30)

Dog.combined.singlet.myeloid <- performNormalization(Dog.combined.singlet.myeloid, 
                                                     assay = "escape.ssGSEA", 
                                                     gene.sets = list, 
                                                     scale.factor = Dog.combined.singlet.myeloid$nFeature_RNA)

all.markers <- FindAllMarkers(Dog.combined.singlet.myeloid, 
                              assay = "escape.ssGSEA_normalized", 
                              min.pct = 0,
                              logfc.threshold = 0)

View(all.markers)


#
rm(tmp)
tmp <- read.csv("Data_Healthy/GO/PDL1_Mo_DEG.csv")

tmp <- tmp %>% #create new categorical column 
  mutate(Trend = ifelse(p_val < 0.05 & avg_log2FC > 0.5, "Up",
                        ifelse(p_val < 0.05 & avg_log2FC < -0.5, "Down", "None")))
filter <- subset(tmp, p_val < 0.05 & abs(avg_log2FC) > 0)
top10 <- filter %>% top_n(n =10, wt = avg_log2FC) 
bottom10 <- filter %>% top_n(n =10, wt = -(avg_log2FC))
tmp$names <- rownames(tmp)
top10$names <- rownames(top10)
bottom10$names <- rownames(bottom10)
tmp %>% dplyr::count(Trend) #obtain Trend counts
tmp %>% dplyr::distinct(Trend) %>% pull() #check Trend categories

FeaturePlot(Dog.combined.singlet.myeloid, "IGLV3-21")

Dog.combined.singlet.myeloid@active.ident<-Dog.combined.singlet.myeloid$seurat_clusters
library(ggrepel)
# volcano plot https://blog.naver.com/pickyu2/223072025364
cluster_marker_mo <- FindMarkers(Dog.combined.singlet.myeloid,
                                 ident.1 = c(20),
                                 # ident.2 = c(0,1,2,3,6,7,11,19,17),
                                 min.pct = 0.3)
write.table(cluster_marker_mo, "Data_Healthy/C20_DEGs.csv", append = F, sep = ",")

cluster_marker_mo_ifn <- FindMarkers(Dog.combined.singlet.myeloid,
                                     ident.1 = c(18),
                                     ident.2 = c(14,16,5,9,8,12,10,4,13),
                                     min.pct = 0.3)
write.table(cluster_marker_mo_ifn, "Data_Healthy/C18_DEGs.csv", append = F, sep = ",")

cluster_marker_mo_ifn <- FindMarkers(Dog.combined.singlet.myeloid,
                                     ident.1 = c(15,21,22),
                                     # ident.2 = c(14,16,5,9,18,8,12,10,4,13),
                                     min.pct = 0.3)


write.table(cluster_marker_mo_ifn, "Data_Healthy/C152122_DEGs.csv", append = F, sep = ",")
head(cluster_marker_mo_ifn,10)
rm(tmp)
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")

DEG <- read.csv("PATH_TO_FILE")

DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")
DEG <- read.csv("PATH_TO_FILE")

tmp <- read.csv("Data_Healthy/C20_DEGs.csv")
tmp <- read.csv("Data_Healthy/C18_DEGs.csv")
tmp <- read.csv("Data_Healthy/C152122_DEGs.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/CD4/cluster14.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/cluster7.csv")
tmp <- read.csv("Data_Healthy/LAG3specific_CD8only_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/Cluster9_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/Cluster6_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/Cluster45_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/Cluster48_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/Cluster5_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/CD8/Tgd_C17_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/plasma458_DEG.csv")
tmp <- read.csv("Data_Healthy/Featureplot_Cluster_Marker/B/plasma5_DEG.csv")
tmp <- read.csv("Data_Healthy/CD4_CTLA4_marker.csv")

tmp <- tmp[-1,]

library(ggrepel)
tmp_updw <- subset(tmp, Trend %in% c("Up","Down"))
updwgene <- intersect(DEG$Genes, tmp_updw$names)
updwgene <- as.data.frame(updwgene)

tmp$names <- rownames(tmp)
tmp <- tmp %>% mutate(Trend = ifelse(p_val<0.05 &  avg_log2FC >1, "Up",
                                     ifelse(p_val<0.05 &avg_log2FC < -1, "Down","None")))
top10 <- subset(tmp, Trend %in% c("Up")) %>% top_n(n = 10, wt = avg_log2FC)
bottom10 <- subset(tmp, Trend %in% c("Down")) %>% top_n(n=10, wt = -avg_log2FC)
ggplot(tmp, aes(x=avg_log2FC, 
                y=-log10(p_val), 
                col = Trend,
                size = Trend)) + 
  geom_point(size =4) +
  scale_color_manual(values = c("Down" = "#0348A6","None" = "grey","Up"= "#FF4B20")) +
  scale_size_manual(values = c("Up" = 2,"Down" = 2, "None" = 1)) + 
  theme_classic(base_size = 15) +
  geom_vline(xintercept = c(-1, 1), lty = 2, col= "black") +   
  geom_hline(yintercept = -log10(0.05), lty = 2, col= "black") +
  labs(x=(Log[2]~Fold~Change), y=(Log[10]~P~value)) + 
  # scale_y_sqrt() +
  # scale_x_continuous(limits = c(-2, 2)) + 
  NoLegend()  + 
  #geom_text_repel(data = subset(tmp, names %in% c(as.character(top10$names), as.character(bottom10$names))), aes(label =names), size = 3, point.padding = unit(2, "lines"), box.padding = unit(0.5, "lines"),max.overlaps = 30, force = 2) # +
  #geom_label_repel(data = subset(tmp, names %in% c(updwgene$updwgene)), aes(label = names), size = 5, colour = "black", nudge_y = 1, nudge_x =1, max.overlaps = 20)  +
  geom_label_repel(data = subset(tmp, names %in% c(DEG$Genes)), aes(label = names), size = 4, colour = "#FF4B20", nudge_y = 0.5, nudge_x = 0.5, max.overlaps = 1000) # +
# + scale_y_continuous(limits = c(-10, 10))
ggsave("Data_Healthy/Volcano_CD4_Cluster14.png",width = 5, height = 5 )
ggsave("Data_Healthy/Volcano_Myeloid_IMHA_MMVD_DEG_LABEL.png",width = 5, height = 5 )
ggsave("Data_Healthy/Volcano_Myeloid_IMHA_DEG.png",width = 15, height = 5)
ggsave("Data_Healthy/Volcano_Myeloid_CD18_MMVD_DEG.png",width = 5, height = 5)
ggsave("Data_Healthy/Volcano_CD4_CTLA4_DEG.png",width = 5, height = 5)

ggsave("Data_Healthy/Volcano_CD8NKT_Lag3specific_FromPD1.png",width = 5, height = 5 )
ggsave("Data_Healthy/Volcano_B_C9.png",width = 5, height = 5 )

VlnPlot(Dog.combined.singlet.CD8NKT, features = c("PMAIP1"), #pt.size = 0.1, 
        cols =  c("#ADDAE6","#E63222"),
        split.by = "PD1_LAG3", 
        label = T) 

DimPlot(Dog.combined.singlet.CD8NKT, reduction = "tsne")
#
Dog.combined.singlet.myeloid$cd274 <- Dog.combined.singlet.myeloid@assays$RNA$data["CD274",]
Dog.combined.singlet.myeloid$PDL1_Sep <- ifelse(Dog.combined.singlet.myeloid$cd274>0, "Pos","Neg")

FeaturePlot(Dog.combined.singlet.myeloid, features = c("PLD4"), #pt.size = 0.1, 
            order = T, cols =  c("#ADDAE6","#E63222"),max.cutoff = 1,
            # split.by = "PDL1_Sep", 
            label = T)
ggsave("Data_Healthy/Feature_Myeloid_group.by_PDL1.png", width = 6, height = 3)
Dog.combined.singlet.myeloid@active.ident <- Dog.combined.singlet.myeloid$seurat_clusters

DimPlot(Dog.combined.singlet.myeloid, group.by = "PDL1_Sep", pt.size = 2,
        cols =  c("#ADDAE6","#E63222"), order = T) + NoLegend()
ggsave("Data_Healthy/Umap_Myeloid_group.by_PDL1.png", width = 5, height = 5)

Dog.combined.singlet.myeloid.PMN <- subset(Dog.combined.singlet.myeloid, 
                                           seurat_clusters%in% c(0,1,2,6,11,7,3,19,20,17))
Dog.combined.singlet.myeloid.PMN$cd274_2 <-Dog.combined.singlet.myeloid.PMN@assays$RNA$data["CD274",]
Dog.combined.singlet.myeloid.PMN$PDL1_Sep_2 <- ifelse(Dog.combined.singlet.myeloid.PMN$cd274_2>0, "Pos","Neg")
Idents(Dog.combined.singlet.myeloid.PMN) <- Dog.combined.singlet.myeloid.PMN$PDL1_Sep_2
tmp <- FindMarkers(Dog.combined.singlet.myeloid.PMN, min.pct = 0.3, ident.1 = "Pos", ident.2 = "Neg")
write.table(tmp, "Data_Healthy/PDL1_PMN_DEG.csv", append = F, sep = ",")

Dog.combined.singlet.myeloid.Mo <- subset(Dog.combined.singlet.myeloid, 
                                          seurat_clusters%in% c(14,16,5,9,8,18,12))
Dog.combined.singlet.myeloid.Mo$cd274_3 <-Dog.combined.singlet.myeloid.Mo@assays$RNA$data["CD274",]
Dog.combined.singlet.myeloid.Mo$PDL1_Sep_3 <- ifelse(Dog.combined.singlet.myeloid.Mo$cd274_3>0, "Pos","Neg")
Idents(Dog.combined.singlet.myeloid.Mo) <- Dog.combined.singlet.myeloid.Mo$PDL1_Sep_3
tmp <- FindMarkers(Dog.combined.singlet.myeloid.Mo, min.pct = 0.3, ident.1 = "Pos", ident.2 = "Neg")
write.table(tmp, "Data_Healthy/PDL1_Mo_DEG.csv", append = F, sep = ",")

#

DimPlot(Dog.combined.singlet.CD8NKT, label = T)
Dog.combined.singlet.CD8NKT$pd1 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["PDCD1",]
Dog.combined.singlet.CD8NKT$lag3 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["LAG3",]
Dog.combined.singlet.CD8NKT$tigit <- Dog.combined.singlet.CD8NKT@assays$RNA$data["TIGIT",]
Dog.combined.singlet.CD8NKT$ctla4 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["CTLA4",]
Dog.combined.singlet.CD8NKT$cd274 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["CD274",]
Dog.combined.singlet.CD8NKT$havcr2 <- Dog.combined.singlet.CD8NKT@assays$RNA$data["HAVCR2",]
Dog.combined.singlet.CD8NKT$tox <- Dog.combined.singlet.CD8NKT@assays$RNA$data["TOX",]

Dog.combined.singlet.CD8NKT$PD1 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PDL1 <- ifelse(Dog.combined.singlet.CD8NKT$cd274>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_CTLA4 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                  Dog.combined.singlet.CD8NKT$ctla4>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_LAG3 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                 Dog.combined.singlet.CD8NKT$lag3>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_TIM3 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                 Dog.combined.singlet.CD8NKT$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_TIGIT <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                  Dog.combined.singlet.CD8NKT$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_TOX <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                Dog.combined.singlet.CD8NKT$tox>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_LAG3_TIM3 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                      Dog.combined.singlet.CD8NKT$lag3>0&
                                                      Dog.combined.singlet.CD8NKT$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_LAG3_TIGIT <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                       Dog.combined.singlet.CD8NKT$lag3>0&
                                                       Dog.combined.singlet.CD8NKT$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_TIM3_TIGIT <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                       Dog.combined.singlet.CD8NKT$havcr2>0&
                                                       Dog.combined.singlet.CD8NKT$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_CTLA4_TIGIT <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                        Dog.combined.singlet.CD8NKT$ctla4>0&
                                                        Dog.combined.singlet.CD8NKT$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_CTLA4_TIM3 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                       Dog.combined.singlet.CD8NKT$ctla4>0&
                                                       Dog.combined.singlet.CD8NKT$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$PD1_CTLA4_LAG3 <- ifelse(Dog.combined.singlet.CD8NKT$pd1>0& 
                                                       Dog.combined.singlet.CD8NKT$ctla4>0&
                                                       Dog.combined.singlet.CD8NKT$lag3>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$Terminal_Exhaustion <- ifelse(Dog.combined.singlet.CD8NKT$T_Cell_Terminal_Differentiation>0& 
                                                            Dog.combined.singlet.CD8NKT$Exhaustion>0,"Pos","Neg")

Dog.combined.singlet.CD8NKT$Terminal <- ifelse(Dog.combined.singlet.CD8NKT$T_Cell_Terminal_Differentiation>0,"Pos","Neg")
Dog.combined.singlet.CD8NKT$Exh <- ifelse(Dog.combined.singlet.CD8NKT$Exhaustion>0,"Pos","Neg")


terexh.marker.from.exh <- subset(Dog.combined.singlet.CD8NKT, Terminal %in% c("Pos"))
DimPlot(terexh.marker.from.exh, group.by = "Exh")

FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("T_Cell_Terminal_Differentiation"), cols = viridis(100))

Idents(clusterExh) <- clusterExh$Terminal_Exhaustion
terexh.marker.from.exh <- FindMarkers(clusterExh, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)

Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$seurat_clusters
DimPlot(Dog.combined.singlet.CD8NKT,)

Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$Terminal_Exhaustion
terexh.marker <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(terexh.marker, "Data_Healthy/CD8NKT_Ter_Exh.csv", sep = ",", append = F)
VlnPlot(Dog.combined.singlet.CD8NKT, features = c("CD160"), cols =  c("#ADDAE6","#E63222"), pt.size = 0)

Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$PD1
terexh.marker <- FindMarkers(subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(3,5,0,1,7,8,14,4,2,6)), ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(terexh.marker, "Data_Healthy/CD8NKT_PDCD1_DEG.csv", sep = ",", append = F)

Idents(Dog.combined.singlet.CD8NKT) <- Dog.combined.singlet.CD8NKT$PD1
terexh.marker <- FindMarkers(subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(9,12)), ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(terexh.marker, "Data_Healthy/T_cycle_PDCD1_DEG.csv", sep = ",", append = F)

FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("PDCD1"), split.by = "PD1", 
            cols =  c("#ADDAE6","#E63222"))

#
Dog.combined.singlet.CD4$pd1 <- Dog.combined.singlet.CD4@assays$RNA$data["PDCD1",]
Dog.combined.singlet.CD4$lag3 <- Dog.combined.singlet.CD4@assays$RNA$data["LAG3",]
Dog.combined.singlet.CD4$tigit <- Dog.combined.singlet.CD4@assays$RNA$data["TIGIT",]
Dog.combined.singlet.CD4$ctla4 <- Dog.combined.singlet.CD4@assays$RNA$data["CTLA4",]
Dog.combined.singlet.CD4$cd274 <- Dog.combined.singlet.CD4@assays$RNA$data["CD274",]
Dog.combined.singlet.CD4$havcr2 <- Dog.combined.singlet.CD4@assays$RNA$data["HAVCR2",]
Dog.combined.singlet.CD4$tox <- Dog.combined.singlet.CD4@assays$RNA$data["TOX",]

Dog.combined.singlet.CD4$PD1 <- ifelse(Dog.combined.singlet.CD4$pd1>0,"Pos","Neg")
Dog.combined.singlet.CD4$PDL1 <- ifelse(Dog.combined.singlet.CD4$cd274>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_CTLA4 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                               Dog.combined.singlet.CD4$ctla4>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_LAG3 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                              Dog.combined.singlet.CD4$lag3>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_TIM3 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                              Dog.combined.singlet.CD4$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_TIGIT <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                               Dog.combined.singlet.CD4$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_TOX <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                             Dog.combined.singlet.CD4$tox>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_LAG3_TIM3 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                   Dog.combined.singlet.CD4$lag3>0&
                                                   Dog.combined.singlet.CD4$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_LAG3_TIGIT <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                    Dog.combined.singlet.CD4$lag3>0&
                                                    Dog.combined.singlet.CD4$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_TIM3_TIGIT <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                    Dog.combined.singlet.CD4$havcr2>0&
                                                    Dog.combined.singlet.CD4$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_CTLA4_TIGIT <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                     Dog.combined.singlet.CD4$ctla4>0&
                                                     Dog.combined.singlet.CD4$tigit>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_CTLA4_TIM3 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                    Dog.combined.singlet.CD4$ctla4>0&
                                                    Dog.combined.singlet.CD4$havcr2>0,"Pos","Neg")
Dog.combined.singlet.CD4$PD1_CTLA4_LAG3 <- ifelse(Dog.combined.singlet.CD4$pd1>0& 
                                                    Dog.combined.singlet.CD4$ctla4>0&
                                                    Dog.combined.singlet.CD4$lag3>0,"Pos","Neg")

#

Dog.combined.singlet.CD4 <- AddMetaData(Dog.combined.singlet.CD4, ES.immune) # add metadata
meta <- Dog.combined.singlet.CD4[[]]
meta$seurat_clusters <- Dog.combined.singlet.CD4@active.ident
View(meta)
meta <- meta[,c(16,162:166,168:172,60,62,63,73:75,183)] 
ncol(meta)
heatmap <- meta %>% 
  group_by(
    PD1_TIM3_TIGIT, PD1_LAG3_TIGIT, #PD1_LAG3_TIM3, 
    PD1_CTLA4_TIGIT, PD1_CTLA4_TIM3, PD1_CTLA4_LAG3,
    PD1_LAG3, PD1_TIM3, PD1_TIGIT, PD1_CTLA4,# PD1_TOX,
    PD1, PDL1,seurat_clusters
  ) %>% summarise(across(1:6, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 1:77) #rownames(heatmap) <- paste0("X", 1:22)
View(heatmap)

headers <- heatmap[,c(1:12)] 
heatmap <- heatmap[,c(13:18)]

PD1 <- unique(Dog.combined.singlet.CD4$PD1) 
PDL1 <- unique(Dog.combined.singlet.CD4$PDL1) 

PD1_TOX <- unique(Dog.combined.singlet.CD4$PD1_TOX) 
PD1_LAG3 <- unique(Dog.combined.singlet.CD4$PD1_LAG3) 
PD1_TIM3 <- unique(Dog.combined.singlet.CD4$PD1_TIM3) 
PD1_TIGIT <- unique(Dog.combined.singlet.CD4$PD1_TIGIT) 
PD1_CTLA4 <- unique(Dog.combined.singlet.CD4$PD1_CTLA4) 

PD1_LAG3_TIM3 <- unique(Dog.combined.singlet.CD4$PD1_LAG3_TIM3) 
PD1_LAG3_TIGIT <- unique(Dog.combined.singlet.CD4$PD1_LAG3_TIGIT) 
PD1_TIM3_TIGIT <- unique(Dog.combined.singlet.CD4$PD1_TIM3_TIGIT) 
PD1_CTLA4_TIGIT <- unique(Dog.combined.singlet.CD4$PD1_CTLA4_TIGIT) 
PD1_CTLA4_TIM3 <- unique(Dog.combined.singlet.CD4$PD1_CTLA4_TIM3) 
PD1_CTLA4_LAG3 <- unique(Dog.combined.singlet.CD4$PD1_CTLA4_LAG3) 

PD1colors <- c("#FCE540", "#3F114E") 
names(PD1colors) <- c("Pos","Neg")
PDL1colors <- c("#FCE540", "#3F114E") 
names(PDL1colors) <- c("Pos","Neg")
PD1_LAG3colors <-c("#FCE540", "#3F114E") 
names(PD1_LAG3colors) <-  c("Pos","Neg") 
PD1_TIM3colors <- c("#FCE540", "#3F114E")  
names(PD1_TIM3colors) <-  c("Pos","Neg") 
PD1_TIGITcolors <-c("#FCE540", "#3F114E") 
names(PD1_TIGITcolors) <-  c("Pos","Neg") 
PD1_CTLA4colors <-c("#FCE540", "#3F114E") 
names(PD1_CTLA4colors) <-  c("Pos","Neg") 
PD1_LAG3_TIM3colors <-c("#FCE540", "#3F114E") 
names(PD1_LAG3_TIM3colors) <- c("Pos","Neg")
PD1_LAG3_TIGITcolors <- c("#FCE540", "#3F114E") 
names(PD1_LAG3_TIGITcolors) <- c("Pos","Neg")
PD1_TIM3_TIGITcolors <- c("#FCE540", "#3F114E") 
names(PD1_TIM3_TIGITcolors) <- c("Pos","Neg")
PD1_CTLA4_TIGITcolors <-c("#FCE540", "#3F114E") 
names(PD1_CTLA4_TIGITcolors) <- c("Pos","Neg")
PD1_CTLA4_TIM3colors <- c("#FCE540", "#3F114E") 
names(PD1_CTLA4_TIM3colors) <- c("Pos","Neg")
PD1_CTLA4_LAG3colors <- c("#FCE540", "#3F114E") 
names(PD1_CTLA4_LAG3colors) <- c("Pos","Neg")

clusterColors <- c(scales::hue_pal()(18))
names(clusterColors) <- 0:17

colors <- list(PD1 = PD1colors,
               PDL1=PDL1colors,
               PD1_LAG3=PD1_LAG3colors,
               PD1_TIM3=PD1_TIM3colors,
               PD1_TIGIT=PD1_TIGITcolors,
               PD1_CTLA4=PD1_CTLA4colors,
               PD1_LAG3_TIM3=PD1_LAG3_TIM3colors,
               PD1_LAG3_TIGIT=PD1_LAG3_TIGITcolors,
               PD1_TIM3_TIGIT=PD1_TIM3_TIGITcolors,
               PD1_CTLA4_LAG3=PD1_CTLA4_LAG3colors,
               PD1_CTLA4_TIM3=PD1_CTLA4_TIM3colors,
               PD1_CTLA4_TIGIT=PD1_CTLA4_TIGITcolors,
               seurat_clusters = clusterColors
)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 1:77) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/GSEA_CD4_Ehx.pdf", width = 10, height =8) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "none", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   #cellwidth = 10, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   ))
dev.off()


Dog.combined.singlet.CD4@active.ident <- Dog.combined.singlet.CD4$seurat_clusters
subset.for.CD4 <- subset(Dog.combined.singlet.CD4, seurat_clusters %in% c(0,1,11,15,
                                                                          4,2,3,
                                                                          14,12,
                                                                          #7,17,
                                                                          5,8,6,13,9))
Idents(subset.for.CD4) <- subset.for.CD4$PD1
marker <- FindMarkers(subset.for.CD4, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/CD4/DEG_PD1_Pos_Neg_selectedclusters_withoutTreg.csv", append = F, sep = ",")
FeaturePlot(Dog.combined.singlet.CD4, features = c("LAG3"))

Idents(subset.for.CD4) <- subset.for.CD4$seurat_clusters
VlnPlot(subset.for.CD4, split.by = "", features = "LAG3")

PD1pos <- subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos"))
DimPlot(PD1pos, split.by = "PD1_LAG3")
FeaturePlot(PD1pos, split.by = "PD1_LAG3", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","LAG3"), )
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_LAG3.png", width = 5, height = 5)

Idents(PD1pos) <- PD1pos$PD1_LAG3
marker <- FindMarkers(PD1pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/LAG3specific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_LAG3)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD4/Cellnumber_PD1_LAG3_CD4.csv", sep = ",", append=F)



TIM3pos <- subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos"))
Idents(TIM3pos) <- TIM3pos$PD1_TIM3
FeaturePlot(TIM3pos, split.by = "PD1_TIM3", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","HAVCR2"))
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_HAVCR2.png", width = 5, height = 5)
marker <- FindMarkers(TIM3pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/TIM3specific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(TIM3pos), TIM3pos$PD1_TIM3)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD4/Cellnumber_PD1_TIM3_CD4.csv", sep = ",", append=F)


CTLA4pos <- subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos"))
Idents(CTLA4pos) <- CTLA4pos$PD1_CTLA4
FeaturePlot(CTLA4pos, split.by = "PD1_CTLA4", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","CTLA4"))
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_CTLA4.png", width = 5, height = 5)

marker <- FindMarkers(CTLA4pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/CTLA4specific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(CTLA4pos), CTLA4pos$PD1_CTLA4)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD4/Cellnumber_PD1_CTLA4_CD4.csv", sep = ",", append=F)



TIGITpos <- subset(Dog.combined.singlet.CD4, PD1 %in% c("Pos"))
Idents(TIGITpos) <- TIGITpos$PD1_TIGIT
FeaturePlot(TIGITpos, split.by = "PD1_TIGIT", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","TIGIT"))
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_TIGIT.png", width = 5, height = 5)

marker <- FindMarkers(TIGITpos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/TIGITspecific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(TIGITpos), TIGITpos$PD1_TIGIT)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD4/Cellnumber_PD1_TIGIT_CD4.csv", sep = ",", append=F)

TOXpos <- subset(Dog.combined.singlet.CD4, TOX %in% c("Pos"))
Idents(TOXpos) <- TOXpos$PD1_TOX
FeaturePlot(TOXpos, split.by = "PD1_TOX", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","rna_TOX"))
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_TOX.png", width = 5, height = 5)
marker <- FindMarkers(TOXpos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/CD4/TOXspecific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(TOXpos), TOXpos$PD1_TOX)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD4/Cellnumber_PD1_TOX_CD4.csv", sep = ",", append=F)


CTLA4marker <- FindMarkers(Dog.combined.singlet.CD4, ident.1 = c(13), ident.2 = c(3,5), min.pct = 0.3)
write.table(CTLA4marker,"Data_Healthy/DEG_CD4_13vs_3_5.csv", append = F, sep = ",")

#
CTLA4pos <- subset(Dog.combined.singlet.CD4, CTLA4 %in% c("Pos"))
Idents(TOXpos) <- TOXpos$PD1_TOX
FeaturePlot(TOXpos, split.by = "PD1_TOX", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","rna_TOX"))
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_TOX.png", width = 5, height = 5)
marker <- FindMarkers(TOXpos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/CD4/TOXspecific_CD4_DEG.csv", append = F, sep = ",")
#
VlnPlot(PD1pos, features = c("MDH1"), group.by = "PD1_LAG3",
        cols = c("#3F114E","#FCE540")) + 
  stat_compare_means(method = "t.test", p.adjust.method = "BH") + 
  NoLegend()
ggsave("Data_Healthy/Lag3specific_MDH1_PD1CD4.png", width = 2, height = 2)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

VlnPlot(PD1pos, features = "MDH1", group.by = "PD1_LAG3", cols = c("#3F114E","#FCE540")) +
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    p.adjust.method = "BH",
    # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  )

VlnPlot(PD1pos, features = "MDH1", group.by = "PD1_LAG3", cols = c("#3F114E","#FCE540")) +
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", step.increase = 0.05) +
  ylim(0, 4000)
Idents(PD1pos) <- PD1pos$seurat_clusters
Idents(TIM3pos) <- PD1pos$seurat_clusters

VlnPlot(TIM3pos, features = c("MDH1"), split.by = "PD1_LAG3", idents = c(7,12,5), 
        sort = T, split.plot = T, cols = c("#3F114E","#FCE540")) +
  #stat_compare_means(method = "t.test") + 
  NoLegend()
ggsave("Data_Healthy/Lag3specific_MDH1_PD1CD4.png", width = 3, height = 2.5)

Idents(TIGITpos) <- TIGITpos$seurat_clusters
Idents(TIGITpos) <- TIGITpos$PD1_TIGIT
TIGITpos@active.ident <-  factor(TIGITpos@active.ident, 
                                 levels = c("Neg", "Pos"))

CTLA4pos@active.ident <-  factor(CTLA4pos@active.ident, 
                                 levels = c("Neg", "Pos"))

VlnPlot(TIGITpos, features = c("OASL"), split.by = "PD1_TIGIT", #idents = c(7,12,5), 
        sort = F, cols = c("#3F114E","#FCE540")) +   NoLegend() #+ stat_compare_means(method = "t.test") + 
NoLegend()
ggsave("Data_Healthy/CD4/TIGITspecific_OASL_PD1CD4.png", width =2, height = 2)

VlnPlot(CTLA4pos, features = c("CTLA4"), group.by = "PD1_CTLA4", #idents = c(7,12,5), 
        sort = F, cols = c("#3F114E","#FCE540")) +   NoLegend() + stat_compare_means(method = "t.test") + 
  NoLegend()

ggsave("Data_Healthy/CD4/CTLA4specific_FUT7_PD1CD4.png", width =2, height = 2)


#
Dog.combined.singlet.CD8NKT@active.ident <- Dog.combined.singlet.CD8NKT$seurat_clusters
meta <- Dog.combined.singlet.CD8NKT[[]]
meta$seurat_clusters <- Dog.combined.singlet.CD8NKT@active.ident
View(meta)
meta <- meta[,c(16,22:24,37:39,137:149)]  #25:27
meta <- meta[,c(16,61:63,65,66,75:77,178,180:184)]  #25:27
ncol(meta)
heatmap <- meta %>% 
  group_by(
    #PD1_TIM3_TIGIT, 
    #  PD1_LAG3_TIGIT,
    #PD1_LAG3_TIM3, 
    # PD1_CTLA4_TIGIT, #PD1_CTLA4_TIM3, PD1_CTLA4_LAG3,
    PD1_LAG3, PD1_TIM3, PD1_TIGIT, PD1_CTLA4, PD1_TOX,
    PD1, seurat_clusters
  ) %>% summarise(across(1:8, mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap)
rownames(heatmap) <- paste0("X", 1:82) 
View(heatmap)

headers <- heatmap[,c(1:7)] 
heatmap <- heatmap[,c(8:15)]
heatmap <- heatmap[c(1:9,15,
                     20:28,32,
                     34:41,47,
                     50:55,
                     62,
                     68,
                     71:73,
                     75,
                     80:83,
                     87:88,
                     93,95),
                   c(11:16)]

PD1 <- unique(Dog.combined.singlet.CD8NKT$PD1) 
PDL1 <- unique(Dog.combined.singlet.CD8NKT$PDL1) 
PD1_TOX <- unique(Dog.combined.singlet.CD8NKT$PD1_TOX) 
PD1_LAG3 <- unique(Dog.combined.singlet.CD8NKT$PD1_LAG3) 
PD1_TIM3 <- unique(Dog.combined.singlet.CD8NKT$PD1_TIM3) 
PD1_TIGIT <- unique(Dog.combined.singlet.CD8NKT$PD1_TIGIT) 
PD1_CTLA4 <- unique(Dog.combined.singlet.CD8NKT$PD1_CTLA4) 

PD1_LAG3_TIM3 <- unique(Dog.combined.singlet.CD8NKT$PD1_LAG3_TIM3) 
PD1_LAG3_TIGIT <- unique(Dog.combined.singlet.CD8NKT$PD1_LAG3_TIGIT) 
PD1_TIM3_TIGIT <- unique(Dog.combined.singlet.CD8NKT$PD1_TIM3_TIGIT) 
PD1_CTLA4_TIGIT <- unique(Dog.combined.singlet.CD8NKT$PD1_CTLA4_TIGIT) 
PD1_CTLA4_TIM3 <- unique(Dog.combined.singlet.CD8NKT$PD1_CTLA4_TIM3) 
PD1_CTLA4_LAG3 <- unique(Dog.combined.singlet.CD8NKT$PD1_CTLA4_LAG3) 

PD1colors <- c("#FCE540", "#3F114E") 
names(PD1colors) <- c("Pos","Neg")
PDL1colors <- c("#FCE540", "#3F114E") 
names(PDL1colors) <- c("Pos","Neg")
PD1_TOXcolors <-c("#FCE540", "#3F114E") 
names(PD1_TOXcolors) <-  c("Pos","Neg") 
PD1_LAG3colors <-c("#FCE540", "#3F114E") 
names(PD1_LAG3colors) <-  c("Pos","Neg") 
PD1_TIM3colors <- c("#FCE540", "#3F114E")  
names(PD1_TIM3colors) <-  c("Pos","Neg") 
PD1_TIGITcolors <-c("#FCE540", "#3F114E") 
names(PD1_TIGITcolors) <-  c("Pos","Neg") 
PD1_CTLA4colors <-c("#FCE540", "#3F114E") 
names(PD1_CTLA4colors) <-  c("Pos","Neg") 
PD1_LAG3_TIM3colors <-c("#FCE540", "#3F114E") 
names(PD1_LAG3_TIM3colors) <- c("Pos","Neg")
PD1_LAG3_TIGITcolors <- c("#FCE540", "#3F114E") 
names(PD1_LAG3_TIGITcolors) <- c("Pos","Neg")
PD1_TIM3_TIGITcolors <- c("#FCE540", "#3F114E") 
names(PD1_TIM3_TIGITcolors) <- c("Pos","Neg")
PD1_CTLA4_TIGITcolors <-c("#FCE540", "#3F114E") 
names(PD1_CTLA4_TIGITcolors) <- c("Pos","Neg")
PD1_CTLA4_TIM3colors <- c("#FCE540", "#3F114E") 
names(PD1_CTLA4_TIM3colors) <- c("Pos","Neg")
PD1_CTLA4_LAG3colors <- c("#FCE540", "#3F114E") 
names(PD1_CTLA4_LAG3colors) <- c("Pos","Neg")

clusterColors <- c(scales::hue_pal()(19))
names(clusterColors) <- 0:18

colors <- list(PD1 = PD1colors,
               PDL1=PDL1colors,
               PD1_TOX=PD1_TOXcolors,
               PD1_LAG3=PD1_LAG3colors,
               PD1_TIM3=PD1_TIM3colors,
               PD1_TIGIT=PD1_TIGITcolors,
               PD1_CTLA4=PD1_CTLA4colors,
               PD1_LAG3_TIM3=PD1_LAG3_TIM3colors,
               PD1_LAG3_TIGIT=PD1_LAG3_TIGITcolors,
               PD1_TIM3_TIGIT=PD1_TIM3_TIGITcolors,
               PD1_CTLA4_LAG3=PD1_CTLA4_LAG3colors,
               PD1_CTLA4_TIM3=PD1_CTLA4_TIM3colors,
               PD1_CTLA4_TIGIT=PD1_CTLA4_TIGITcolors,
               seurat_clusters = clusterColors
)

normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}
heatmap2 <- sapply(heatmap, normalize)
heatmap2  <- heatmap2 [,colSums(is.na(heatmap2 ))<nrow(heatmap2)]
rownames(heatmap2) <- paste0("X", 1:97) #rownames(heatmap2) <- paste0("X", 1:22)
rownames(heatmap2) <- paste0("X", 1:82) #rownames(heatmap2) <- paste0("X", 1:22)

rownames(heatmap2) <- paste0("X",c(1:9,15,
                                   20:28,32,
                                   34:41,47,
                                   50:55,
                                   62,
                                   68,
                                   71:73,
                                   75,
                                   80:83,
                                   87:88,
                                   93,95)) #rownames(heatmap2) <- paste0("X", 1:22)

pdf("Data_Healthy/Featureplot_Cluster_Marker/CD8/SEA_CD8_Ehx.pdf", width = 20, height =8) #5.5 x 3.8, 9x4
pheatmap::pheatmap(t(heatmap2), scale = "none", show_colnames = T,
                   annotation_col = headers,
                   cluster_cols= T,
                   annotation_colors = colors,
                   fontsize = 12, cluster_rows =T,
                   legend = T, legend_labels = T, annotation_legend = T,
                   cellwidth = 10, 
                   cellheight = 12,
                   # color = viridis(100)
                   color = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(100)
                   ))
dev.off()

TIM3pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos") & seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) )
Idents(TIM3pos) <- TIM3pos$PD1_TIM3
FeaturePlot(TIM3pos, split.by = "PD1_TIM3", cols =  c("#ADDAE6","#E63222"), order = T, 
            features = c("PDCD1","HAVCR2"))
ggsave("Data_Healthy/Featureplot_CD8_PDCD1_HAVCR2.png", width = 5, height = 5)
Cellnumbers.subset <- table(Idents(TIM3pos), TIM3pos$PD1_TIM3)

LAG3pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos") & seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) )
Idents(LAG3pos) <- LAG3pos$PD1_LAG3
FeaturePlot(LAG3pos, split.by = "PD1_LAG3", cols =  c("#ADDAE6","#E63222"), order = T, 
            features = c("PDCD1","LAG3"))
ggsave("Data_Healthy/Featureplot_CD8_PDCD1_LAG3.png", width = 5, height = 5)
Cellnumbers.subset <- table(Idents(LAG3pos), LAG3pos$PD1_LAG3)

TIGITpos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos") & seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) )
Idents(TIGITpos) <- TIGITpos$PD1_TIGIT
FeaturePlot(TIGITpos, split.by = "PD1_TIGIT", cols =  c("#ADDAE6","#E63222"), order = T, 
            features = c("PDCD1","TIGIT"))
ggsave("Data_Healthy/Featureplot_CD8_PDCD1_TIGIT.png", width = 5, height = 5)
Cellnumbers.subset <- table(Idents(TIGITpos), TIGITpos$PD1_TIGIT)


CTLA4pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos") & seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) )
Idents(CTLA4pos) <- CTLA4pos$PD1_CTLA4
FeaturePlot(CTLA4pos, split.by = "PD1_CTLA4", cols =  c("#ADDAE6","#E63222"), order = T, 
            features = c("PDCD1","CTLA4"))
ggsave("Data_Healthy/Featureplot_CD8_PDCD1_CTLA4.png", width = 5, height = 5)
Cellnumbers.subset <- table(Idents(CTLA4pos), CTLA4pos$PD1_CTLA4)


cluster912 <- FindMarkers(Dog.combined.singlet.CD8NKT, ident.1 = c(12,9), min.pct = 0.3)
write.table(cluster912,"Data_Healthy/Cluster9_12_DEG.csv", append = F, sep = ",")

PD1pos <- subset(Dog.combined.singlet.CD8NKT,seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) & PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_LAG3
DimPlot(PD1pos, split.by = "PD1_LAG3", 
        cols =  c("#ADDAE6","#E63222"),
        order = T) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/CD8_LAG3sep.png", width = 3, height = 3)
marker <- FindMarkers(PD1pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/LAG3specific_CD8only_DEG.csv", append = F, sep = ",")

PD1pos <- subset(Dog.combined.singlet.CD8NKT,seurat_clusters %in% c(9,12) & PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_LAG3
DimPlot(PD1pos, split.by = "PD1_LAG3", 
        cols =  c("#ADDAE6","#E63222"),
        order = T) + NoLegend()
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/C9-12_CD8_LAG3sep.png", width = 3, height = 3)
marker <- FindMarkers(PD1pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/LAG3specific_C9-12_only_DEG.csv", append = F, sep = ",")

FeaturePlot(PD1pos, features = c("LAG3"), split.by = "PD1_LAG3", cols =  c("#ADDAE6","#E63222"))
ggsave("Data_Healthy/Feature_CD8only_PD1_LAG3_Sep.png", width = 6, height = 3)

Idents(PD1pos) <- PD1pos$seurat_clusters
FeaturePlot(PD1pos, features = c("PDCD1","LAG3"), split.by = "PD1_LAG3", label = T, #pt.size = 2,
            cols =  c("#ADDAE6","#E63222"), order = T)
ggsave("Data_Healthy/Feature_CD8only_PD1_LAG3_Sep.png", width = 5, height = 5)


VlnPlot(PD1pos, features = c("LAG3"), group.by = "PD1_LAG3",cols =  c("#ADDAE6","#E63222")) +NoLegend()
ggsave("Data_Healthy/Vln_CD8only_PD1_LAG3_Sep.png", width = 1.5, height = 3)


VlnPlot(PD1pos, features = c("PMAIP1","ANKRD17","TBX21","IRF5","OASL"), pt.size = 0, stack = T,flip = T,
        cols =  c("#E63222","#E63222","#ADDAE6","#E63222","#E63222")
)  +NoLegend()
ggsave("Data_Healthy/Vln_LAG3specific_CD8only_ResToVirus_5genes.png", width = 2.5, height = 4)

TIM3pos <- subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) &PD1 %in% c("Pos"))
Idents(TIM3pos) <- TIM3pos$PD1_TIM3
DimPlot(TIM3pos, split.by = "PD1_TIM3",cols =  c("#ADDAE6","#E63222"),
        order = T)
FeaturePlot(TIM3pos, features = c("PDCD1","HAVCR2"), split.by = "PD1_TIM3", label = T, #pt.size = 2,
            cols =  c("#ADDAE6","#E63222"), order = T)
marker <- FindMarkers(TIM3pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/TIM3specific_CD8only_DEG.csv", append = F, sep = ",")

CTLA4pos <- subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) &PD1 %in% c("Pos"))
Idents(CTLA4pos) <- CTLA4pos$PD1_CTLA4
DimPlot(CTLA4pos, split.by = "PD1_CTLA4",cols =  c("#ADDAE6","#E63222"),
        order = T)
FeaturePlot(CTLA4pos, features = c("PDCD1","CTLA4"), split.by = "PD1_CTLA4", label = T, #pt.size = 2,
            cols =  c("#ADDAE6","#E63222"), order = T)
marker <- FindMarkers(CTLA4pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/CTLA4specific_CD8only_DEG.csv", append = F, sep = ",")

TIGITpos <- subset(Dog.combined.singlet.CD8NKT,seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) & PD1 %in% c("Pos"))
Idents(TIGITpos) <- TIGITpos$PD1_TIGIT
DimPlot(TIGITpos, split.by = "PD1_TIGIT",cols =  c("#ADDAE6","#E63222"),
        order = T)
marker <- FindMarkers(TIGITpos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/TIGITspecific_CD8only_DEG.csv", append = F, sep = ",")

TOXpos <- subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(3,5,0,1,2,4,6,7,8,14) & PD1 %in% c("Pos"))
Idents(TOXpos) <- TIM3pos$PD1_TOX
DimPlot(TOXpos, split.by = "PD1_TOX",cols =  c("#ADDAE6","#E63222"),
        order = T)
marker <- FindMarkers(TOXpos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/TOXspecific_CD8only_DEG.csv", append = F, sep = ",")

FeaturePlot(Dog.combined.singlet.CD8NKT, features = "PDCD1", split.by = "PD1",  order = T,
            cols =  c("#ADDAE6","#E63222") )
ggsave("PDCD1_Pos_Neg.png", width = 8, height = 5)
DimPlot(Dog.combined.singlet.CD8NKT, group.by = "PD1",  pt.size = 0.8,
        cols =  c("#ADDAE6","#E63222"), order = T) +NoLegend()
ggsave("Data_Healthy/CD8+PDCD1_Pos_Neg.png", width = 5, height = 5)

DimPlot(Dog.combined.singlet.CD8NKT, split.by = "PD1",  pt.size = 0.8,
        cols =  c("#ADDAE6","#E63222"), order = T) +NoLegend()
ggsave("Data_Healthy/CD8+PDCD1_Pos_Neg.png", width = 4, height = 4)
#
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet, '0' = "C_0", '1' = "C_1", '2' = "C_2", '3' = "C_3",
                                     '4' = "C_4", '5' = "C_5", '6' = "C_6", '7' = "C_7",
                                     '8' = "C_8", '9'= "C_9", '10' = "C_10", '11' = "C_11",'12' = "C_12",
                                     '13' = "C_13", '14' ="C_14", '15' ="C_15", '16' ="C_16",'17' ="C_17",
                                     '18' ="C_18", '19' ="C_19", '20' ="C_20",  '21' ="C_21", '22' ="C_22",
                                     '23' ="C_23", '24' ="C_24", '25' ="C_25",
                                     '26' = "C_26", '27' = "C_27", '28' = "C_28", '29' = "C_29",
                                     '30' = "C_30", '31' = "C_31", '32' = "C_32", '33' = "C_33",
                                     '34' = "C_34", '35'= "C_35", '36' = "C_36", '37' = "C_37",'38' = "C_38",
                                     '39' = "C_39", '40' = "C_40", '41' = "C_41", '42' = "C_42",
                                     '43' = "C_43", '44' = "C_44", '45' = "C_45", '46' = "C_46",
                                     '47' = "C_47", '48'= "C_48", '49' = "C_49",'50' = "C_50")
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters
DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) #+NoLegend()

Dog.combined.singlet <- RenameIdents(Dog.combined.singlet, '0' = "C_0", '1' = "C_1", '2' = "C_2", '3' = "C_3",
                                     '4' = "C_4", '5' = "C_5", '6' = "C_6", '7' = "C_7",
                                     '8' = "C_8", '9'= "C_9", '10' = "C_10", '11' = "C_11",'12' = "C_12",
                                     '13' = "C_13", '14' ="C_14", '15' ="C_15", '16' ="C_16",'17' ="C_17",
                                     '18' ="Mo", '19' ="Mo", '20' ="C_20",  '21' ="Treg", '22' ="C_22",
                                     '23' ="Mo", '24' ="C_24", '25' ="C_25",
                                     '26' = "C_26", '27' = "Tcycle", '28' = "C_28", '29' = "C_29",
                                     '30' = "C_30", '31' = "C_31", '32' = "C_32", '33' = "C_33",
                                     '34' = "C_34", '35'= "Mo", '36' = "C_36", '37' = "C_37",'38' = "DC",
                                     '39' = "C_39", '40' = "C_40", '41' = "C_41", '42' = "C_42",
                                     '43' = "C_43", '44' = "C_44", '45' = "C_45", '46' = "C_46",
                                     '47' = "C_47", '48'= "C_48", '49' = "pDC",'50' = "C_50")
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters

DimPlot(Dog.combined.singlet, reduction = "tsne", label = T) +NoLegend()
Dog.combined.singlet <- RenameIdents(Dog.combined.singlet, '0' = "CD4", '1' = "CD8", '2' = "CD4", '3' = "CD4_nv",
                                     '4' = "CD4", '5' = "CD4", '6' = "CD8_nv", '7' = "PMN",
                                     '8' = "CD4_nv", '9'= "PMN", '10' = "Mo2", '11' = "PMN",'12' = "CD4_nv",
                                     '13' = "B", '14' ="PMN", '15' ="CD4", '16' ="CD8",'17' ="PMN",
                                     '18' ="Mo", '19' ="Mo", '20' ="PMN",  '21' ="Treg", '22' ="CD8",
                                     '23' ="Mo", '24' ="CD8", '25' ="Un",
                                     '26' = "B", '27' = "Tcycle", '28' = "CD4", '29' = "CD4",
                                     '30' = "PMN", '31' = "CD4_nv", '32' = "CD8", '33' = "Mo2",
                                     '34' = "PMN", '35'= "Mo", '36' = "PMN", '37' = "P",'38' = "DC",
                                     '39' = "Tgd", '40' = "CD4", '41' = "P", '42' = "Mo2",
                                     '43' = "TDN", '44' = "DC2", '45' = "PMN", '46' = "PMN",
                                     '47' = "PMN", '48'= "Tgd", '49' = "pDC",'50' = "Un")
Dog.combined.singlet@active.ident <- Dog.combined.singlet$seurat_clusters

#
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
cellchat.Dog <- readRDS("Data_Healthy/cellchat/cellchat.Dog.rds")
saveRDS(cellchat.Dog, "cellchat.Dog.rds")

data.input <- GetAssayData(Dog.combined.singlet, assay = "RNA", slot = "data")  # Dog.combined.singlet.cellchat
labels <- Idents(Dog.combined.singlet) # Dog.combined.singlet.cellchat
meta <- data.frame(group = labels, row.names = names(labels)) 
colnames(meta) <-  "labels" 
cellchat.Dog <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat.Dog <- addMeta(cellchat.Dog, meta = meta, meta.name = "labels")
cellchat.Dog <- setIdent(cellchat.Dog, ident.use = "labels") 
levels(cellchat.Dog@idents) 
groupSize <- as.numeric(table(cellchat.Dog@idents))


CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat.Dog@DB <- CellChatDB.use
CellChatDB.use
cellchat.Dog <- subsetData(cellchat.Dog) 
future::plan("multisession", workers = 4 ) 
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.Dog <- identifyOverExpressedGenes(cellchat.Dog)
cellchat.Dog <- identifyOverExpressedInteractions(cellchat.Dog)
cellchat.Dog <- projectData(cellchat.Dog, PPI.human)
cellchat.Dog <- computeCommunProb(cellchat.Dog)
cellchat.Dog <- filterCommunication(cellchat.Dog, min.cells = 10) 
cellchat.Dog <- computeCommunProbPathway(cellchat.Dog)
cellchat.Dog <- aggregateNet(cellchat.Dog)
groupSize <- as.numeric(table(cellchat.Dog@idents))
netVisual_circle(cellchat.Dog@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.Dog@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

df.net1 <- subsetCommunication(cellchat.Dog) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.table(df.net1, file="Data_Healthy/cellchat/All_Interaction.csv", sep = ",", append=F)
significantpathways <- as.data.frame(cellchat.Dog@netP$pathways)
write.table(significantpathways, file="Data_Healthy/cellchat/Significant_Pathways.csv", sep = ",", append=F)

cellchat.Dog <- netAnalysis_computeCentrality(cellchat.Dog, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat.Dog, "Data_Healthy/cellchat/cellchat.Dog.label.T.rds")

#
pdf("Data_Healthy/cellchat/Diff_Number.pdf", width = 12, height = 12)
netVisual_heatmap(cellchat.Dog,  font.size = 14, 
                  color.heatmap = "OrRd", measure = "count",
                  cluster.rows = T, cluster.cols = T, #sources.use = 
                  remove.isolate = T,
                  title.name = "Differential number of interactions")
dev.off()

pdf("Data_Healthy/cellchat/Diff_Number.pdf", width = 12, height = 12)
netVisual_heatmap(cellchat.Dog,  font.size = 14, 
                  color.heatmap = "OrRd", measure = "count",
                  cluster.rows = F, cluster.cols = T, 
                  #    sources.use = c(39,50,24,20,36,19),
                  remove.isolate = T,
                  title.name = "Differential number of interactions")
dev.off()

pdf("Data_Healthy/cellchat/Diff_Stregnth.pdf", width = 12, height = 12) #6x4
netVisual_heatmap(cellchat.Dog,font.size = 14, color.heatmap = "OrRd",
                  cluster.rows= T, cluster.cols = T, 
                  #  targets.use = c(42,38,51,14,27),
                  remove.isolate = F,slot.name = "net",
                  measure = "weight",
                  title.name = "Differential interaction strength")
dev.off()

#
groupSize <- as.numeric(table(cellchat.Dog@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.Dog@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.Dog@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#
mat <- cellchat.Dog@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#
pathways.show <- c("FN1","SPP1","COLLAGEN","APP")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,10) # a numeric vector. 
netVisual_aggregate(cellchat.Dog, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.dog, signaling = pathways.show, layout = "circle")

#
pdf("Chord.pdf", width = 5, height = 5)
netVisual_aggregate(cellchat.Dog, signaling = pathways.show, layout = "chord", 
                    sources.use = c(10,20,24),
                    targets.use = c(6,11,13,19,23)
)
dev.off()
#
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.Dog, signaling = "CD6", color.heatmap = "Reds", remove.isolate = T,
                  cluster.cols = T)

#
pairLR.use <- extractEnrichedLR(cellchat.Dog, signaling = c("MIF","CD99"))
netVisual_bubble(cellchat.Dog, font.size = 12, pairLR.use = pairLR.use, 
                 sort.by.source = T, sort.by.target = T, 
                 #sources.use = c(39,50,24,20,36,19),
                 #targets.use = c(22,28),
                 targets.use = c(39,50,24,20,36,19),
                 angle.x = 45, thresh = 0.0001, 
                 remove.isolate = T)
ggsave("Data_Healthy/cellchat/L-R_Myeloid.png", width = 5.75, height = 4)

pairLR.use <- extractEnrichedLR(cellchat.Dog, signaling = c("MIF","FN1","CD99","APP"))
netVisual_bubble(cellchat.Dog, font.size = 12, #pairLR.use = pairLR.use, 
                 sort.by.source = T, sort.by.target = T, 
                 targets.use = c(19,36,47,21,26), #21,26,
                 sources.use = c(19,36,47), #21,26,
                 angle.x = 45, thresh = 0.001, grid.on = F,#min.quantile = 0.25,
                 remove.isolate = T)
ggsave("Data_Healthy/cellchat/L-R_Myeloid.png", width = 5.75, height = 9)

pairLR.use <- extractEnrichedLR(cellchat.Dog, signaling = c("VISFATIN","CD45"))
netVisual_bubble(cellchat.Dog, font.size = 12, #pairLR.use = pairLR.use, 
                 sort.by.source = T, sort.by.target = T, 
                 targets.use = c(1,2), # PMN
                 sources.use = c(-3,-4,-5,-6,-7,-8,-10,-16,-17), 
                 angle.x = 45, thresh = 0.001, grid.on = F,#min.quantile = 0.25,
                 remove.isolate = T)
ggsave("Data_Healthy/cellchat/L-R_CD4CD8.png", width = 7, height = 7)

pairLR.use <- extractEnrichedLR(cellchat.Dog, signaling = c("VISFATIN","CD45"))
netVisual_bubble(cellchat.Dog, font.size = 12, #pairLR.use = pairLR.use, 
                 sort.by.source = T, sort.by.target = T, 
                 sources.use = c(22), # PMN
                 # targets.use = c(-3,-4,-5,-6,-7,-8,-10,-16,-17), 
                 angle.x = 45, thresh = 0.001, grid.on = F,#min.quantile = 0.25,
                 remove.isolate = T)
ggsave("Data_Healthy/cellchat/L-R_CD4CD8.png", width = 7, height = 7)

pairLR.use <- extractEnrichedLR(cellchat.Dog, signaling = c("VISFATIN","CD45"))
netVisual_bubble(cellchat.Dog, font.size = 12, #pairLR.use = pairLR.use, 
                 sort.by.source = T, sort.by.target = T, 
                 targets.use = c(3,4,7,15), # PMN
                 sources.use = c(-1,-2,-5,-6,-8,-3,-4,-7,-15,-16),
                 angle.x = 45, thresh = 0.001, grid.on = F,#min.quantile = 0.25,
                 remove.isolate = T)
ggsave("Data_Healthy/cellchat/L-R_Minor.png", width = , height = 7)
#
pdf("incoming.pdf", width = 12, height = 10)
netAnalysis_signalingRole_heatmap(cellchat.Dog, pattern = "incoming", font.size = 10, cluster.cols = T,
                                  width = 12, height = 20)
dev.off()

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.Dog@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.Dog@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pdf("Data_Healthy/cellchat/Overall_Patterns.pdf", width = 20, height = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.Dog, pattern = "outgoing", height = 18, width = 20, font.size = 12) # , height = 16
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.Dog, pattern = "incoming", height = 18,width = 20, font.size = 12) # , height = 16
ht1 + ht2
dev.off()

pdf("Data_Healthy/cellchat/Scatter_In_Out_Patterns_label.pdf", width = 5, height = 5)
netAnalysis_signalingRole_scatter(cellchat.Dog, color.use = colorcode,
                                  label.size = 3, show.legend = T,
)
dev.off()
gg2 <- netAnalysis_signalingRole_scatter(cellchat.Dog, signaling = c("MIF"))

netVisual_heatmap(cellchat.Dog, signaling = c("MHC-II"), color.heatmap = "Reds", remove.isolate = T,cluster.rows = T, cluster.cols = T)

#
pdf("Data_Healthy/cellchat/MHC-II_network.pdf", width = 15, height = 10)
netAnalysis_signalingRole_network(cellchat.Dog, signaling = c("MHC-II"), 
                                  cluster.rows = T, cluster.cols = T,
                                  color.heatmap = "Blues",
                                  width = 30, height = 2, font.size = 10)
dev.off()

netAnalysis_contribution(cellchat.Dog, signaling = c("MHC-II"), thresh = 0.01, return.data = T, font.size = 12)
ggsave("Data_Healthy/cellchat/Contribution_MIF_CD99_MHCII.png", width = 3, height = 3)

plotGeneExpression(cellchat.Dog, signaling = c("CD99"),
                   enriched.only = T, type = "violin")

FeaturePlot(Dog.combined.singlet.myeloid, features = c("CEACAM1"),order = F,label = F, # ,"CD74","CD44","CD99","P"
            reduction = "umap", cols =  c("#ADDAE6","#E63222") ) + NoLegend()
ggsave("Data_Healthy/Myeloid/CEACAM1.png", width = 5, height = 5)

plot_density(Dog.combined.singlet.myeloid, 
             features = c("MX1"), size = 0.5, pal = "plasma",# "cividis", ""plasma", "magma"
             reduction = "umap") +NoLegend()
ggsave("Data_Healthy/Myeloid/MX1_density.png", width = 5, height = 5)


library(NMF)

library(ggalluvial)
selectK(cellchat.Dog, pattern = "outgoing")
nPatterns = 3

pdf("Data_Healthy/cellchat/patterns_outgoing.pdf", width = 12, height = 12)
cellchat.Dog <- identifyCommunicationPatterns(cellchat.Dog, pattern = "outgoing", k = nPatterns, 
                                              width = 3, height = 20, font.size = 12)
dev.off()

pdf("Data_Healthy/cellchat/River_outgoing.pdf", width = 12, height = 12)
netAnalysis_river(cellchat.Dog, pattern = "outgoing", 
                  color.use.pattern = c("#0432FF","#FF0000","black"), do.order = F)
dev.off()

selectK(cellchat.Dog, pattern = "incoming")
nPatterns = 3
pdf("Data_Healthy/cellchat/Patterns_incoming.pdf", width = 12, height = 12)
cellchat.Dog <- identifyCommunicationPatterns(cellchat.Dog, pattern = "incoming", k = 2, 
                                              width = 1, height = 20, font.size = 12)
dev.off()

netAnalysis_river(cellchat.Dog, pattern = "incoming", 
                  color.use.pattern = c("#0432FF","#FF0000"), do.order = F)
#
cc.genes <- Seurat::cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Dog.combined.singlet.CD8NKT <- CellCycleScoring(Dog.combined.singlet.CD8NKT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

freq_table <- Dog.combined.singlet.CD8NKT[[]]
freq_table <- freq_table[,c("Subset", "seurat_clusters", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(9,12))

freq_table <- freq_table %>%
  group_by(Subset, seurat_clusters) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$seurat_clusters <- factor(freq_table$seurat_clusters, levels = c(9,12))
ggplot(freq_table, aes(x=seurat_clusters, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() #+
geom_text(aes(label = percent),
          position = position_stack(vjust = .5), size = 4.5)
ggsave("PATH_TO_FILE", height=5, width=15)
ggsave("PATH_TO_FILE", height=6, width=20)
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD8NKT/Cluster_9_12_cellcycle.png", height=2, width=3)

FeatureScatter(Dog.combined.singlet.CD8NKT,jitter = T, log = T,
               feature1 = "T_Cell_Proliferation" , feature2 = "LAG3")

meta <- Dog.combined.singlet.CD8NKT[[]]
View(meta)
#
freq_table <- Dog.combined.singlet.CD8NKT[[]]
freq_table <- freq_table[,c("Subset", "old.ident", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided")
freq_table <- freq_table %>%
  group_by(Subset, old.ident, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$old.ident <- factor(freq_table$old.ident, levels = c("T_CD4","T_CD8"))

freq_table <- freq_table %>%
  group_by(Subset, old.ident) %>%
  mutate(sum = sum(n))
freq_table$percent <- round((freq_table$n/freq_table$sum)*100, 1)
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "S", "G2M"))
freq_table$old.ident <- factor(freq_table$old.ident, levels = c("T_CD4","T_CD8"))
ggplot(freq_table, aes(x=old.ident, y=percent, fill=Phase)) +
  geom_bar(stat="identity", color="black", lwd=0.25) +
  theme(axis.title.x = element_blank())+
  facet_grid(Subset ~.) +
  scale_fill_manual(values=c("#095096", "#F0642F", "#E8F0FE")) +    #BCABC1", "#B4D5D3", "#FEF4B7"     "#F67770", "#1BB941", "#649FFC"   /# # #
  theme_classic() +
  geom_text(aes(label = percent),
            position = position_stack(vjust = .5), size = 3)
ggsave("Data_Healthy/cellcycle.png", height=3, width=15)


#

1
install_github("velocyto-team/velocyto.R")

git clone https://github.com/velocyto-team/velocyto.R

library(DcjComm)
library(velocyto.R)

library(pcaMethods)
BiocManager::install("pcaMethods")
# setwd("PATH_TO_PROJECT")
devtools::install_local("velocyto.R")
library(slingshot, quietly = FALSE)
library(tradeSeq)
n
# Begins
sce <- as.SingleCellExperiment(Dog.combined.singlet.B, assay = "RNA")  
sce <- slingshot(sce, clusterLabels = 'seurat_clusters',reducedDim = "UMAP",
                 allow.breaks = FALSE)
sce$slingPseudotime_1
summary(sce$slingPseudotime_3)
lnes <- getLineages(SlingshotDataSet(sce),reducedDim(sce,"UMAP"),
                    sce$seurat_clusters)
lnes <- getLineages(SlingshotDataSet(sce), reducedDim(sce,"UMAP"),
                    sce$seurat_clusters, start.clus = "0")

plotSmoothers(sce, counts, gene = "IRF4")

dimred <- Dog.combined.singlet.B@reductions$lmds@cell.embeddings
clustering <- Dog.combined.singlet.B@meta.data$spliced_snn_res.0.25

set.seed(1)
lineages <- getLineages(dimred, clustering)


plot(reducedDims(sce)$UMAP,
     pch=16, 
     asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = c("black"))

slingshot_df <- colData(sce)[[]]
levels(colData(sce))

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)

rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(slingshot)
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2
plot(rd2, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')

summary(sce$slingPseudotime_1)

lnes <- getLineages(reducedDim(sce,"UMAP"), sce$seurat_clusters)

sce <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP",
                 allow.breaks = FALSE, start.clus="0")

lnes <- getLineages(reducedDim(sce,"UMAP"),
                    sce$seurat_clusters, start.clus = "0")
lnes <- getLineages(SlingshotDataSet(sce), reducedDim(sce,"UMAP"),
                    sce$seurat_clusters, start.clus = "0")

plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

plot(reducedDim(sce,"UMAP"), pch=16, asp = 3, col = brewer.pal(9,'Set1')[sce$seurat_clusters])
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

dev.off()
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)

my_color <- createPalette(length(levels(sce$cell_type)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$cell_type))

slingshot_df <- data.frame(colData(sce))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = cell_type, 
                         colour = cell_type)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

library(grDevices)

lnes <- getLineages(reducedDim(sce,"UMAP"),
                    sce$seurat_clusters)

sce
plot(reducedDims(sce)$PCA, col = my_color[as.character(sce$cell_type2)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(deng_SCE$cell_type2)]),  
       fill = my_color[levels(deng_SCE$cell_type2)])
lines(SlingshotDataSet(deng_SCE), lwd=2, type = 'lineages', col = c("black"))
library(S4Vectors)


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')


plotSmoothers(sce, assays(sce)$counts, gene = "CD79B")

library(tradeSeq)

set.seed(1)
lineages <- getLineages(SlingshotDataSet(sce), reducedDim(sce,"UMAP"),
                        sce$seurat_clusters, start.clus = "0")
lineages
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

sce <- fitGAM(sce)
ATres <- associationTest(sce)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]
View(heatdata)

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

heatmap(as.matrix(log1p(heatdata)), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

ColSideColors = brewer.pal(9,"Set1")(ncol(as.matrix(heatdata)))

sce$slingPseudotime_2
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#downloading_dataset
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(igraph)
  library(Seurat)
  library(cowplot)
})
BiocParallel::register(BiocParallel::SerialParam())
library(dplyr)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

ATres <- associationTest(sce)
ATres$fdr <- p.adjust(ATres$pvalue, method = "fdr")
ATres <- ATres[order(ATres$pvalue), ]
ATres$feature_id <- rownames(ATres)
feature_id <- ATres %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

pseudotime_start_end_association <- startVsEndTest(sce, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)
feature_id <- pseudotime_start_end_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)

branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)
feature_id <- branch_point_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)

#
seu <- Dog.combined.singlet.B
sds <- slingshot(Embeddings(seu, 'umap'), clusterLabels = seu$seurat_clusters,
                 start.clus = 0, stretch = 0)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

reducedDim()
cell_colors <- cell_pal(seu$seurat_clusters, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())

plot(reducedDim(sds), col = cell_colors_clust, 
     pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
SlingshotDataSet(sce)

nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}

#
table(sce$cell_type)
sce <- as.SingleCellExperiment(Dog.combined.singlet.B)
sce <- slingshot(sce, clusterLabels = cell_type, reducedDim = "PCA",
                 allow.breaks = FALSE)

#
BiocManager::install("destiny")
library(destiny)
library(Seurat)
library(dplyr)

seurat_dorsal <- subset(Dog.combined.singlet.CD8NKT, seurat_clusters %in% c(3,5,1,0,7,8,14,2,4,6))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 3000)
VariableFeatures(Dog.combined.singlet.CD8NKT) <- setdiff(VariableFeatures(Dog.combined.singlet.CD8NKT), unlist(cc.genes))
seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:30)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)
seurat_dorsal <- CellCycleScoring(seurat_dorsal,
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score", "G2M.Score"))
seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:30)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)

dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:30])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(-dpt$dpt) # if you find the constructed pseudotime goes to the wrong direction, 
# flip it (e.g. seurat_dorsal$dpt <- rank(-dpt$dpt))

FeaturePlot(seurat_dorsal, c("dpt","PDCD1","LAG3","TIGIT"), ncol=4)

seurat_dorsal$dpt <- max(seurat_dorsal$dpt) - seurat_dorsal$dpt # FLIP IT

if (is(seurat_dorsal[['RNA']], 'Assay5')){
  expr <- LayerData(seurat_dorsal, assay = "RNA", layer = "data")
} else{
  expr <- seurat_dorsal[['RNA']]@data
}

seurat_dorsal$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(seurat_dorsal[[]]), 
       aes(x = pseudotime_diffusionmap, 
           y = seurat_clusters, colour = seurat_clusters)) +
  # geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color)  + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")
library(ggplot2)
plot1 <- qplot(seurat_dorsal$dpt, as.numeric(expr["PDCD1",]),
               xlab="Dpt", ylab="Expression", main="PDCD1") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot2 <- qplot(seurat_dorsal$dpt, as.numeric(expr["LAG3",]),
               xlab="Dpt", ylab="Expression", main="LAG3") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot3 <- qplot(seurat_dorsal$dpt, as.numeric(expr["HAVCR2",]),
               xlab="Dpt", ylab="Expression", main="HAVCR2") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot1 + plot2 + plot3


#
PD1pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos"))
DimPlot(PD1pos, group.by = "PD1")
FeaturePlot(PD1pos, split.by = "PD1", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","CTLA4"), order = T)
ggsave("Data_Healthy/CD4/Featureplot_PDCD1_LAG3.png", width = 5, height = 5)

Idents(PD1pos) <- PD1pos$PD1_LAG3
marker <- FindMarkers(PD1pos, ident.1 = "Pos", ident.2 = "Neg", min.pct = 0.3)
write.table(marker,"Data_Healthy/LAG3specific_CD4_DEG.csv", append = F, sep = ",")

Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_LAG3)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD8/Cellnumber_PD1_LAG3_CD8.csv", sep = ",", append=F)

PD1pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_CTLA4
FeaturePlot(PD1pos, features = c("PDCD1"))
FeaturePlot(PD1pos, split.by = "PD1_CTLA4", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","CTLA4"), order = T)
Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_CTLA4)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD8/Cellnumber_PD1_LAG3_CD8.csv", sep = ",", append=F)

PD1pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_TIGIT
FeaturePlot(PD1pos, split.by = "PD1_CTLA4", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","CTLA4"), order = T)
Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_TIGIT)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD8/Cellnumber_PD1_LAG3_CD8.csv", sep = ",", append=F)

PD1pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_TIM3
FeaturePlot(PD1pos, split.by = "PD1_CTLA4", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","CTLA4"), order = T)
Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_TIM3)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD8/Cellnumber_PD1_LAG3_CD8.csv", sep = ",", append=F)

PD1pos <- subset(Dog.combined.singlet.CD8NKT, PD1 %in% c("Pos"))
Idents(PD1pos) <- PD1pos$PD1_TIGIT
FeaturePlot(PD1pos, split.by = "PD1_TIGIT", cols =  c("#ADDAE6","#E63222"),
            features = c("PDCD1","TIGIT"), order = T)
Cellnumbers.subset <- table(Idents(PD1pos), PD1pos$PD1_TIGIT)
write.table(Cellnumbers.subset,
            "Data_Healthy/CD8/Cellnumber_PD1_LAG3_CD8.csv", sep = ",", append=F)


# What would be the C13 in CD8NKT cells
treg_genes <- c("FOXP3", "IL2RA", "IKZF2", "CTLA4")
exhausted_genes <- c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX")
early_activation_genes <- c("CD69", "IL2", "NR4A1", "EGR1", "FOS")

Dog.combined.singlet.CD8NKT <- AddModuleScore(
  object = Dog.combined.singlet.CD8NKT,
  features = list(treg_genes, exhausted_genes, early_activation_genes),
  name = c("Treg_Score", "Exhausted_Score", "Activation_Score")
)
FeaturePlot(Dog.combined.singlet.CD8NKT, features = c("Treg_Score1", "Exhausted_Score2", "Activation_Score3"), 
            min.cutoff = "q10", max.cutoff = "q90", order = T)
DotPlot(Dog.combined.singlet.CD8NKT, 
        features = c("Treg_Score1", "Exhausted_Score2", "Activation_Score3"))
ggsave("Data_Healthy/Featureplot_Cluster_Marker/CD413/Scores.png", width = 5, height = 10)

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(RColorBrewer)
library(viridis)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce <- as.SingleCellExperiment(Dog.combined.singlet.CD8NKT)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
reducedDims(sce)$UMAP <- Dog.combined.singlet.CD8NKT@reductions$umap@cell.embeddings

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce$cluster <- as.character(Dog.combined.singlet.CD8NKT@meta.data$seurat_clusters)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
summary(sce$slingPseudotime_1)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(reducedDims(sce)$UMAP,
     col = viridis(100)[cut(sce$slingPseudotime_1, breaks = 100)],
     pch = 16, asp = 1,
     main = "Slingshot Pseudotime (Trajectory 1)")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
text(reducedDims(sce)$UMAP, labels = sce$cluster, cex = 0.6, col = "black")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce <- slingshot(
  sce,
  clusterLabels = 'cluster',
  reducedDim = 'UMAP',
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
)

plot(reducedDims(sce)$UMAP,
     col = viridis::viridis(100)[cut(sce$slingPseudotime_1, breaks = 100)],
     pch = 16, asp = 1,
     main = "Pseudotime from Cluster 3 and 5 (Naive T cells)")
text(reducedDims(sce)$UMAP, labels = sce$cluster, cex = 0.6, col = "black")

plot(sce$slingPseudotime_1, 
     log1p(assay(sce)["CTLA4",]), 
     col = "darkred", pch = 16,
     xlab = "Pseudotime", ylab = "CTLA4 expression")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(reducedDims(sce)$UMAP, 
     col = "lightgray", pch = 16, asp = 1,
     main = "Pseudotime from Cluster 3 and 5 (Naive T cells)")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
is_valid <- !is.na(sce$slingPseudotime_1)
points(reducedDims(sce)$UMAP[is_valid, ], 
       col = viridis::viridis(100)[cut(sce$slingPseudotime_1[is_valid], breaks = 100)], 
       pch = 16)


library(slingshot)
library(viridis)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
is_valid <- !is.na(sce$slingPseudotime_1)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
point.size <- 0.8

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(reducedDims(sce)$UMAP,
     col = "lightgray", 
     pch = 16, cex = point.size,
     asp = 1,
     main = "Pseudotime with Trajectory Curve")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
points(reducedDims(sce)$UMAP[is_valid, ],
       col = viridis(100)[cut(sce$slingPseudotime_1[is_valid], breaks = 100)],
       pch = 16, cex = point.size)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
lines(SlingshotDataSet(sce), lwd = 2, col = "black", type = 'line')


# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
curve_data <- SlingshotDataSet(sce)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
start_nodes <- lapply(curve_data@curves, function(crv) crv$s[1, ])
end_nodes <- lapply(curve_data@curves, function(crv) crv$s[nrow(crv$s), ])

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(reducedDims(sce)$UMAP,
     col = "lightgray", pch = 16, cex = 0.8, asp = 1,
     main = "Trajectory with Start (green) and End (red)")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
points(reducedDims(sce)$UMAP[is_valid, ],
       col = viridis(100)[cut(sce$slingPseudotime_1[is_valid], breaks = 100)],
       pch = 16, cex = 0.8)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
lines(curve_data, lwd = 2, col = "black", type = 'line')

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (pt in start_nodes) {
  points(pt[1], pt[2], col = "forestgreen", pch = 16, cex = 2)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (pt in end_nodes) {
  points(pt[1], pt[2], col = "red", pch = 16, cex = 2)
}

plot(reducedDims(sce)$UMAP,
     col = "lightgray", 
     # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
     asp = 1,
     main = "Pseudotime with Trajectory Curve")

points(reducedDims(sce)$UMAP[is_valid, ],
       col = viridis::viridis(100)[cut(sce$slingPseudotime_1[is_valid], breaks = 100)],
       pch = 16, cex = 0.4)

lines(SlingshotDataSet(sce), lwd = 2, col = "black", type = 'line')

curves <- SlingshotDataSet(sce)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_nodes <- lapply(curves@curves, function(crv) crv$s[nrow(crv$s), ])
length(end_nodes)

library(slingshot)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
curves <- SlingshotDataSet(sce)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_coords <- lapply(curves@curves, function(crv) crv$s[nrow(crv$s), ])

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
umap_coords <- reducedDims(sce)$UMAP

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_cells <- sapply(end_coords, function(pt) {
  dists <- sqrt((umap_coords[,1] - pt[1])^2 + (umap_coords[,2] - pt[2])^2)
  which.min(dists)
})

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_clusters <- sce$cluster[end_cells]
names(end_clusters) <- paste0("End_", seq_along(end_clusters))
end_clusters

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cells_in_13 <- which(sce$cluster == "13")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
included_in_curves <- sapply(SlingshotDataSet(sce)@curves, function(crv) {
  any(rownames(reducedDims(sce))[cells_in_13] %in% rownames(crv$s))
})
included_in_curves

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce2 <- slingshot(
  sce,
  clusterLabels = 'cluster',
  reducedDim = 'UMAP',
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
)

included_in_curves <- sapply(SlingshotDataSet(sce)@curves, function(crv) {
  any(rownames(reducedDims(sce))[cells_in_13] %in% rownames(crv$s))
})
included_in_curves


#
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
subset_13 <- sce[, sce$cluster == "13"]

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
subset_13$cluster <- "13"

subset_13 <- slingshot(
  subset_13,
  clusterLabels = 'cluster',
  reducedDim = 'UMAP'
)

plot(reducedDims(subset_13)$UMAP,
     col = viridis(100)[cut(subset_13$slingPseudotime_1, breaks = 100)],
     pch = 16, asp = 1,
     main = "Pseudotime within Cluster 13")
lines(SlingshotDataSet(subset_13), lwd = 2, col = "black")

#
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
length(end_nodes)
library(ggplot2)
library(mgcv)

df <- data.frame(
  pseudotime = sce$slingPseudotime_1,
  CTLA4_expr = log1p(assay(sce)["CTLA4", ])
)

ggplot(df, aes(x = pseudotime, y = CTLA4_expr)) +
  geom_point(color = "darkred", alpha = 0.6) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", se = TRUE) +
  labs(x = "Pseudotime", y = "CTLA4 expression (log1p)", 
       title = "CTLA4 expression along pseudotime") +
  theme_minimal()

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

plot(reducedDims(sce)$UMAP, 
     col = ifelse(CTLA4_high, "firebrick", "lightgray"), 
     pch = 16, main = "CTLA4-high cells on UMAP")

genes <- c("FOXP3", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4")
par(mfrow = c(3, 2))
for (g in genes) {
  plot(sce$slingPseudotime_1,
       log1p(assay(sce)[g, ]),
       col = "darkblue", pch = 16,
       xlab = "Pseudotime", ylab = paste0(g, " expression"))
}
library(dplyr)
library(ggplot2)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
df <- data.frame(
  pseudotime = sce$slingPseudotime_1,
  CTLA4 = log1p(assay(sce)["CTLA4", ])
)

#
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cells_13 <- which(sce$cluster == "13")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ctla4_expr_13 <- log1p(assay(sce)["CTLA4", cells_13])

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
pct_ctla4_13 <- sum(ctla4_expr_13 > 0) / length(ctla4_expr_13) * 100

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

#
DotPlot(Dog.combined.singlet.CD8NKT, features = c("CTLA4"))
DimPlot(Dog.combined.singlet.CD8NKT)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
target_clusters <- c("9", "12", "13")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ctla4_stats <- data.frame(
  Cluster = character(),
  Percent_Expressed = numeric(),
  Average_Expression = numeric(),
  stringsAsFactors = FALSE
)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (clus in target_clusters) {
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  cells <- which(sce$cluster == clus)
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  expr_vals <- log1p(assay(sce)["CTLA4", cells])
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  pct_expressed <- sum(expr_vals > 0) / length(expr_vals) * 100
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  avg_expr <- mean(expr_vals)
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  ctla4_stats <- rbind(ctla4_stats, data.frame(
    Cluster = clus,
    Percent_Expressed = round(pct_expressed, 1),
    Average_Expression = round(avg_expr, 3)
  ))
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
print(ctla4_stats)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
df_summary <- df %>%
  filter(!is.na(pseudotime)) %>%
  mutate(pt_bin = cut(pseudotime, breaks = 50)) %>%
  group_by(pt_bin) %>%
  summarise(
    pt_mid = mean(pseudotime, na.rm = TRUE),
    CTLA4_mean = mean(CTLA4, na.rm = TRUE)
  )

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ggplot(df_summary, aes(x = pt_mid, y = CTLA4_mean)) +
  geom_line(color = "darkred", size = 1.2) +
  labs(x = "Pseudotime", y = "Mean CTLA4 expression (log1p)", 
       title = "Average CTLA4 expression over pseudotime") +
  theme_minimal()

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
df_box <- df %>%
  filter(!is.na(pseudotime)) %>%
  mutate(pt_group = cut(pseudotime, breaks = 6, labels = paste0("Bin", 1:6)))

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ggplot(df_box, aes(x = pt_group, y = CTLA4)) +
  geom_boxplot(fill = "firebrick", alpha = 0.7, outlier.color = "black", outlier.size = 1) +
  labs(x = "Pseudotime bin", y = "CTLA4 expression (log1p)", 
       title = "CTLA4 expression distribution by pseudotime bin") +
  theme_minimal()


pdf("Data_Healthy/Featureplot_Cluster_Marker/CD413/Figure_Pseudotime_Trajectory.pdf", width = 7, height = 6)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

pt_vals <- sce$slingPseudotime_1
pt_colors <- ifelse(!is.na(pt_vals),
                    viridis(100)[cut(pt_vals, breaks = 100)],
                    "lightgray")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(reducedDims(sce)$UMAP,
     col = ifelse(!is.na(sce$slingPseudotime_1),
                  viridis(100)[cut(sce$slingPseudotime_1, breaks = 100)],
                  "lightgray"),
     pch = 16, cex = 0.4,
     asp = 1,
     main = "Pseudotime with Start/End Points and Cluster 13")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cells_13 <- which(sce$cluster == "13")
points(reducedDims(sce)$UMAP[cells_13, ],
       col = "#7E96FF", pch = 1, lwd = 1.2, cex = 0.6)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
curves <- SlingshotDataSet(sce)@curves
for (i in seq_along(curves)) {
  pts <- curves[[i]]$s
  points(pts[,1], pts[,2], col = "gray70", pch = 16, cex = 0.2)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
start_points <- lapply(curves, function(crv) crv$s[1, ])
for (pt in start_points) {
  points(pt[1], pt[2], col = "blue", pch = 16, cex = 1)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_points <- lapply(curves, function(crv) crv$s[nrow(crv$s), ])
for (pt in end_points) {
  points(pt[1], pt[2], col = "firebrick", pch = 16, cex = 1)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
image.plot(legend.only = TRUE,
           zlim = range(pt_vals, na.rm = TRUE),
           col = viridis(100),
           legend.args = list(text = "Pseudotime", side = 3, line = 1))

dev.off()

end_points <- lapply(curves, function(crv) crv$s[nrow(crv$s), ])
for (pt in end_points) {
  points(pt[1], pt[2], col = "firebrick", pch = 16, cex = 1)
}



# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
curves <- SlingshotDataSet(sce)@curves
umap_coords <- reducedDims(sce)$UMAP

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_cells <- sapply(curves, function(crv) {
  pt <- crv$s[nrow(crv$s), ]
  dists <- sqrt((umap_coords[,1] - pt[1])^2 + (umap_coords[,2] - pt[2])^2)
  which.min(dists)
})

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_coords <- umap_coords[end_cells, ]
end_clusters <- sce$cluster[end_cells]

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
points(end_coords[,1], end_coords[,2], col = "firebrick", pch = 16, cex = 1.6)
text(end_coords[,1], end_coords[,2],
     labels = paste0("Clus ", end_clusters),
     col = "firebrick", pos = 3, cex = 0.7)

library(ggplot2)
library(tidyr)
library(dplyr)

pdf("Data_Healthy/Featureplot_Cluster_Marker/CD413/Figure_Pseudotime_CTLA4.pdf", width = 7, height = 6)

# long-format pseudotime + expression matrix
df <- data.frame(
  cell = colnames(sce),
  cluster = sce$cluster,
  CTLA4 = log1p(assay(sce)["CTLA4", ])
)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
lineage_cols <- grep("^slingPseudotime_", colnames(colData(sce)), value = TRUE)
df <- cbind(df, as.data.frame(colData(sce)[, lineage_cols]))

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
df_long <- df %>%
  pivot_longer(
    cols = all_of(lineage_cols),
    names_to = "lineage",
    values_to = "pseudotime"
  ) %>%
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  ggplot(df_long, aes(x = pseudotime, y = CTLA4)) +
  geom_point(alpha = 0.3, size = 0.4, color = "darkred") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, color = "black") +
  facet_wrap(~ lineage, scales = "free_x") +
  labs(title = "CTLA4 expression across pseudotime by lineage",
       x = "Pseudotime", y = "CTLA4 expression (log1p)") +
  theme_minimal(base_size = 12)

dev.off()

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
curves <- SlingshotDataSet(sce)@curves
umap_coords <- reducedDims(sce)$UMAP

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_cells <- sapply(curves, function(crv) {
  pt <- crv$s[nrow(crv$s), ]
  dists <- sqrt((umap_coords[,1] - pt[1])^2 + (umap_coords[,2] - pt[2])^2)
  which.min(dists)
})

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
end_clusters <- sce$cluster[end_cells]

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
lineage_map <- data.frame(
  Lineage = paste0("slingPseudotime_", seq_along(end_cells)),
  Terminal_Cell = colnames(sce)[end_cells],
  Terminal_Cluster = end_clusters
)

print(lineage_map)


#
install.packages("clustree")
library(clustree)

DefaultAssay(Dog.combined.singlet) <- "integrated"
Dog.combined.singlet <-ScaleData(Dog.combined.singlet, verbose = FALSE)
Dog.combined.singlet <-RunPCA(Dog.combined.singlet, npcs = 30, verbose = FALSE)
Dog.combined.singlet <-RunUMAP(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- RunTSNE(Dog.combined.singlet, dims = 1:30, reduction = "pca")
Dog.combined.singlet <- FindNeighbors(Dog.combined.singlet, reduction = "pca", dims = 1:30)
Dog.combined.singlet <- FindClusters(Dog.combined.singlet,
                                     resolution = c(0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4))
DefaultAssay(Dog.combined.singlet) <- "RNA"
clustree(Dog.combined.singlet, node_size = 3.5,
         prefix = "integrated_snn_res.")
ggsave("Data_Healthy/dendrogram.png", width = 12, height = 7)

FeaturePlot(Dog.combined.singlet.CD8NKT, order = T, features = c("CTLA4","FOXP3","IL10", "GATA3","IL2RA"))


# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
features <- c("CTLA4", "FOXP3", "IL10", "GATA3", "IL2RA")
cells_13 <- WhichCells(Dog.combined.singlet.CD8NKT, idents = "13")
expr_percent_cluster13 <- sapply(features, function(gene) {
  expr <- FetchData(Dog.combined.singlet.CD8NKT, vars = gene, cells = cells_13)
  percent <- sum(expr[, 1] > 0) / length(cells_13) * 100
  return(round(percent, 2))
})
expr_percent_cluster13







# setwd("PATH_TO_PROJECT")# age
unique(Dog.combined.singlet$Subset5) 


#

#


Dog.combined.singlet$Age_num <- as.numeric(gsub("y", "", Dog.combined.singlet$Subset5))
unique(Dog.combined.singlet$Age_num)

DimPlot(Dog.combined.singlet, group.by="Age_num", cols=viridis::viridis(5), shuffle = T)

tab <- table(Dog.combined.singlet$Age_num, Dog.combined.singlet$seurat_clusters)
prop <- prop.table(tab, 1)
pdf("Age_Cell_Composition.pdf", width =7, height = 2)
pheatmap::pheatmap(prop, cluster_rows=FALSE, cluster_cols=TRUE)
dev.off()


cd4 <- subset(Dog.combined.singlet, idents = c(1,24,16,32,22,6,27,25,39,48,43))
DefaultAssay(cd4) <- "RNA"
mat  <- GetAssayData(cd4, slot = "data")
age  <- cd4$Age_num
fit_list <- apply(mat, 1, function(y) {
  co <- tryCatch(summary(lm(y ~ age))$coefficients, error = function(e) NA)
  if (is.matrix(co)) c(beta = co["age","Estimate"], p = co["age","Pr(>|t|)"]) else c(beta = NA, p = NA)
})
df <- as.data.frame(t(fit_list))
df$padj <- p.adjust(df$p, method = "BH")
df <- df[order(df$padj), ]
head(df)

hist(df$beta, breaks=60); summary(df$beta[is.finite(df$beta)])

DefaultAssay(cd4) <- "RNA"
mat <- GetAssayData(cd4, slot="data")
age <- cd4$Age_num
umi <- cd4$nCount_RNA
mt  <- cd4$percent.mt

fit <- apply(mat, 1, function(y){
  co <- tryCatch(summary(lm(y ~ age + umi + mt))$coefficients, error=function(e) NA)
  if (is.matrix(co)) c(beta_age = co["age","Estimate"], p_age = co["age","Pr(>|t|)"])
  else c(beta_age=NA,p_age=NA)
})

df_adj <- as.data.frame(t(fit))
df_adj$padj <- p.adjust(df_adj$p_age, "BH")

sig <- subset(df_adj, abs(beta_age) >= 0.15 & padj < 0.05)
nrow(sig)

library(Matrix); library(edgeR); library(limma)


###
library(Seurat)
library(Matrix)
library(edgeR)
library(limma)

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
cd8 <- subset(
  Dog.combined.singlet,
  idents = c(29,14,17,30,20,15,18,26,
             43,32,8,23,
             34,42,41,46,
             16,27)
)

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
cd8$Age_num <- as.numeric(gsub("y", "", cd8$Subset5))
unique(cd8$Age_num)

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
cd8$SampleID <- cd8$orig.ident
unique(cd8$SampleID)

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
cts  <- GetAssayData(cd8, slot="counts")
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
pb_counts <- t(sapply(
  split(seq_len(ncol(cts)), samp),
  function(ix) Matrix::rowSums(cts[, ix, drop=FALSE])
))

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
meta <- unique(cd8@meta.data[, c("SampleID", "Age_num")])
rownames(meta) <- meta$SampleID
meta <- meta[rownames(pb_counts), , drop=FALSE]

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
y <- DGEList(counts = t(pb_counts))
y <- calcNormFactors(y)

design <- model.matrix(~ Age_num, data = meta)

v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res_cd8 <- topTable(fit, coef="Age_num", number=Inf, sort.by="P")

#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
#########################################################
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sig_cd8 <- subset(res_cd8, adj.P.Val < 0.05 & abs(logFC) > log2(1.2))
sig_cd8


#
unique(Dog.combined.singlet$Subset6)
Dog.combined.singlet$Breed <- Dog.combined.singlet$Subset6   # "Mixed", "Maltese", "Poodle"
table(Dog.combined.singlet$Breed)

clusters_list <- list(
  CD4     = c(25,13,6,19,22,24,3,0,11,10,21,38,9,2,12),
  CD8     = c(1,24,16,32,22,6,27,25,39,48,43),
  Myeloid = c(29,14,17,30,20,15,18,26,43,32,8,23,34,42,41,46,16,27),
  B       = c(13,26,41,37,50)
)

library(Matrix); library(edgeR); library(limma)

library(Matrix)
library(edgeR)
library(limma)

pseudobulk_breed <- function(obj, clusters) {
  
  sub <- subset(obj, idents = clusters)
  
  counts <- GetAssayData(sub, slot="counts")
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  breed  <- sub$Breed       # Maltese / Mixed / Poodle
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  pb <- t(sapply(split(seq_len(ncol(counts)), samp),
                 function(ix) Matrix::rowSums(counts[, ix, drop=FALSE])))
  
  # metadata
  meta <- unique(sub@meta.data[, c("orig.ident","Breed")])
  rownames(meta) <- meta$orig.ident
  meta <- meta[rownames(pb), , drop=FALSE]
  
  meta$Breed <- factor(meta$Breed, levels=c("Mixed","Maltese","Poodle"))
  
  y <- DGEList(counts = t(pb))
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ 0 + Breed, data=meta)
  colnames(design) <- levels(meta$Breed)
  
  v   <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)
  
  L <- makeContrasts(
    Maltese_vs_Mixed = Maltese - Mixed,
    Poodle_vs_Mixed  = Poodle  - Mixed,
    Maltese_vs_Poodle= Maltese - Poodle,
    levels=design
  )
  
  fit2 <- contrasts.fit(fit, L)
  fit2 <- eBayes(fit2)
  
  list(
    Maltese_vs_Mixed = topTable(fit2, coef="Maltese_vs_Mixed", number=Inf),
    Poodle_vs_Mixed  = topTable(fit2, coef="Poodle_vs_Mixed", number=Inf),
    Maltese_vs_Poodle= topTable(fit2, coef="Maltese_vs_Poodle", number=Inf)
  )
}


results <- lapply(names(clusters_list), function(nm) {
  pseudobulk_breed(Dog.combined.singlet, clusters_list[[nm]])
})
names(results) <- names(clusters_list)

summary_df <- do.call(rbind, lapply(names(results), function(comp) {
  res <- results[[comp]]
  
  data.frame(
    Compartment = comp,
    Contrast    = c("Maltese vs Mixed", "Poodle vs Mixed", "Maltese vs Poodle"),
    FDR_0.05    = c(
      sum(res$Maltese_vs_Mixed$adj.P.Val < 0.05),
      sum(res$Poodle_vs_Mixed$adj.P.Val  < 0.05),
      sum(res$Maltese_vs_Poodle$adj.P.Val < 0.05)
    ),
    FDR_Eff12   = c(
      sum(res$Maltese_vs_Mixed$adj.P.Val < 0.05 & abs(res$Maltese_vs_Mixed$logFC) > log2(1.2)),
      sum(res$Poodle_vs_Mixed$adj.P.Val  < 0.05 & abs(res$Poodle_vs_Mixed$logFC)  > log2(1.2)),
      sum(res$Maltese_vs_Poodle$adj.P.Val < 0.05 & abs(res$Maltese_vs_Poodle$logFC) > log2(1.2))
    )
  )
}))

library(Matrix)
library(edgeR)
library(limma)

pseudobulk_breed_fixed <- function(obj, clusters, covars = c("Age_num")) {
  sub <- subset(obj, idents = clusters)
  counts <- GetAssayData(sub, slot="counts")
  samp   <- sub$orig.ident
  breed  <- sub$Breed
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  pb <- t(sapply(split(seq_len(ncol(counts)), samp),
                 function(ix) Matrix::rowSums(counts[, ix, drop=FALSE])))
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  md <- unique(sub@meta.data[, c("orig.ident","Breed", covars[covars %in% colnames(sub@meta.data)])])
  rownames(md) <- md$orig.ident
  md <- md[rownames(pb), , drop=FALSE]
  md$Breed <- factor(md$Breed, levels = c("Mixed","Maltese","Poodle"))
  
  # 3) DGEList
  y <- DGEList(counts = t(pb))
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  y <- calcNormFactors(y)
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  form <- as.formula(paste("~ 0 + Breed", 
                           if (any(colnames(md) %in% covars)) paste("+", paste(covars[covars %in% colnames(md)], collapse="+")) else "" ))
  design <- model.matrix(form, data = md)
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  v <- voomWithQualityWeights(y, design = design, plot = FALSE)
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  fit <- lmFit(v, design)
  L <- makeContrasts(
    Maltese_vs_Mixed  = Maltese - Mixed,
    Poodle_vs_Mixed   = Poodle  - Mixed,
    Maltese_vs_Poodle = Maltese - Poodle,
    levels = design
  )
  fit2 <- contrasts.fit(fit, L)
  fit2 <- eBayes(fit2)
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  get_sum <- function(tbl) {
    data.frame(
      n_FDR0.05 = sum(tbl$adj.P.Val < 0.05),
      n_FDR_Eff12 = sum(tbl$adj.P.Val < 0.05 & abs(tbl$logFC) > log2(1.2))
    )
  }
  res <- list(
    Maltese_vs_Mixed  = topTable(fit2, coef="Maltese_vs_Mixed",  number=Inf),
    Poodle_vs_Mixed   = topTable(fit2, coef="Poodle_vs_Mixed",   number=Inf),
    Maltese_vs_Poodle = topTable(fit2, coef="Maltese_vs_Poodle", number=Inf),
    md = md, v = v
  )
  res$summary <- rbind(
    cbind(Contrast="Maltese vs Mixed",  get_sum(res$Maltese_vs_Mixed)),
    cbind(Contrast="Poodle vs Mixed",   get_sum(res$Poodle_vs_Mixed)),
    cbind(Contrast="Maltese vs Poodle", get_sum(res$Maltese_vs_Poodle))
  )
  res
}

clusters_list <- list(
  CD4     = c(25,13,6,19,22,24,3,0,11,10,21,38,9,2,12),
  CD8     = c(1,24,16,32,22,6,27,25,39,48,43),
  Myeloid = c(29,14,17,30,20,15,18,26,43,32,8,23,34,42,41,46,16,27),
  B       = c(13,26,41,37,50)
)

res_all <- lapply(names(clusters_list), function(nm) {
  out <- pseudobulk_breed_fixed(Dog.combined.singlet, clusters_list[[nm]], covars = c("Age_num"))
  attr(out, "name") <- nm
  out
})

summary_df2 <- do.call(rbind, lapply(res_all, function(x){
  cbind(Compartment = attr(x,"name"), x$summary)
}))
summary_df2

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
b_idx <- which(names(res_all) == "B")
b_res <- res_all[[b_idx]]

res_all
length(res_all)
str(res_all)
names(res_all)

b_res <- res_all[[4]]
names(b_res)
b_sig <- lapply(b_res, function(x){
  subset(x, adj.P.Val < 0.05)
})

sig_B_MvMix <- subset(
  b_res$Maltese_vs_Mixed,
  adj.P.Val < 0.05 & abs(logFC) > log2(1.2)
)

head(sig_B_MvMix[, c("logFC","adj.P.Val","AveExpr")])

loo_counts <- lapply(unique(b_res$md$orig.ident), function(drop_id){
  md_loo <- subset(b_res$md, orig.ident != drop_id)
  keep_ids <- md_loo$orig.ident
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
})

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]



# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
obj <- Dog.combined.singlet
DefaultAssay(obj) <- "RNA"

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cd14_cells <- WhichCells(obj, expression = CD14 > 0)
obj_cd14 <- subset(obj, cells = cd14_cells)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
triple_pos <- WhichCells(obj_cd14, expression = S100A4 > 0 & CCL23 > 0 & TNFSF13 > 0)
triple_neg <- WhichCells(obj_cd14, expression = S100A4 == 0 & CCL23 == 0 & TNFSF13 == 0)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
obj_cd14$TripleGroup <- "OtherCD14"
obj_cd14$TripleGroup[colnames(obj_cd14) %in% triple_pos] <- "TriplePositive"
obj_cd14$TripleGroup[colnames(obj_cd14) %in% triple_neg] <- "TripleNegative"

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
table(obj_cd14$TripleGroup)

Idents(obj_cd14) <- "TripleGroup"

deg_triplepos_vs_tripleneg <- FindMarkers(
  obj_cd14,
  ident.1 = "TriplePositive",
  ident.2 = "TripleNegative",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0,
  test.use = "wilcox"
)
write.csv(deg_triplepos_vs_tripleneg, "DEG_CD14_TriplePos_vs_TripleNeg.csv")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
DimPlot(obj_cd14, reduction = "umap", group.by = "TripleGroup", pt.size = 0.5) + NoLegend()
ggsave("umap_myeloid.png", width = 5, height = 5)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
FeaturePlot(obj_cd14, features = c("S100A4","CCL23","TNFSF13","CD14"),
            order = TRUE, reduction = "umap")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
obj_cd14_sub <- subset(obj_cd14, cells = colnames(obj_cd14)[obj_cd14$TripleGroup %in% c("TriplePositive","TripleNegative")])
Idents(obj_cd14_sub) <- "TripleGroup"
table(obj_cd14_sub$TripleGroup)

DimPlot(obj_cd14_sub, cols = c(c("#ADDAE6","#E63222")), shuffle = T) + NoLegend()
ggsave("UMAP_TP_TN.png", width = 5, height = 5)

deg_direct <- FindMarkers(
  obj_cd14_sub,
  ident.1 = "TriplePositive",
  ident.2 = "TripleNegative",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0,
  test.use = "wilcox"
)
write.csv(deg_direct, "DEG_CD14_TriplePos_vs_TripleNeg_direct.csv")

set.seed(123)
n_neg <- sum(Idents(obj_cd14_sub) == "TripleNegative")  # 189
pos_cells <- WhichCells(obj_cd14_sub, idents = "TriplePositive")
neg_cells <- WhichCells(obj_cd14_sub, idents = "TripleNegative")

pos_down <- sample(pos_cells, size = n_neg, replace = FALSE)
obj_ds <- subset(obj_cd14_sub, cells = c(pos_down, neg_cells))
Idents(obj_ds) <- "TripleGroup"
table(Idents(obj_ds))

deg_downsampled <- FindMarkers(
  obj_ds,
  ident.1 = "TriplePositive",
  ident.2 = "TripleNegative",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0,
  test.use = "wilcox"
)
write.csv(deg_downsampled, "DEG_CD14_TriplePos_vs_TripleNeg_downsampled.csv")

library(Matrix)
DefaultAssay(obj_cd14_sub) <- "RNA"
counts <- GetAssayData(obj_cd14_sub, slot = "counts")

meta <- obj_cd14_sub@meta.data[, c("orig.ident","TripleGroup")]
meta$cell <- rownames(meta)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
meta$sample_group <- paste(meta$orig.ident, meta$TripleGroup, sep = "__")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
pb_list <- lapply(split(meta$cell, meta$sample_group), function(cells) {
  Matrix::rowSums(counts[, cells, drop = FALSE])
})
pb_mat <- do.call(cbind, pb_list)  # genes x (sample_group)
colnames(pb_mat) <- names(pb_list)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
design <- do.call(rbind, strsplit(colnames(pb_mat), "__"))
design <- data.frame(sample = design[,1], group = design[,2])
rownames(design) <- colnames(pb_mat)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
keep_samples <- names(which(table(design$sample) >= 2))
keep_cols <- rownames(design)[design$sample %in% keep_samples]
pb_mat <- pb_mat[, keep_cols]
design <- design[keep_cols, , drop = FALSE]

library(edgeR)
dge <- DGEList(counts = pb_mat, group = design$group)
dge <- calcNormFactors(dge, method = "TMM")

design_matrix <- model.matrix(~ 0 + group, data = design)
colnames(design_matrix) <- gsub("group", "", colnames(design_matrix))

dge <- estimateDisp(dge, design_matrix)
fit <- glmQLFit(dge, design_matrix)

# TriplePositive vs TripleNegative
contr <- makeContrasts(TriplePositive - TripleNegative, levels = design_matrix)
qlf <- glmQLFTest(fit, contrast = contr)
res_pb <- topTags(qlf, n = Inf)$table
write.csv(res_pb, "DEG_pseudobulk_TriplePos_vs_TripleNeg_edgeR.csv")

DimPlot(obj_cd14_sub, reduction="umap", group.by="TripleGroup", pt.size=0.5)
prop.table(table(obj_cd14_sub$orig.ident, obj_cd14_sub$TripleGroup), 1)


tnf_genes <- c("NFKB1","RELA","NFKBIA","TNFAIP3","JUNB","FOS")
chemokine_genes <- c("CCL2","CCL3","CCL4","CCL7","CCL20","CXCL8","CCL23")
ecm_genes <- c("S100A4","MMP9","ITGA4","ITGB2","LGALS3")

obj_cd14_sub <- AddModuleScore(obj_cd14_sub, features = list(tnf_genes), name = "TNF_Score")
obj_cd14_sub <- AddModuleScore(obj_cd14_sub, features = list(chemokine_genes), name = "Chemokine_Score")
obj_cd14_sub <- AddModuleScore(obj_cd14_sub, features = list(ecm_genes), name = "ECM_Score")

VlnPlot(obj_cd14_sub, features = c("TNF_Score1","Chemokine_Score1","ECM_Score1"),
        group.by="TripleGroup", pt.size=0.1, cols = c("#ADDAE6","#E63222"))
ggsave("TP_TN_TNF_Chemokine_ECM_score.png", width = 5, height = 4)

library(ggpubr)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
features <- c("TNF_Score1","Chemokine_Score1","ECM_Score1")

for (f in features) {
  p <- VlnPlot(obj_cd14_sub, features = f, group.by = "TripleGroup", pt.size = 0) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif", 
                       ref.group = "TripleNegative") + 
    ggtitle(paste0(f, " (TriplePositive vs TripleNegative)"))
  print(p)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
features <- c("TNF_Score1","Chemokine_Score1","ECM_Score1")
for (f in features) {
  df <- FetchData(obj_cd14_sub, vars = c(f, "TripleGroup"))
  cat("\n=== ", f, " ===\n")
  print(pairwise.wilcox.test(df[[f]], df$TripleGroup, p.adjust.method = "BH"))
}


features <- c("TNF_Score1","Chemokine_Score1","ECM_Score1")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (f in features) {
  df <- FetchData(obj_cd14_sub, vars = c(f,"TripleGroup"))
  cat("\n=== ", f, " ===\n")
  print(wilcox.test(df[df$TripleGroup=="TriplePositive", f],
                    df[df$TripleGroup=="TripleNegative", f]))
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
library(effsize)
for (f in features) {
  df <- FetchData(obj_cd14_sub, vars = c(f,"TripleGroup"))
  cat("\nCliff's delta - ", f, "\n")
  print(cliff.delta(df[df$TripleGroup=="TriplePositive", f],
                    df[df$TripleGroup=="TripleNegative", f]))
}



#
Dog.combined.singlet.T <- readRDS("PATH_TO_FILE")
unique(Dog.combined.singlet.T$Subset)
unique(Dog.combined.singlet.T$Subset3)

library(Seurat)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cells_osa_tumor <- WhichCells(
  Dog.combined.singlet.T,
  expression = Subset == "Tumor" & Subset3 == "OSA"
)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Dog.OSA.Tumor <- subset(Dog.combined.singlet.T, cells = cells_osa_tumor)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
table(Dog.OSA.Tumor$Subset, Dog.OSA.Tumor$Subset3)

DefaultAssay(Dog.OSA.Tumor) <- "RNA"

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cd8_cells <- WhichCells(Dog.OSA.Tumor, expression = CD8A > 0)
Dog.OSA.Tumor.CD8 <- subset(Dog.OSA.Tumor, cells = cd8_cells)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
length(cd8_cells)

# PDCD1 cutoff: >0
Dog.OSA.Tumor.CD8$PD1_group <- ifelse(
  FetchData(Dog.OSA.Tumor.CD8, vars = "PDCD1")[,1] > 0,
  "PDCD1_pos",
  "PDCD1_neg"
)

Idents(Dog.OSA.Tumor.CD8) <- "PD1_group"
table(Dog.OSA.Tumor.CD8$PD1_group)


# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
exhaustion_genes   <- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1")
cytotoxic_genes    <- c("GZMB","PRF1","GNLY","NKG7","GZMK")
activation_genes   <- c("IFNG","TNF","IL2RA","CD69","HLA-DRA")
proliferation_genes<- c("MKI67","TOP2A","PCNA","TYMS","UBE2C")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Dog.OSA.Tumor.CD8 <- AddModuleScore(Dog.OSA.Tumor.CD8, list(exhaustion_genes),   name="Exhaustion")
Dog.OSA.Tumor.CD8 <- AddModuleScore(Dog.OSA.Tumor.CD8, list(cytotoxic_genes),    name="Cytotoxic")
Dog.OSA.Tumor.CD8 <- AddModuleScore(Dog.OSA.Tumor.CD8, list(activation_genes),   name="Activation")
Dog.OSA.Tumor.CD8 <- AddModuleScore(Dog.OSA.Tumor.CD8, list(proliferation_genes),name="Prolif")

library(ggpubr)

features <- c("Exhaustion1","Cytotoxic1","Activation1","Prolif1")

for (f in features) {
  p <- VlnPlot(Dog.OSA.Tumor.CD8, features=f, group.by="PD1_group", pt.size=0) +
    stat_compare_means(method="wilcox.test", label="p.format", ref.group="PDCD1_neg") +
    ggtitle(f)
  print(p)
}

VlnPlot(Dog.OSA.Tumor.CD8, features = c("Exhaustion1","Cytotoxic1","Prolif1"), ncol = 4,
        group.by="PD1_group", pt.size=0, cols = c("#ADDAE6","#E63222"))
ggsave("CD8_TP_TN_Modules_score.png", width = 6, height = 3)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
features <- c("Exhaustion1","Cytotoxic1","Activation1","Prolif1")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (f in features) {
  df <- FetchData(Dog.OSA.Tumor.CD8, vars = c(f, "PD1_group"))
  w <- wilcox.test(df[[f]] ~ df$PD1_group)
  cat(f, ": p =", signif(w$p.value, 3), "\n")
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
FeaturePlot(
  Dog.OSA.Tumor.CD8,
  cols = c("#ADDAE6","#E63222"),
  features = c("PDCD1","CD8A"),
  reduction = "umap",
  pt.size = 1,
  order = TRUE
)
ggsave("UMAP_TP_TN_CD8A_CD8A_OSA.png", width = 7, height = 4)



library(Seurat)
DefaultAssay(Dog.OSA.Tumor.CD8) <- "RNA"

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
pdcd1_cd8_cells <- WhichCells(
  Dog.OSA.Tumor.CD8,
  expression = CD8A > 0 & PDCD1 > 0
)

Dog.OSA.Tumor.CD8.PD1pos <- subset(Dog.OSA.Tumor.CD8, cells = pdcd1_cd8_cells)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ncol(Dog.OSA.Tumor.CD8.PD1pos)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Dog.OSA.Tumor.CD8.PD1pos$LAG3_group <- ifelse(
  FetchData(Dog.OSA.Tumor.CD8.PD1pos, vars = "LAG3")[,1] > 0,
  "LAG3_pos",
  "LAG3_neg"
)

Idents(Dog.OSA.Tumor.CD8.PD1pos) <- "LAG3_group"
table(Dog.OSA.Tumor.CD8.PD1pos$LAG3_group)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
DimPlot(
  Dog.OSA.Tumor.CD8.PD1pos,
  reduction = "umap",
  group.by = "LAG3_group",
  cols = c("#ADDAE6","#E63222"),
  pt.size = 1
)  + NoLegend()
ggsave("UMAP_TP_TN_CD8A_PDCD1_LAG3_OSA.png", width = 4, height = 4)
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
exhaustion_genes <- c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1")
Dog.OSA.Tumor.CD8.PD1pos <- AddModuleScore(Dog.OSA.Tumor.CD8.PD1pos, list(exhaustion_genes), name="Exhaustion")
Dog.OSA.Tumor.CD8.PD1pos <- AddModuleScore(
  Dog.OSA.Tumor.CD8.PD1pos,
  list(response_to_virus_genes, exhaustion_genes),
  name = c("ResponseToVirus","Exhaustion1")
)
VlnPlot(
  Dog.OSA.Tumor.CD8.PD1pos,
  features = c("Exhaustion1","ResponseToVirus"),
  group.by = "LAG3_group",
  cols = c("#ADDAE6","#E63222"),
  pt.size = 0.1
)


response_to_virus_genes <- c(
  "ADAR","ARF1","ATF2","BAX","CALR","CCL19","CCL5","CHUK","CXCL10","DDX3X",
  "EIF5A","EXT1","FMR1","GBF1","GSDME","HIF1A","HSP90AA1","IFI6","IFIH1",
  "IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21",
  "IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNE",
  "IFNG","IFNGR1","IFNGR2","IFNK","IFNL1","IFNL2","IFNL3","IFNL4","IFNLR1",
  "IFNW1","IKBKE","IL10RB","IL12A","IL21","IL6","IRF3","IRF5","IRGM","JAK1",
  "JAK2","LGALS8","LGALS9","MAPK11","MAPK14","MIR130A","MIR146A","MIR21",
  "MIR29B1","MIR302A","MIR30C1","MIR675","MIR758","MMP12","NFKB1","NLRP3",
  "OAS1","PENK","POU2AF1","POU2F2","PYDC5","RIOK3","RRP1B","SMAD3","TGFB1",
  "TLR3","TLR7","TOMM70","TRIM6","TYK2","VWCE","WDFY4","ZC3H12A"
)

# Exhaustion-related genes
exhaustion_genes <- c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","TOX","ENTPD1")

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
response_to_virus_genes <- c(
  "ADAR","ARF1","ATF2","BAX","CALR","CCL19","CCL5","CHUK","CXCL10","DDX3X",
  "EIF5A","EXT1","FMR1","GBF1","GSDME","HIF1A","HSP90AA1","IFI6","IFIH1",
  "IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21",
  "IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNE",
  "IFNG","IFNGR1","IFNGR2","IFNK","IFNL1","IFNL2","IFNL3","IFNL4","IFNLR1",
  "IFNW1","IKBKE","IL10RB","IL12A","IL21","IL6","IRF3","IRF5","IRGM","JAK1",
  "JAK2","LGALS8","LGALS9","MAPK11","MAPK14","MMP12","NFKB1","NLRP3","OAS1",
  "SMAD3","TGFB1","TLR3","TLR7","TOMM70","TRIM6","TYK2","WDFY4","ZC3H12A"
)

Dog.OSA.Tumor.CD8.PD1pos <- AddModuleScore(
  Dog.OSA.Tumor.CD8.PD1pos,
  list(exhaustion_genes),
  name = "Exhaustion"
)

Dog.OSA.Tumor.CD8.PD1pos <- AddModuleScore(
  Dog.OSA.Tumor.CD8.PD1pos,
  list(response_to_virus_genes),
  name = "ResponseToVirus"
)

library(ggpubr)

features <- c("Exhaustion1","ResponseToVirus1")

for (f in features) {
  df <- FetchData(Dog.OSA.Tumor.CD8.PD1pos, vars = c(f, "LAG3_group"))
  w <- wilcox.test(df[[f]] ~ df$LAG3_group)
  cat(f, ": p =", signif(w$p.value, 3), "\n")
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
VlnPlot(
  Dog.OSA.Tumor.CD8.PD1pos,
  features = features,
  group.by = "LAG3_group",
  pt.size = 0.1,
  cols = c("#ADDAE6","#E63222"),
  ncol = 4
) 
ggsave("PD1pos_CD8_LAG3_Modules.png", width = 6, height = 3)

####
cellchat <- createCellChat(object = Dog.combined.singlet, group.by = "celltype") 
Dog.combined.singlet$celltype <- Idents(Dog.combined.singlet)


##########################################################################################

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
library(Seurat); library(SingleCellExperiment); library(slingshot)
library(ggplot2); library(dplyr); library(mgcv)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Bplasma_clust <- c(0,1,2,3,4,5,6,8,9)
Bobj <- subset(Dog.combined.singlet.B, idents = Bplasma_clust)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
common.cells <- intersect(
  rownames(Dog.combined.singlet.B@reductions$umap@cell.embeddings),
  colnames(Bobj)
)
umap_coords <- Dog.combined.singlet.B@reductions$umap@cell.embeddings[common.cells, ]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

Bobj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",
  assay = DefaultAssay(Bobj)
)
DimPlot(Bobj)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Bobj$B_subtype <- "Naive"
Bobj$B_subtype[Idents(Bobj) %in% c(9)] <- "IFN_B"
Bobj$B_subtype[Idents(Bobj) %in% c(4,5,8)] <- "Plasma"
Bobj$B_subtype <- factor(Bobj$B_subtype, levels = c("Naive","IFN_B","Plasma"))

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
if(is.null(Bobj@reductions$umap)){
  Bobj <- FindVariableFeatures(Bobj)
  Bobj <- ScaleData(Bobj, verbose = FALSE)
  Bobj <- RunPCA(Bobj, npcs = 30, verbose = FALSE)
  Bobj <- RunUMAP(Bobj, dims = 1:30, verbose = FALSE)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
sce <- as.SingleCellExperiment(Bobj)
colLabels(sce) <- Bobj$B_subtype
rd <- as.matrix(Bobj@reductions$umap@cell.embeddings)
reducedDims(sce)$UMAP <- rd

set.seed(1)
sce <- slingshot(sce,
                 clusterLabels = 'B_subtype',
                 reducedDim    = 'UMAP',
                 start.clus    = 'Naive',
                 end.clus      = c('Plasma','IFN_B'))

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
pt <- slingPseudotime(sce, na = FALSE)[,1]
Bobj$pseudotime <- as.numeric(pt)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
pt_scaled <- scales::rescale(pt, to = c(0, 1))
Bobj$pseudotime <- pt_scaled

pt <- slingPseudotime(sce, na = FALSE)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
Bobj$pseudotime <- scales::rescale(rowMeans(pt, na.rm = TRUE), to = c(0, 1))

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
p1 <-FeaturePlot(Bobj, features = "pseudotime", reduction = "umap", label = T, repel = T,
                 cols = c("#ADDAE6","#E63222")
) +
  ggtitle("B-cell trajectory pseudotime (Naive → Plasma)")
p1
ggsave("B-cell trajectory pseudotime (Naive → Plasma).png", width = 4, height = 4)
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
mat <- t(as.matrix(GetAssayData(Bobj, slot = "data")[genes_to_show, , drop = FALSE]))
df  <- as.data.frame(mat); df$pseudotime <- Bobj$pseudotime
long <- tidyr::pivot_longer(df, cols = all_of(genes_to_show),
                            names_to = "gene", values_to = "expr")
fitpred <- long %>%
  group_by(gene) %>%
  do({
    m <- gam(expr ~ s(pseudotime, k = 5), data = .)
    data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 200),
               fit = predict(m, newdata = data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 200))))
  })

ggplot(long, aes(pseudotime, expr)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_line(data = fitpred, aes(pseudotime, fit), size = 1.0) +
  facet_wrap(~gene, scales = "free_y") +
  labs(title = "Key B-cell differentiation genes along pseudotime",
       x = "Pseudotime (Naive → Plasma)", y = "Normalized expression")

library(mgcv)
library(dplyr)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
gam_stats <- long %>%
  group_by(gene) %>%
  do({
    dd <- .
    # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
    dd <- dd[is.finite(dd$pseudotime) & is.finite(dd$expr), , drop = FALSE]
    m  <- gam(expr ~ s(pseudotime, k = 5), data = dd, method = "REML")
    st <- summary(m)$s.table[1, , drop = FALSE]
    data.frame(
      edf      = unname(st[,"edf"]),
      F        = unname(st[,"F"]),
      p_smooth = unname(st[,"p-value"])
    )
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_smooth, method = "BH")) %>%
  arrange(FDR)

gam_stats
#>   gene    edf      F   p_smooth      FDR
#>  IRF4   ...    ...      ...         ...
#>  PRDM1  ...    ...      ...         ...
#>  XBP1   ...    ...      ...         ...
lab_df <- gam_stats %>%
  mutate(label = paste0(gene, " (FDR=", signif(FDR, 2), ")")) %>%
  select(gene, label)

long2 <- long %>% left_join(lab_df, by = "gene")

library(ggplot2)
p <- ggplot(long2, aes(x = pseudotime, y = expr, color = gene)) +
  geom_point(alpha = 0.25, size = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, size = 1.1) +
  facet_wrap(~label, scales = "free_y") +
  scale_color_manual(values = c(IRF4="firebrick", PRDM1="darkorange", XBP1="steelblue")) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime (Naive → Plasma)", y = "Normalized expression",
       title = "Plasma differentiation genes along pseudotime")

p



# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
stat_tbl <- Bobj@meta.data %>%
  group_by(B_subtype) %>%
  summarise(n = n(), median_pt = median(pseudotime), mean_pt = mean(pseudotime))

print(stat_tbl)
p1; p2

slingLineages(sce)

pdf("Bcell_pseudotime_boxplot.pdf", width = 2.8, height = 3)
boxplot(pseudotime ~ B_subtype, data = Bobj@meta.data,
        col = c("#2B83BA", "#F7E967", "#D7191C"),
        ylab = "Pseudotime", xlab = "B subtype",
        main = "Pseudotime distribution across B-cell subtypes")
dev.off()

library(RColorBrewer)
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

plot(reducedDims(sce)$UMAP, col = 'grey80', pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 3, col = colors)


# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
library(slingshot)
library(RColorBrewer)

# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
umap <- reducedDims(sce)$UMAP
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]

# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
sds <- SlingshotDataSet(sce)

# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
plot(
  umap[,1], umap[,2],
  col = "grey80", pch = 16, asp = 1,
  xlab = "UMAP_1", ylab = "UMAP_2",
  main = "B-cell differentiation trajectories"
)
lines(sds, lwd = 3, col = c("red3", "dodgerblue3"))  # Plasma = red, IFN_B = blue

# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
umap_df <- as.data.frame(umap)
umap_df$B_subtype <- colLabels(sce)
centers <- aggregate(umap_df[,1:2], by = list(umap_df$B_subtype), FUN = mean)
colnames(centers) <- c("B_subtype", "x", "y")

# -------------------------------
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
# -------------------------------
# Naive
naive <- subset(centers, B_subtype == "Naive")
points(naive$x, naive$y, pch = 19, col = "navy", cex = 2)
text(naive$x, naive$y, "Naive", pos = 4, col = "navy", font = 2, cex = 1.3)

# IFN-B
ifn <- subset(centers, B_subtype == "IFN_B")
points(ifn$x, ifn$y, pch = 17, col = "dodgerblue3", cex = 2)
text(ifn$x, ifn$y, "IFN-B", pos = 2, col = "dodgerblue3", font = 2, cex = 1.3)

# Plasma
plasma <- subset(centers, B_subtype == "Plasma")
points(plasma$x, plasma$y, pch = 17, col = "red3", cex = 2)
text(plasma$x, plasma$y, "Plasma", pos = 4, col = "red3", font = 2, cex = 1.3)



library(slingshot)
library(grid)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
umap <- reducedDims(sce)$UMAP
sds <- SlingshotDataSet(sce)
dev.off()
png("Bcell_Trajectory_with_arrows.png", width = 1000, height = 1500, res = 300)
# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plot(
  umap[,1], umap[,2],
  col = "grey85", pch = 16, asp = 1,
  xlab = "UMAP_1", ylab = "UMAP_2",
  main = "B-cell differentiation trajectories"
)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
line_cols <- c("#D7191C", "#F7E967")


# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
for (i in seq_along(sds@curves)) {
  crv <- sds@curves[[i]]
  lines(crv$s[,1], crv$s[,2], lwd = 3, col = line_cols[i])
  
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  end_x <- crv$s[nrow(crv$s), 1]
  end_y <- crv$s[nrow(crv$s), 2]
  prev_x <- crv$s[nrow(crv$s)-1, 1]
  prev_y <- crv$s[nrow(crv$s)-1, 2]
  
  arrows(prev_x, prev_y, end_x, end_y,
         col = line_cols[i], lwd = 2, length = 0.15, angle = 25)
}

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
naive_cells <- which(colLabels(sce) == "Naive")
points(mean(umap[naive_cells,1]),
       mean(umap[naive_cells,2]),
       pch = 19, col = "#2B83BA", cex = 1)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
ifn_cells <- which(colLabels(sce) == "IFN_B")
points(mean(umap[ifn_cells,1]),
       mean(umap[ifn_cells,2]),
       pch = 19, col = "#F7E967", cex = 1)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
plasma_cells <- which(colLabels(sce) == "Plasma")
points(mean(umap[plasma_cells,1]),
       mean(umap[plasma_cells,2]),
       pch = 19, col = "#D7191C", cex = 1)
dev.off()



# Plasma differentiation genes
plasma_genes <- c("IRF4", "PRDM1", "XBP1")

# IFN-related B module
ifn_genes <- c("MX1", "ISG15", "IFI6", "IFIT3", "OAS1", "STAT1")

# Add IFN module score
Bobj <- AddModuleScore(Bobj, features = list(ifn_genes), name = "IFN_Module")

# Plot GAM for plasma genes
genes_to_show <- plasma_genes
mat <- t(as.matrix(GetAssayData(Bobj, slot = "data")[genes_to_show, , drop = FALSE]))
df  <- as.data.frame(mat)
df$pseudotime <- Bobj$pseudotime

long <- tidyr::pivot_longer(df, cols = all_of(genes_to_show),
                            names_to = "gene", values_to = "expr")

library(mgcv)
library(dplyr)

fitpred <- long %>%
  group_by(gene) %>%
  do({
    m <- gam(expr ~ s(pseudotime, k = 5), data = .)
    data.frame(
      pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 200),
      fit = predict(m, newdata = data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 200)))
    )
  })

ggplot(fitpred, aes(x = pseudotime, y = fit, color = gene)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("firebrick", "darkorange", "steelblue")) +
  theme_classic(base_size = 14) +
  ylab("Normalized expression") +
  ggtitle("Plasma differentiation genes along pseudotime")
ggsave("Plasma differentiation genes along pseudotime.png", width = 4, height = 3)


library(ggplot2)

ggplot(Bobj@meta.data, aes(x = pseudotime, y = IFN_Module1)) +
  geom_point(size = 0.7, alpha = 0.4, color = "grey60") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, color = "#0072B2", fill = "#BBD7EC", size = 1.2) +
  theme_classic(base_size = 14) +
  xlab("Pseudotime (Naive → IFN-related B)") +
  ylab("IFN module score") +
  ggtitle("Dynamic activation of IFN-related transcriptional program") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


library(ggplot2)

# [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
cols <- c("Naive" = "#0072B2", "IFN_B" = "#F7E967", "Plasma" = "#D7191C")

ggplot(Bobj@meta.data, aes(x = pseudotime, y = IFN_Module1)) +
  geom_point(
    data = subset(Bobj@meta.data, B_subtype != "IFN_B"),
    aes(color = B_subtype),
    size = 1.8, alpha = 0.4
  ) +
  # [KOREAN REMOVED - NEED ENGLISH TRANSLATION]
  geom_point(
    data = subset(Bobj@meta.data, B_subtype == "IFN_B"),
    aes(color = B_subtype),
    size = 3.5, alpha = 0.9, shape = 21, fill = "#F7E967", color = "black", stroke = 0.4
  ) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "black", fill = "#BBD7EC", se = TRUE, size = 1.2) +
  scale_color_manual(values = cols) +
  theme_classic(base_size = 14) +
  ylab("IFN module score") +
  xlab("Pseudotime (Naive → IFN-related B)") +
  ggtitle("Selective activation of IFN-related module in IFN_B subset") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave("Selective_activation_IFN_B_emphasized.png", width = 4, height = 4, dpi = 300)
