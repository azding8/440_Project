### title: "Endometrioma vs unaffected Ovary scRNA-seq Analysis"
### author: "Andrew Ding"
### output: Figures from written report

# Need this version of Seurat to use MAST
remotes::install_version("Seurat", version = "4.1.1")

#Load Libraries
library(Seurat)
library(splitstackshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MAST)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(1)

# Load Data
allcells<-readRDS(file = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Data/Labeled_data/epithelial.annotated.shared.rds")
#table(allcells@meta.data$Major.Class)
Idents(allcells)<-"Major.Class"
# Subset Data
Ovary <- subset(allcells, idents = c("Endometrioma","Unaffected ovary"))
#table(Ovary@meta.data$Major.Class)

#####################
### Figure 1 Left ###
#####################
gg<-DimPlot(Ovary, reduction = "umap", group.by = "active.cluster")+ggtitle("Ovarian Epithelial Cell Subset")
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/umaps/Ovary_EpiAnno_ActiveCluster.png",
       width = 2500, height = 2000, units = "px")
#table(Ovary@meta.data$"active.cluster)

######################
### Figure 1 Right ###
######################
gg<-DimPlot(Ovary, reduction = "umap", group.by = "Major.Class")+ggtitle("Ovarian Epithelial Cell Subset")
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/umaps/Ovary_EpiAnno_MajorClass.png",
       width = 2500, height = 2000, units = "px")
#table(Ovary@meta.data$Major.Class)

#####################
### Figure 2 Left ###
#####################
Idents(object = Ovary) <- "Major.Class"
#head(Idents(Ovary))
# Find DEGs
all.markers <- FindAllMarkers(object = Ovary, test.use="MAST", return.thresh = 0.05)
write.csv(all.markers, file ="/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/CSVs/Epianno_ovary_major_class_MAST.csv")
all.markers<-read.csv("/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/CSVs/Epianno_ovary_major_class_MAST.csv")

all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  ungroup() -> DEGs
gg<-DoHeatmap(Ovary,
              features = DEGs$gene,
              group.by = "Major.Class")+ 
  theme(axis.text.y=element_text(size=8),legend.text=element_text(size=14))
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/heatmaps/Ovary_EpiAnno_MajorClass_Heatmap_by_MajorClass.png",
       width = 3000, height = 3500, units = "px")

#######################
### Figure 2 Middle ###
#######################
Idents(Ovary)<-"Major.Class"
Endometrioma <- subset(Ovary, idents = "Endometrioma")
Unaffected <- subset(Ovary, idents = "Unaffected ovary")

gg<-DoHeatmap(Endometrioma,
              features = DEGs$gene,
              group.by = "active.cluster")+ 
  theme(axis.text.y=element_text(size=8),legend.text=element_text(size=14))
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/heatmaps/Ovary_EpiAnno_Endometrioma_Heatmap_by_ActiveCluster.png",
       width = 3000, height = 3500, units = "px")

######################
### Figure 1 Right ###
######################
gg<-DoHeatmap(Unaffected,
              features = DEGs$gene,
              group.by = "active.cluster")+ 
  theme(axis.text.y=element_text(size=8),legend.text=element_text(size=14))
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/heatmaps/Ovary_EpiAnno_Unaffected_Heatmap_by_ActiveCluster.png",
       width = 3000, height = 3500, units = "px")

################
### Figure 3 ###
################
all.markers<-read.csv("/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/CSVs/Epianno_ovary_major_class_MAST.csv")

all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  ungroup() -> absDEGs
gene <- absDEGs$gene
# Map to ENTREZ IDs
symbol_to_entrez <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Extract the mapped ENTREZ IDs
gene <- as.character(symbol_to_entrez$ENTREZID)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

# Over-representation analysis
ego <- enrichGO(gene          = gene,
                universe      = names(gene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
# Plot
gg<-goplot(ego)
ggsave(plot = gg,
       filename = "/Volumes/Extreme Pro/Endometrium scRNAseq Analysis/Plots/epithelial.annotated.shared/enrichment/Ovary_EpiAnno_MajorClass_Enrichment.png",
       width = 3500, height = 3000, units = "px")
