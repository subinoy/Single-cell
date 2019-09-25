## Load Seurat library
library(Seurat)
## Load other libraries needed fr Seurat
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(tibble)

# PBLs_TILs_HD_filt
cluster_2 <- SubsetData(object = PBLs_TILs_HD_filt , ident.use = 2, do.clean = TRUE, do.scale = T)


table(cluster_2@ident)

cluster_2_df <- as.data.frame(as.matrix(cluster_2@data))
#write.csv(cluster_2_df, file="lung_cluster_2.csv")
cluster_2_df
head(cluster_2@meta.data)

is.numeric(cluster_2_df)
str(cluster_2_df)
colnames(cluster_2_df)
dim(cluster_2_df)


filt_clust_2_df <- cluster_2_df %>%
  rownames_to_column(., var="genes") %>%
  mutate(gene_sum=rowSums(.)) %>%
  column_to_rownames(., var= "genes") %>%
  filter(gene_sum >0)


write.csv(filt_clust_0_df, file="cluster_2_all_genes_expressed.csv")

dim(filt_clust_2_df)

cluster_2 <- NormalizeData(cluster_2)

cluster_2 <- FindVariableGenes(cluster_2, num.bin=10)

length(cluster_2@var.genes)
#[1] 206

str(cluster_2@var.genes)

slotNames(cluster_2)

# Select variable genes from object@hvg.info
selected.var.genes <- cluster_2@hvg.info[cluster_2@var.genes,]
head(selected.var.genes)

# Order variable gene data by scaled gene dispersion
selected.var.genes <- selected.var.genes[order(selected.var.genes$gene.dispersion.scaled, decreasing = T),]
selected.var.genes
dim(selected.var.genes)

write.csv(selected.var.genes, file="cluster_2_variable_genes.csv")

top_200 <- selected.var.genes %>%
  rownames_to_column(., var= "row_name") %>%
  top_n(500) %>%
  column_to_rownames(., var="row_name")

head(top_200)

#cluster_2.var.genes <- rownames(top_200)


cluster_2.var.genes <- rownames(cluster_2@hvg.info)[1:206]
cluster_2.var.genes

# mitochondria genes conveniently start with MT
cluster_2_mito.genes <- grep(pattern = "^MT-", x = rownames(x = cluster_2@data), value = TRUE)

length(cluster_2_mito.genes)

cluster_2_percent.mito <- Matrix::colSums(cluster_2@raw.data[mito.genes, ]) / Matrix::colSums(cluster_2@raw.data)

str(cluster_2_percent.mito)
summary(cluster_2_percent.mito)

# check out the meta data
head(cluster_2@meta.data)
dim(cluster_2@meta.data)

# add some more meta data
cluster_2 <- AddMetaData(object = cluster_2,
                        metadata = percent.mito,
                        col.name = "cluster_2_percent.mito")

cluster_2=ScaleData(cluster_2,genes.use=cluster_2.var.genes, vars.to.regress = c("nGene","cluster_2_percent.mito"))

cluster_2 <- RunPCA(cluster_2,pc.genes=cluster_2@var.genes)

PCElbowPlot(cluster_2)
PCHeatmap(cluster_2, pc.use=1:10,cells.use=500,do.balanced=T)


PCAPlot(cluster_2, dim.1=1, dim.2=2)
#Run TSNE

cluster_2 <- RunTSNE(cluster_2,reduction.use="pca",dims.use=1:15,do.fast=T)

TSNEPlot(cluster_2)


#Clustering

cluster_2 <- FindClusters(cluster_2,reduction.type="pca",dims.use=1:15,res=c(0.3,0.5, 0.7,1))

# cluster_2 <- FindClusters(cluster_2,reduction.type="pca",dims.use=1:15,
#                           res=c(0.3,0.5, 0.7,1), save.SNN = T, force.recalc = T)

TSNEPlot(cluster_2,do.label=T,group.by="res.0.5")

TSNEPlot(cluster_2,do.label=T,group.by="res.0.7")

TSNEPlot(cluster_2,do.label=T,group.by="res.1")


cluster_2.markers <- FindAllMarkers(object = cluster_2, only.pos = TRUE, min.pct = 0.25,
                                   thresh.use = 0.25)

write.csv(cluster_2.markers, file="cluster_2.markers_lung_genes_April_8.csv")

cluster_2.markers_20 <-cluster_2.markers  %>% group_by(cluster) %>% top_n(20, avg_logFC)
write.csv(cluster_2.markers_20, file="cluster_2.markers_20_lung_genes_April_8.csv")


cluster_2.markers_top10 <- cluster_2.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

# ****************
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
# *****************
DoHeatmap(object = cluster_2, genes.use = cluster_2.markers_top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# Takes all the unique cell type specific genes
cluster_2_cell_type_genes <- unique(cluster_2.markers)
