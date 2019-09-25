# Author: Subinoy Biswas

## Load Seurat library
library(Seurat)
## Load other libraries needed fr Seurat
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(RColorBrewer)

TILs_PBMC_HD_pbl <- MergeSeurat(object1 = TILs_PBMC,
                            object2 = HD_pbl,
                            add.cell.id1 = "TILs_PBMC",
                            add.cell.id2 = "HD_pbl",
                            project = "TIL_PBMC_HD_pbl")


# Checking identity
head(x = TILs_PBMC_HD_pbl@cell.names)
tail(x = TILs_PBMC_HD_pbl@cell.names)


mito.genes <- grep(pattern = "^MT-", x=rownames(x=TILs_PBMC_HD_pbl@raw.data), value = TRUE)

percent.mito <- Matrix::colSums(TILs_PBMC_HD_pbl@raw.data[mito.genes, ])/Matrix::colSums(TILs_PBMC_HD_pbl@raw.data)

TILs_PBMC_HD_pbl <- AddMetaData(object = TILs_PBMC_HD_pbl, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = TILs_PBMC_HD_pbl,
        features.plot = c("nGene", "nUMI", "percent.mito"),cols.use = NULL,nCol = 3)


TILs_PBMC_HD_pbl <- NormalizeData(object = TILs_PBMC_HD_pbl,
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4)

TILs_PBMC_HD_pbl <- FindVariableGenes(object = TILs_PBMC_HD_pbl,
                                  mean.function = ExpMean,
                                  dispersion.function = LogVMR,
                                  x.low.cutoff = 0.0125,
                                  x.high.cutoff = 3,
                                  y.cutoff = 0.8)

length(TILs_PBMC_HD_pbl@var.genes)

TILs_PBMC_HD_pbl <- ScaleData(object = TILs_PBMC_HD_pbl,
                          vars.to.regress = c("nUMI", "percent.mito"))


TILs_PBMC_HD_pbl <- RunPCA(object = TILs_PBMC_HD_pbl,
                       pc.genes = TILs_PBMC_HD_pbl@var.genes,
                       do.print = TRUE,
                       pcs.print = 1:5,
                       genes.print = 5)

PCAPlot(object = TILs_PBMC_HD_pbl, dim.1 = 1, dim.2 = 2)



TILs_PBMC_HD_pbl <- RunTSNE(object = TILs_PBMC_HD_pbl,
                        dims.use = 1:10,
                        do.fast = TRUE)

TSNEPlot(object = TILs_PBMC_HD_pbl, do.label = TRUE)

TILs_PBMC_HD_pbl <- FindClusters(object = TILs_PBMC_HD_pbl, reduction.type = "pca", dims.use = 1:10,
                             resolution = c(0.3, 0.5, 0.7, 1), print.output = 0, save.SNN = TRUE)

TILs_PBMC_HD_pbl <- SetAllIdent(TILs_PBMC_HD_pbl, id="res.0.5")
levels(TILs_PBMC_HD_pbl@ident)

TSNEPlot(object = TILs_PBMC_HD_pbl, do.label = TRUE)


# First lets stash our identities for later
TILs_PBMC_HD_pbl <- StashIdent(object = TILs_PBMC_HD_pbl, save.name = "ClusterNames_0.5")


TILs_PBMC_HD_pbl <- FindClusters(object = TILs_PBMC_HD_pbl, reduction.type = "pca", dims.use = 1:10, print.output = FALSE)


# Creating Feature plot with Known Markers
FeaturePlot(object = TILs_PBMC_HD_pbl, features.plot = c("IL7R","CD4", "CD3D","CD8A","CD14", "GNLY","LYZ", "HLA-DRB1",
                                                     "MS4A7", "MS4A1", "FCGR3A","MS4A7", "GZMA", "GZMB",
                                                     "GNLY","NKG7","FCER1A", "CST3", "CD1A","CD1C","IL3RA"),
            cols.use = c("grey", "blue")
            ,nCol = 4,reduction.use = "tsne")


# Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria

my_TILs_PBMC <- grep("TILs_PBMC", row.names(TILs_PBMC_HD_pbl@meta.data), value = TRUE)
my_HD_pbl <- grep("HD_pbl", row.names(TILs_PBMC_HD_pbl@meta.data), value = TRUE)
my_TILs_PBMC_tsne <- TSNEPlot(object = TILs_PBMC_HD_pbl, do.label = TRUE, cells.use = my_TILs_PBMC, do.return = TRUE)
my_HD_pbl_tsne <- TSNEPlot(object = TILs_PBMC_HD_pbl, do.label = TRUE, cells.use = my_HD_pbl, do.return = TRUE)

plot_grid(
  my_TILs_PBMC_tsne,
  my_HD_pbl_tsne, labels=c("TILs_PBMC", "HD_pbl")
)



# find markers for every cluster compared to all remaining cells, report
# only the positive ones
TILs_PBMC_HD_pbl.markers <- FindAllMarkers(object = TILs_PBMC_HD_pbl, only.pos = TRUE, min.pct = 0.25,
                                       thresh.use = 0.25)

write.csv(TILs_PBMC_HD_pbl.markers, file="TILs_PBMC_HD_pbl.markers.csv")

# Top 20 markers
TILs_PBMC_HD_pbl.markers_20 <-TILs_PBMC_HD_pbl.markers  %>%
                              group_by(cluster) %>%
                              top_n(20, avg_logFC)

write.csv(TILs_PBMC_HD_pbl.markers_20, file="TILs_PBMC_HD_pbl_cluster_genes_20.csv")

#Visualize a heatmap of markers

TILs_PBMC_HD_pbl.markers_top10 <- TILs_PBMC_HD_pbl.markers %>%
                                  group_by(cluster) %>%
                                  top_n(10, avg_logFC)

## Heatmap
DoHeatmap(object = TILs_PBMC_HD_pbl, genes.use = TILs_PBMC_HD_pbl.markers_top10$gene, slim.col.label = TRUE, remove.key = TRUE)

levels(TILs_PBMC_HD_pbl@ident)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13, 14, 15, 16)

# Assigning new Cluster name
new.cluster.ids <- c("CD4 T cells", "1", "2","NK cells", "B cells", "CD8 T Cells",
                     "6", "7","8 ", "9","CD14+ Monocytes","11", "Dendritic Cells", "13", "FCGR3A+ Monocytes", "15", "16" )

TILs_PBMC_HD_pbl@ident <- plyr::mapvalues(x = TILs_PBMC_HD_pbl@ident, from = current.cluster.ids, to = new.cluster.ids)

# Redrawing tSNE plot
TSNEPlot(object = TILs_PBMC_HD_pbl, do.label = TRUE, pt.size = 0.7)

# Saving seurat object for future use
saveRDS(TILs_PBMC_HD_pbl, file= "TILs_PBMC_HD_pbl.rds")
