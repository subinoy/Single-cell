#load("lung_tsne_cluster.RData")

## Load Seurat library
library(Seurat)
## Load other libraries needed fr Seurat
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(tibble)
library(topGO)
library(org.Hs.eg.db)
library(GO.db)

#PBLs_TILs_HD_filt
cluster_2

## manual check; I already know all cells have >200 genes
table(cluster_2@meta.data$percent.mito < 0.05 & cluster_2@meta.data$nGene<500)

# ****************************************
# Showing number of cells in each cluster
table(cluster_2@meta.data$res.0.7)
# OR
table(cluster_2@ident)

# VIEWING SLOTS
str(cluster_2)

table(cluster_2@meta.data$res.0.7,cluster_2@meta.data$sampleid)

# Store cluster identities in object@meta.data$my.clusters
cluster_2_lung <- Seurat::StashIdent(object = cluster_2,
   save.name = "cluster_2_lung.clusters")

counts <- as.data.frame(table(cluster_2_lung@meta.data$cluster_2_lung.clusters,
   cluster_2_lung@ident))

head(counts)

counts


##  Adding metadata sampleid *************************************
cluster_2@meta.data$sampleid <- cluster_2@cell.names

meta_df <- as.data.frame(cluster_2@meta.data)

meta_df <- meta_df%>% select(-orig.ident)

head(meta_df)
dim(meta_df)

## Replacing long sampleid with the sample id name ******************
#gsub(x = names(a), pattern = "\\.", replacement = "#")

meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "PnTs_PBLS_P607.*", replacement ="PnTs_PBLS_P607" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "PnTs_TILS_T607.*", replacement ="PnTs_TILS_T607" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "PnTs_TILS_T213.*", replacement ="PnTs_TILS_T213" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "PnTs_PBLS_P213.*", replacement ="PnTs_PBLS_P213" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "HDps_HD4_PBL.*", replacement ="HDps_HD4_PBL" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "HDps_HD3_PBL.*", replacement ="HDps_HD3_PBL" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "HDps_HD2_PBL.*", replacement ="HDps_HD2_PBL" )
meta_df$sampleid <- gsub(x=meta_df$sampleid, pattern= "HDps_HD1_PBL.*", replacement ="HDps_HD1_PBL" )
head(meta_df)
tail(meta_df)

## Extracting cluster specific cell distribution from each sample *********
## ************************************************************************
sample_cell_dist_cluster_2_tSNE <- table(meta_df$res.0.7, meta_df$sampleid)

write.csv(sample_cell_dist_cluster_2_tSNE, file= "cluster_2_sample_cell_dist_tSNE.csv")

head(sample_cell_dist_cluster_2_tSNE)

# Takes all the unique cell type specific genes
cluster_2_cell_type_genes <- unique(cluster_2.markers)

cluster_2_cell_type_genes

GOterms.cluster_2 = topGO::topGOterms(fg.genes = cluster_2_cell_type_genes,
   bg.genes = rownames(cluster_2@data), organism = "Human")


# Add gene full name
cluster_2_cell_type_genes$description <- mapIds(x = org.Hs.eg.db,
                                             keys = row.names(cluster_2_cell_type_genes),
                                             column = "GENENAME",
                                             keytype = "SYMBOL",
                                             multiVals = "first")

head(cluster_2_cell_type_genes)

# Add gene symbol
#cluster_2_cell_type_genes$symbol <- row.names(cluster_2_cell_type_genes)

# Add ENTREZ ID
cluster_2_cell_type_genes$entrez <- mapIds(x = org.Hs.eg.db,
                                        keys = row.names(cluster_2_cell_type_genes),
                                        column = "ENTREZID",
                                        keytype = "SYMBOL",
                                        multiVals = "first")

head(cluster_2_cell_type_genes)

# Add ENSEMBL
cluster_2_cell_type_genes$ensembl <- mapIds(x = org.Hs.eg.db,
                                         keys = row.names(cluster_2_cell_type_genes),
                                         column = "ENSEMBL",
                                         keytype = "SYMBOL",
                                         multiVals = "first")


head(cluster_2_cell_type_genes)

dim(cluster_2_cell_type_genes)

write.csv(cluster_2_cell_type_genes, file="cluster_2_lung_gene_annotation.csv")

#____END ____
