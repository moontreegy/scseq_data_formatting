mat<-SCSeq.mtx
# df <- mat %>% as.matrix() %>% as.data.frame() %>% t()
# df[, 1:3] %>% head()
# 
# all.equal(cell2clusterAssignment$Barcode, row.names(df))
# tail(cell2clusterAssignment$Barcode)
# tail(row.names(df))
# 
# df <- cbind(cell2clusterAssignment[ , c("Cluster"), drop = FALSE], df)
# df[, 1:3] %>% head()
# 
# # sum gene counts by each celltype
# df_group_by_celltype <- df %>%
#   group_by(Cluster) %>%
#   summarise_all(sum) %>%
#   as.data.frame()
# 
# class(df_group_by_celltype)
# 
# row.names(df_group_by_celltype) <- df_group_by_celltype$Cluster
# df_group_by_celltype$Cluster <- NULL
# 
# df_group_by_celltype <- df_group_by_celltype %>% t()
# 
# df_group_by_celltype[1:3, 1:3]
# str(df_group_by_celltype)
# 
# # export pseudoBulkMatrix
# pseudoBulkMatrix <- df_group_by_celltype
# write.csv(pseudoBulkMatrix,
#           file = paste0(path.result, PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))



### dot plot
seuratObj <- CreateSeuratObject(counts = mat)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(seuratObj)
seuratObj <- ScaleData(seuratObj, features = all_genes)
seuratObj$celltype <- as.character(cell2clusterAssignment$Cluster)
seurat_metadata <- seuratObj@meta.data
str(seurat_metadata)
Idents(seuratObj) <- "celltype"

clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
colnames(clusterMetadataTable) <- c("celltype_id", "count")

# cluster.anno<-read.delim(paste0(path.result,"annotate_cluster_name.txt"),
#                          header = F)
# names(cluster.anno)<-c("celltype_id","cell_name")
# clusterMetadataTable<-merge(clusterMetadataTable, cluster.anno)

# export clusterMetadataTable
write.csv(clusterMetadataTable,
          file = paste0(path.result, PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)


# extract data matrix from DotPlot function
dot <- DotPlot(object = seuratObj, features = all_genes)
gene2clusterTable <- dot$data
colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")

# export gene2clusterTable
write.csv(gene2clusterTable,
          file = paste0(path.result, PMID, "_gene2clusterTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)
