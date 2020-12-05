### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyr)
library(Seurat)
# library(tidyverse)

PMID<-"32855340"
Species<-"Mosquito"
Tissue<-"Blood"
paper.title<-"Mosquito cellular immunity at single-cell resolution"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")

### create directory
# path.use <- "C:/Users/gaoyu/Downloads/"
path.use <- "/n/groups/flyrnai/Yue/scseq_data_formatting/processing_data/"
foldername<-paste(PMID,Species,Tissue,sep="_")
system("tree")

path.raw.data <- paste0(path.use, foldername, "/raw_data/")
path.result <- paste0(path.use, foldername, "/results/")
path.script <- paste0(path.use, foldername, "/script/")
if (!dir.exists(path.raw.data)) {dir.create(path.raw.data, recursive = T)}
if (!dir.exists(path.result)) {dir.create(path.result, recursive = T)}
if (!dir.exists(path.script)) {dir.create(path.script, recursive = T)}

file.rmd.name<-paste0(path.use, foldername, "/README.md")
file.script.name<-paste0(path.script, "formatting_",PMID,".R")

if (!file.exists(file.rmd.name)) {file.create(file.rmd.name)}
if (!file.exists(file.script.name)) {file.create(file.script.name)}


write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Count matrix is from Sudhir email 12/05/2020, metadata could be found in below link
https://zenodo.org/record/3882128#.X8ugzs1Kg2w"),file=file.rmd.name,append=TRUE)



### count matix
count.mtx.files<-list.files(path = path.raw.data, "rds", full.names = T)
count.mtx.files
count.mtx.list<-list()

i<-1
for (i in 1:length(count.mtx.files)) {
  tmp<-readRDS(count.mtx.files[i])
  tmp<-UpdateSeuratObject(tmp)
  count.mtx.list[[i]]<-tmp@assays$RNA@counts
}
names(count.mtx.list)<-c("Aedes", "Anopheles")
dim(count.mtx.list[[1]])
dim(count.mtx.list[[2]])



### metadata
metadata.files<-list.files(path = path.raw.data, "meta", full.names = T)
metadata.files
metadata.list<-list()

for (i in 1:length(metadata.files)) {
  tmp<-read.delim(metadata.files[i])
  tmp$Barcode<-rownames(tmp)
  metadata.list[[i]]<-tmp
}
names(metadata.list)<-c("Aedes", "Anopheles")

unique(metadata.list$Aedes$name.ident) %>% length() #15
# View(metadata.list[["Aedes"]])
metadata.list[[1]]<-metadata.list[[1]][,10:11]

unique(metadata.list$Anopheles$name.types.ident) %>% length()
# View(metadata.list[["Anopheles"]])
metadata.list[[2]]<-metadata.list[[2]][,c(8,19)]


colnames(count.mtx.list[[1]]) %>% head()
metadata.list[[1]]$Barcode %>% head()

colnames(count.mtx.list[[2]]) %>% head()
metadata.list[[2]]$Barcode %>% head()


SCSeq.mtx<-list()
i<-2
for (i in 1:length(count.mtx.list)) {
  SCSeq.mtx[[i]] <- count.mtx.list[[i]][ , metadata.list[[i]]$Barcode]
  dim(SCSeq.mtx[[i]]) %>% print()
  all.equal(metadata.list[[i]]$Barcode, colnames(SCSeq.mtx[[i]])) %>% print()
}
names(SCSeq.mtx)<-c("Aedes", "Anopheles")

head(SCSeq.mtx[[1]])
head(SCSeq.mtx[[2]])

metadata<-list()
metadata[[1]]<-metadata.list[[1]][,c(2,1)]
metadata[[2]]<-metadata.list[[2]]
names(metadata)<-names(metadata.list)


i<-1
for (i in 1:length(SCSeq.mtx)) {
  # logNormalizeMatrix <- log1p(sweep(SCSeq.mtx[[i]], 2, Matrix::colSums(SCSeq.mtx[[i]]), FUN = "/") * 10000)
  # 
  # mean(SCSeq.mtx[[i]][,1])
  # mean(logNormalizeMatrix[,1])
  # write.csv(logNormalizeMatrix,
  #           file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
  #           row.names = F)
  
  
  ###annotate cluster to cell
  str(metadata[[i]])
  cell2clusterAssignment<-metadata[[i]]
  dim(cell2clusterAssignment)
  names(cell2clusterAssignment)<-c("Barcode", "Cluster")
  write.csv(cell2clusterAssignment,
            file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)
  
  ###pseudoBulkMatrix
  mat<-SCSeq.mtx[[i]]
  head(mat[,1:5])
  # df<- mat %>% as.matrix() %>% as.data.frame() %>% t()
  # head(df[,1:5])
  # 
  # all.equal(cell2clusterAssignment$Barcode, row.names(df))
  # 
  # tail(cell2clusterAssignment$Barcode)
  # tail(row.names(df))
  # 
  # df <- cbind(cell2clusterAssignment[ , c("Cluster"), drop = FALSE], df)
  # df[, 1:3] %>% head()
  # 
  # # sum gene counts by each celltype
  # df_group_by_celltype <- df[,1:1000] %>%
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
  #           file = paste0(path.result,names(SCSeq.mtx)[i],"_", PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))
  # 
  
  
  ### dot plot
  seuratObj <- CreateSeuratObject(counts = mat)
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  all_genes <- rownames(seuratObj)
  all_genes
  seuratObj <- ScaleData(seuratObj, features = all_genes)
  seuratObj$celltype <- as.character(cell2clusterAssignment$Cluster)
  seurat_metadata <- seuratObj@meta.data
  str(seurat_metadata)
  Idents(seuratObj) <- "celltype"
  
  clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
  colnames(clusterMetadataTable) <- c("celltype_id", "count")
  
  
  # export clusterMetadataTable
  write.csv(clusterMetadataTable,
            file = paste0(path.result,names(SCSeq.mtx)[i],"_",PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
            row.names = FALSE)
  
  
  # extract data matrix from DotPlot function
  dot <- DotPlot(object = seuratObj, features = all_genes)
  gene2clusterTable <- dot$data
  colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")
  
  # export gene2clusterTable
  write.csv(gene2clusterTable,
            file = paste0(path.result,names(SCSeq.mtx)[i],"_",  PMID, "_gene2clusterTable_", Species, "_", Tissue,".csv"),
            row.names = FALSE)
}
