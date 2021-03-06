### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyr)
library(Seurat)
library(tidyverse)

PMID<-"31753849"
Species<-"Human"
Tissue<-"Intestine"
paper.title<-"Single-cell transcriptome analysis reveals differential nutrient absorption functions in human intestine"
GSE.id<-"GSE125970"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")
GEO.link<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id)
  
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
file.script.name<-paste0(path.script, "formatting_",GSE.id,".R")

if (!file.exists(file.rmd.name)) {file.create(file.rmd.name)}
if (!file.exists(file.script.name)) {file.create(file.script.name)}


write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("GEO link: ", GEO.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)

# 
getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
system(paste0("gunzip ",path.raw.data,"*.gz"))



# metadata
metadata.raw <- read.delim(file = paste0(path.raw.data,
                                         "GSE125970_cell_info.txt"),
                           stringsAsFactors = F)
unique(metadata.raw$Sample_ID)

metadata.raw$Group<-gsub(metadata.raw$Sample_ID, 
                         pattern = "-\\d", 
                         replacement = "")
unique(metadata.raw$Group)
unique(metadata.raw$CellType)

metadata.raw$UniqueCell_ID<-gsub(metadata.raw$UniqueCell_ID, 
                         pattern = "-| ", 
                         replacement = ".", perl = T)

metadata<-split(metadata.raw, metadata.raw$Group)
metadata[["merged"]]<-metadata.raw


list.files(path=path.raw.data,recursive = TRUE,full.names = T)
SCSeq.mtx.raw<-read.delim(file = paste0(path.raw.data,
                                        "GSE125970_raw_UMIcounts.txt"),
                          stringsAsFactors = F)
head(SCSeq.mtx.raw[,1:2])
names(SCSeq.mtx.raw) %>% head()
rownames(SCSeq.mtx.raw)<-SCSeq.mtx.raw$GENE
setdiff(names(SCSeq.mtx.raw), metadata.raw$UniqueCell_ID)
setdiff(metadata.raw$UniqueCell_ID, names(SCSeq.mtx.raw))

SCSeq.mtx<-list()
i<-1
for (i in 1:length(metadata)) {
  SCSeq.mtx[[i]] <- SCSeq.mtx.raw[ , metadata[[i]]$UniqueCell_ID]
  # SCSeq.mtx[[i]] %>% str() %>% print()
  all.equal(metadata[[i]]$UniqueCell_ID, colnames(SCSeq.mtx[[i]])) %>% print()
}
names(SCSeq.mtx)<-names(metadata)
head(SCSeq.mtx[[1]][,1:5])
dim(SCSeq.mtx[[1]])
dim(SCSeq.mtx[[2]])
dim(SCSeq.mtx[[3]])



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
  cell2clusterAssignment<-metadata[[i]][,c(1,3)]
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
