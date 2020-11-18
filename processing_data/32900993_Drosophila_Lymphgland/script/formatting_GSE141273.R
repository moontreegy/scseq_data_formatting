### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)

PMID<-"32900993"
Species<-"Drosophila"
Tissue<-"Lymphgland"
paper.title<-"Single-cell transcriptome maps of myeloid blood cell lineages in Drosophila"
GSE.id<-"GSE141273"
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


# write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("GEO link: ", GEO.link,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Barcode2clutser info is obtained from: https://github.com/sangho1130/Dmel_Dropseq/tree/master/3.Revision/updated_figures/Normal_LG/rdata
# \n"),file=file.rmd.name,append=TRUE)
# write(paste0("The cell numbers are updated in their provided barcode2cluster info, so it is different from paper, but the cell types are consistant\n
#              There are 2227 cells (72h), 9383 cells (96h), 7829 cells (120h), 19439 cells in total"),file=file.rmd.name,append=TRUE)

# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))
# 


list.files(path=path.raw.data,recursive = TRUE,full.names = T)

# metadata
metadata.raw <- readRDS(file = paste0(path.raw.data,"/label.Rds"))
# tmp<-subset(metadata.raw, orig.ident=="Drop-seq 72hr tmp")
metadata.raw$Barcode<-gsub(rownames(metadata.raw), 
                           pattern = "-", replacement = "\\.")
rownames(metadata.raw) <- NULL
unique(metadata.raw$anno_simple)
unique(metadata.raw$orig.ident)

names(metadata.raw)[7]<-c("Cluster")
metadata.raw$Cluster<-as.character(metadata.raw$Cluster)

group.names<-c("AEL72hr", "AEL96hr", "AEL120hr", "merged")

metadata<-list()
i<-1
for (i in 1:3) {
  tmp<-subset(metadata.raw, timepoint==group.names[i])
  metadata[[i]]<-tmp[, c("Barcode", "Cluster")]
}
metadata[[4]]<-metadata.raw[, c("Barcode", "Cluster")]
names(metadata)<-group.names




use.file.path<-list.files(path=path.raw.data,
                          pattern = "Dropseq.LymphGland.Normal",
                          recursive = TRUE,
                          full.names = T)
use.file.path

file.name<-list.files(path=path.raw.data,
                      pattern = "Dropseq.LymphGland.Normal",
                      recursive = TRUE,
                      full.names = F)
file.name
file.name<-gsub(file.name, pattern = "\\.lib", replacement = "_lib")
file.name <-
  lapply(file.name, 
                    FUN = function(x){strsplit(x, split = "\\.") %>% 
                        unlist() %>% .[4] }) %>% unlist()
file.name



tmp.mtx<-list()
# i<-1
for (i in 1:length(use.file.path)) {
  tmp.mtx[[i]]<-read.delim(use.file.path[i],  stringsAsFactors = F)
}
names(tmp.mtx)<-file.name
str(tmp.mtx)


SCSeq.mtx<-list()
SCSeq.mtx[[1]]<-tmp.mtx[1:5] %>% reduce(left_join, by = "Symbol")
SCSeq.mtx[[2]]<-tmp.mtx[6:10] %>% reduce(left_join, by = "Symbol")
SCSeq.mtx[[3]]<-tmp.mtx[11:14] %>% reduce(left_join, by = "Symbol")
SCSeq.mtx[[4]]<-tmp.mtx %>% reduce(left_join, by = "Symbol")
names(SCSeq.mtx)<-group.names

head(SCSeq.mtx[[1]][,1:5])
head(SCSeq.mtx[[2]][,1:5])
# str(SCSeq.mtx)

i<-1
for (i in 1:length(SCSeq.mtx)) {
  rownames(SCSeq.mtx[[i]])<-SCSeq.mtx[[i]]$Symbol
  SCSeq.mtx[[i]] <- SCSeq.mtx[[i]][ , metadata[[i]]$Barcode]
  # SCSeq.mtx[[i]] %>% str() %>% print()
  all.equal(metadata[[i]]$Barcode, colnames(SCSeq.mtx[[i]])) %>% print()
}
head(SCSeq.mtx[[2]][,1:5])

# save(SCSeq.mtx, metadata ,file = paste0(path.script, "/mtx.RData"))
# load(paste0(path.script, "/mtx.RData"))







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
