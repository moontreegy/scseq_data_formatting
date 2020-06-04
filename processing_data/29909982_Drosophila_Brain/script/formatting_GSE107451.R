###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"29909982"
Species<-"Drosophila"
Tissue<-"Brain"
GSE.id<-"GSE107451"
paper.title<-"A Single-Cell Transcriptome Atlas of the Aging Drosophila Brain"
paper.link<-"https://pubmed.ncbi.nlm.nih.gov/29909982/"
GEO.link<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451"
  
### create directory
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


# write(paste0("PMID: ", PMID,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Species: ", Species,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Tissue: ", Tissue,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper title: ", paper.title,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper link: ", paper.link,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("GEO link: ", GEO.link,"\n\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Note: ", "file 'FBP_meta_annotation_formated.tsv' in rawdata is from Sudir, and I use the metadata from supplementary","\n\n"),
#       file=file.rmd.name,append=TRUE)
# write(paste0("Note: ", "I processed 57k data, and they also have raw 157k datasets to do filtering","\n\n"),
#       file=file.rmd.name,append=TRUE)
# write(paste0("Note: ", "the dataset has 8 different age stages, I seprated them into different group to generate results","\n\n"),
#       file=file.rmd.name,append=TRUE)


# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("gunzip ",path.raw.data,"*.gz"))
# system(paste0("tar -xvf  ",path.raw.data,"*.tar -C", path.raw.data))
# 
list.files(path=path.raw.data,recursive = TRUE,full.names = T)


###read matrix
library(Matrix)
SCSeq.mtx<- readMM(paste0(path.raw.data,"matrix.mtx"))

SCSeq.colnames.barcode<- read.delim(paste0(path.raw.data,"barcodes.tsv"),
                                    stringsAsFactors = F,header = F)
head(SCSeq.colnames.barcode)
SCSeq.rownames.gene<- read.delim(paste0(path.raw.data,"genes.tsv"),
                                 stringsAsFactors = F,header = F)
head(SCSeq.rownames.gene)

dim(SCSeq.mtx) #17473 56902
head(SCSeq.mtx[,1:3])
colnames(SCSeq.mtx) <- SCSeq.colnames.barcode$V1
rownames(SCSeq.mtx) <- SCSeq.rownames.gene$V2
head(SCSeq.mtx[, 1:3])

# metadata
metadata <- read.delim(paste0(
  path.raw.data, "GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv"), 
  header = TRUE,stringsAsFactors = F)
str(metadata)
unique(metadata$res.2) %>% sort
dim(metadata) #56902    31
intersect(metadata$new_barcode, SCSeq.colnames.barcode$V1) %>% length()
SCSeq.colnames.barcode$V1 %>% length()


class(SCSeq.mtx)

# save.image(file=paste0(path.script,GSE.id,".mtx.RData"))
# load(paste0(path.script,GSE.id,".mtx.RData"))

###normallize matrix
# logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)
# mean(SCSeq.mtx[,1])
# mean(logNormalizeMatrix[,1])
# write.csv(logNormalizeMatrix,
#           file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
#           row.names = F)


###annotate cluster to cell
str(metadata)
cell2clusterAssignment<-metadata[,c("new_barcode","res.2","Age")]
dim(cell2clusterAssignment)
names(cell2clusterAssignment)<-c("Barcode","Cluster","Age")
# write.csv(cell2clusterAssignment,
#           file=paste0(path.result,PMID,"_cell2clusterAssignment_",Species,"_",Tissue,"_ALL",".csv"),
#           row.names = F)

cell2clusterAssignment.age.group<-split(cell2clusterAssignment,cell2clusterAssignment$Age)
str(cell2clusterAssignment.age.group)
# # $ 0 :'data.frame':     12105 obs. of  3 variables:
# # $ 1 :'data.frame':     9173 obs. of  3 variables:
# # $ 3 :'data.frame':     7962 obs. of  3 variables:



###pseudoBulkMatrix
for (i in 1:length(cell2clusterAssignment.age.group)) {
  
  names(cell2clusterAssignment.age.group[[i]])<-c("Barcode","Cluster","Age")
  write.csv(cell2clusterAssignment.age.group[[i]],
            file=paste0(path.result,"Age_",names(cell2clusterAssignment.age.group)[i],"_", PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)  
  
  
  mat<-SCSeq.mtx[, cell2clusterAssignment.age.group[[i]]$Barcode]
  df<- mat %>% as.matrix() %>% as.data.frame() %>% t()
  str(df)

  all.equal(cell2clusterAssignment.age.group[[i]]$Barcode, row.names(df))
  
  tail(cell2clusterAssignment.age.group[[i]]$Barcode)
  tail(row.names(df))
  
  df <- cbind(cell2clusterAssignment.age.group[[i]][ , c("Cluster"), drop = FALSE], df)
  df[, 1:3] %>% head()
  
  # sum gene counts by each celltype
  df_group_by_celltype <- df %>%
    group_by(Cluster) %>%
    summarise_all(sum) %>%
    as.data.frame()
  
  class(df_group_by_celltype)
  
  row.names(df_group_by_celltype) <- df_group_by_celltype$Cluster
  df_group_by_celltype$Cluster <- NULL
  
  df_group_by_celltype <- df_group_by_celltype %>% t()
  
  df_group_by_celltype[1:3, 1:3]
  str(df_group_by_celltype)
  
  # export pseudoBulkMatrix
  pseudoBulkMatrix <- df_group_by_celltype
  write.csv(pseudoBulkMatrix,
            file = paste0(path.result,"Age_",names(cell2clusterAssignment.age.group)[i], "_", PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))
  
  
  
  ### dot plot
  seuratObj <- CreateSeuratObject(counts = mat)
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  all_genes <- rownames(seuratObj)
  seuratObj <- ScaleData(seuratObj, features = all_genes)
  seuratObj$celltype <- as.character(cell2clusterAssignment.age.group[[i]]$Cluster)
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
            file = paste0(path.result,"Age_",names(cell2clusterAssignment.age.group)[i],  "_",PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
            row.names = FALSE)
  
  
  # extract data matrix from DotPlot function
  dot <- DotPlot(object = seuratObj, features = all_genes)
  gene2clusterTable <- dot$data
  colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")
  
  # export gene2clusterTable
  write.csv(gene2clusterTable,
            file = paste0(path.result,"Age_",names(cell2clusterAssignment.age.group)[i],"_",  PMID, "_gene2clusterTable_", Species, "_", Tissue,".csv"),
            row.names = FALSE)
}




