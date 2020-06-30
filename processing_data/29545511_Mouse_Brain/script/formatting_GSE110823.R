###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"29545511"
Species<-"Mouse"
Tissue<-"Brain"
GSE.id<-"GSE110823"
paper.title<-"Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding"
paper.link<-"https://pubmed.ncbi.nlm.nih.gov/29545511/"
GEO.link<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110823"
  
### create directory
# path.use<-"/Users/gaoyu/Downloads/"
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
# write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)


# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"*.tar -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))


list.files(path=path.raw.data,recursive = TRUE,full.names = T)


###read matrix
library(rmatio)
raw.data<- read.mat(paste0(path.raw.data,"GSM3017261_150000_CNS_nuclei.mat"))
str(raw.data)

SCSeq.colnames.barcode<- raw.data$barcodes %>% as.character()
head(SCSeq.colnames.barcode)

SCSeq.rownames.gene<- raw.data$genes
SCSeq.rownames.gene<-trimws(SCSeq.rownames.gene)
head(SCSeq.rownames.gene)

SCSeq.mtx<- raw.data$DGE
dim(SCSeq.mtx) 
head(SCSeq.mtx[,1:3])
SCSeq.mtx@Dimnames<-list(SCSeq.colnames.barcode, SCSeq.rownames.gene)
dim(SCSeq.mtx)
head(SCSeq.mtx[,1:3])
class(SCSeq.mtx)
### *** this matrix is a transposed matrix ***  ###



###normallize matrix
# logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)
# mean(SCSeq.mtx[,1])
# mean(logNormalizeMatrix[,1])
# write.csv(logNormalizeMatrix,
#           file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
#           row.names = F)



###annotate cluster to cell
# metadata
metadata <- cbind(SCSeq.colnames.barcode, raw.data$cluster_assignment) %>% as.data.frame()
str(metadata)

cell2clusterAssignment<-metadata
names(cell2clusterAssignment)<-c("Barcode","Cluster")
cell2clusterAssignment$Barcode<-cell2clusterAssignment$Barcode %>% as.character()
cell2clusterAssignment$Cluster<-cell2clusterAssignment$Cluster %>% as.character()
cell2clusterAssignment$Cluster<-trimws(cell2clusterAssignment$Cluster)
cell2clusterAssignment$Cluster
write.csv(cell2clusterAssignment,
          file=paste0(path.result,PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
          row.names = F)

# save(SCSeq.mtx, cell2clusterAssignment, file=paste0(path.script,GSE.id,".mtx.RData"))
# rm(list = ls())
# load(paste0(path.script,GSE.id,".mtx.RData"))
ls()


###pseudoBulkMatrix
# #using loop convert dgCMatrix2 matrix
# Rcpp::sourceCpp(code='
# #include <Rcpp.h>
# using namespace Rcpp;
# 
# 
# // [[Rcpp::export]]
# IntegerMatrix asMatrix(NumericVector rp,
#                        NumericVector cp,
#                        NumericVector z,
#                        int nrows,
#                        int ncols){
# 
#   int k = z.size() ;
# 
#   IntegerMatrix  mat(nrows, ncols);
# 
#   for (int i = 0; i < k; i++){
#       mat(rp[i],cp[i]) = z[i];
#   }
# 
#   return mat;
# }
# ' )
# 
# as_matrix <- function(mat){
#   
#   row_pos <- mat@i
#   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
#   
#   tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
#                   nrows =  mat@Dim[1], ncols = mat@Dim[2])
#   
#   row.names(tmp) <- mat@Dimnames[[1]]
#   colnames(tmp) <- mat@Dimnames[[2]]
#   return(tmp)
# }
# 
# df<-as_matrix(SCSeq.mtx)
# df <- df  %>% as.data.frame()
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
mat<-SCSeq.mtx %>% t()
dim(mat)
rm(SCSeq.mtx)
seuratObj <- CreateSeuratObject(counts = mat)
rm(mat)
ls()
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", 
                                  nfeatures = 2000)
all_genes <- rownames(seuratObj)
seuratObj <- ScaleData(seuratObj, features = all_genes)
seuratObj$celltype <- as.character(cell2clusterAssignment$Cluster)
seurat_metadata <- seuratObj@meta.data
str(seurat_metadata)
Idents(seuratObj) <- "celltype"
# saveRDS(seuratObj, file = file=paste0(path.script,GSE.id,".seuratObj.Rds"))

clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
colnames(clusterMetadataTable) <- c("celltype_id", "count")
# export clusterMetadataTable
write.csv(clusterMetadataTable,
          file = paste0(path.result, PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)



# extract data matrix from DotPlot function
# dot <- DotPlot(object = seuratObj, features = all_genes)
# ### Warning: Could not find ... in the default search locations, found in RNA assay instead

# i<-1
gene2clusterTable<-NULL
for (i in 1:length(all_genes)) {
  tmp<-DotPlot(object = seuratObj, features = all_genes[i]) 
  gene2clusterTable <- rbind(gene2clusterTable,tmp$data)
}

colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")
gene2clusterTable<-arrange(gene2clusterTable,desc(celltype_id))

# export gene2clusterTable
write.csv(gene2clusterTable,
          file = paste0(path.result, PMID, "_gene2clusterTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)


