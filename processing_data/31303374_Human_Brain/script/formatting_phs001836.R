###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)

PMID<-"31303374"
Species<-"Human"
Tissue<-"Brain"
GSE.id<-"phs001836"
paper.title<-"A Single-Cell Transcriptomic Atlas of Human Neocortical Development during Mid-gestations"
paper.link<-"https://pubmed.ncbi.nlm.nih.gov/31303374/"
GEO.link<-"http://solo.bmap.ucla.edu/shiny/webapp/"
  
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


# write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Data source link: ", GEO.link,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Note: metadata have two configuration, cluster and subcluster, generating results of two version\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)


list.files(path=path.raw.data,recursive = TRUE,full.names = T)


###read matrix
library(Matrix)

# load in sparse matrix
load(paste0(path.raw.data,"sc_dev_cortex_geschwind/raw_counts_mat.rdata"))
SCSeq.mtx <- as.matrix(raw_counts_mat)

dim(SCSeq.mtx)
# [1] 35543 33986

head(SCSeq.mtx[, 1:3])

# metadata
metadata <- read.csv(paste0(path.raw.data, "sc_dev_cortex_geschwind/cell_metadata.csv"), header = TRUE,stringsAsFactors = F)
str(metadata)
dim(metadata)
# [1] 33976    14

setdiff(metadata$Cell, colnames(SCSeq.mtx))
setdiff( colnames(SCSeq.mtx), metadata$Cell)
# [1] "AAAAAAAAAAAA.y.2" "AAAAAAAAAAAA.x.3" "AAAAAAAAAAAA.y.5" "AAAAAAAAAAAA.x.6"
# [5] "AAAAAAAAAAAA.y.6" "AAAAAAAAAAAA.x.7" "AAAAAAAAAAAA.y.7" "AAAAAAAAAAAA.x.8"
# [9] "AAAAAAAAAAAA.y.8" "AAAAAAAAAAAA"



SCSeq.mtx <- SCSeq.mtx[ , metadata$Cell]
dim(SCSeq.mtx)


class(SCSeq.mtx)
head(SCSeq.mtx[,1:3])

# save.image(file=paste0(path.script,GSE.id,".mtx.RData"))
load(paste0(path.script,GSE.id,".mtx.RData"))

# ###normallize matrix
# logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)
# 
# mean(SCSeq.mtx[,1])
# mean(logNormalizeMatrix[,1])
# write.csv(logNormalizeMatrix,
#           file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
#           row.names = F)



# rm(logNormalizeMatrix)

###annotate cluster to cell
cell2clusterAssignment.sub.group<-list()
cell2clusterAssignment.sub.group[[1]]<-metadata[,c("Cell","Cluster")]
cell2clusterAssignment.sub.group[[2]]<-metadata[,c("Cell","Subcluster")]
names(cell2clusterAssignment.sub.group)<-c("Cluster", "Subcluster")

###pseudoBulkMatrix
for (i in 2:length(cell2clusterAssignment.sub.group)) {
  prefix.name<-names(cell2clusterAssignment.sub.group)[i]

  names(cell2clusterAssignment.sub.group[[i]])<-c("Barcode","Cluster")
  
  dim(cell2clusterAssignment.sub.group[[i]])
  cell2clusterAssignment.sub.group[[i]]$Cluster<-replace( x<-cell2clusterAssignment.sub.group[[i]]$Cluster, 
                                                          is.na(x) ,"Mic_no_subclusters" ) 
  
  write.csv(cell2clusterAssignment.sub.group[[i]],
            file=paste0(path.result,prefix.name,"_", PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)  
  
  
  
  mat<-SCSeq.mtx[, cell2clusterAssignment.sub.group[[i]]$Barcode]
  df<- mat %>% as.matrix() %>% as.data.frame() %>% t()
  str(df)
  
  all.equal(cell2clusterAssignment.sub.group[[i]]$Barcode, row.names(df))
  
  tail(cell2clusterAssignment.sub.group[[i]]$Barcode)
  tail(row.names(df))
  
  df <- cbind(cell2clusterAssignment.sub.group[[i]][ , c("Cluster"), drop = FALSE], df)
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
            file = paste0(path.result,prefix.name,"_", PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))
  
  
  
  ### dot plot
  seuratObj <- CreateSeuratObject(counts = mat)
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  all_genes <- rownames(seuratObj)
  seuratObj <- ScaleData(seuratObj, features = all_genes)
  seuratObj$celltype <- as.character(cell2clusterAssignment.sub.group[[i]]$Cluster)
  seurat_metadata <- seuratObj@meta.data
  str(seurat_metadata)
  Idents(seuratObj) <- "celltype"
  
  clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
  colnames(clusterMetadataTable) <- c("celltype_id", "count")
  

    # export clusterMetadataTable
  write.csv(clusterMetadataTable,
            file = paste0(path.result,prefix.name,"_", PMID,"_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
            row.names = FALSE)
  
  
  # extract data matrix from DotPlot function
  dot <- DotPlot(object = seuratObj, features = all_genes)
  gene2clusterTable <- dot$data
  colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")
  
  # export gene2clusterTable
  write.csv(gene2clusterTable,
            file = paste0(path.result,prefix.name,"_", PMID, "_gene2clusterTable_", Species, "_", Tissue,".csv"),
            row.names = FALSE)
}







