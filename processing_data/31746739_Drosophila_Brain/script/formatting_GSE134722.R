###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
rm(list=ls())

PMID<-"31746739"
Species<-"Drosophila"
Tissue<-"Brain"
GSE.id<-"GSE134722"
paper.title<-"Single cell transcriptome atlas of the Drosophila larval brain"
paper.link<-"https://pubmed.ncbi.nlm.nih.gov/31746739/"
GEO.link<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134722"
  
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
# write(paste0("GEO link: ", GEO.link,"\n"),file=file.rmd.name,append=TRUE)
# write(paste0("Note:\n
# \tSupplementary providing two conditions in this paper, see raw data folder.\n
# \t\tStarvationCondition_finalaggr: 17493 genes, 4645 cells \n
# \t\tNormalCondition_finalaggr: 17493 genes, 4708 cells\n
# \tSudhir contacted author getting two metadata folders, also contain normalized matrix, only using metadata for processing, see raw data folder.\n
# \t\tCondition based data processing\n
# \t\t\tUsing NormalMetadata for NormalCondition_finalaggr, getting 4407 cells (the number is consistent with paper) \n
# \t\t\tStarvationMetadata only contains information about merged condition, so using StarvationMetadata match StarvationCondition_finalaggr, getting 4347 cells\n
# \t\tMerged condition\n
# \t\t\tUsing StarvationMetadata to match StarvationCondition_finalaggr and NormalCondition_finalaggr, getting 4347 and 4349 cells, 8696 in total (which is consistent with paper)\n\n"),
#       file=file.rmd.name,append=TRUE)
# write(paste0("Result location: ", path.result,"\n\n"),file=file.rmd.name,append=TRUE)



# 
# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# # system(paste0("tar -xvf  ",path.raw.data,"*.tar -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))



###test
# 
# # file.bc<-list.files(path=path.raw.data,pattern = "*.barcodes.tsv",recursive = TRUE,full.names = T)
# # file.bc
# # file.gene<-list.files(path=path.raw.data,pattern = "*genes.tsv",recursive = TRUE,full.names = T)
# # file.gene
# # file.mtx<-c(list.files(path=path.raw.data,pattern = "*matrix.mtx",recursive = TRUE,full.names = T),
# #             list.files(path=path.raw.data,pattern = "*exprMatrix.tsv",recursive = T,full.names = T))
# # file.mtx
# # file.meta<-list.files(path=path.raw.data,pattern = "*annotation.tsv",recursive = T,full.names = T)
# # file.meta
# # 
# # samplenames.1 <- lapply(file.mtx[1:2],
# #                       FUN = function(x){gsub(x = x,pattern = "_10X_matrix.mtx|GSE134722_FirstInstarLarvalBrain","", ignore.case = T) %>%
# #                           strsplit(split = "//") %>%
# #                           unlist() %>% .[2] }) %>% unlist()
# # 
# # 
# # samplenames.2<- lapply(file.mtx[3:4],
# #                       FUN = function(x){gsub(x = x,pattern = "_exprMatrix.tsv","", ignore.case = T) %>%
# #                           strsplit(split = "Metadata/") %>%
# #                           unlist() %>% .[2] }) %>% unlist()
# # 
# # samplenames<-c(samplenames.1,samplenames.2)
# # samplenames
# # # [1] "NormalCondition_finalaggr"     "StarvationCondition_finalaggr"
# # # [3] "Normal"                        "Merge"
# # 
# # 
# # sample.meta.names<-samplenames[3:4]
# # # [1] "Normal" "Merge"
# # 
# # 
# # 
# # ###read matrix
# # library(Matrix)
# # SCSeq.mtx<-list()
# # SCSeq.colnames.barcode<-list()
# # SCSeq.rownames.gene<-list()
# # 
# # for(i in 1:length(file.mtx)){
# #   if(i<=2){
# #     SCSeq.mtx[[i]]<- readMM(file.mtx[i])
# #     SCSeq.colnames.barcode[[i]]<- read.delim(file.bc[i],
# #                                         stringsAsFactors = F,header = F)
# #     SCSeq.rownames.gene[[i]]<- read.delim(file.gene[i],
# #                                      stringsAsFactors = F,header = F)
# #     
# #     colnames(SCSeq.mtx[[i]]) <- SCSeq.colnames.barcode[[i]]$V1
# #     rownames(SCSeq.mtx[[i]]) <- SCSeq.rownames.gene[[i]]$V2
# #   }else{
# #     SCSeq.mtx[[i]] <- read.delim(file.mtx[i],
# #                                       stringsAsFactors = F,header = T)
# #     rownames(SCSeq.mtx[[i]])<-SCSeq.mtx[[i]]$gene
# #     SCSeq.mtx[[i]]<-SCSeq.mtx[[i]][,-1]
# #     names(SCSeq.mtx[[i]])<-gsub(names(SCSeq.mtx[[i]]), pattern = "\\.", replacement = "-")
# #   }
# #   names(SCSeq.mtx)[i]<-samplenames[i]
# #   dim(SCSeq.mtx[[i]]) %>% print()
# # }
# # 
# # #     genes  cells
# # # [1] 17493  4708
# # # [1] 17493  4645
# # # [1] 12942  4407 ***normalized count matrix in the paper (Condition normal)
# # # [1] 14064  8696
# # 
# # # 8696-4407=4289
# # 
# # 
# # metadata<-list()
# # for(i in 1:length(file.meta)){
# #   metadata[[i]] <- read.delim(file.meta[i],
# #                             stringsAsFactors = F,header = T)
# #   names(metadata)[i]<-sample.meta.names[i]
# #   if(i==2){
# #     metadata[[i]]$id<-gsub(metadata[[i]]$Cell, pattern = "_.", replacement = "")
# #   }
# #   dim(metadata[[i]]) %>% print()
# # }
# # metadata[[2]]$id %>% unique()
# # 
# # #     cells
# # # [1] 4407    5
# # # [1] 8696    5
# # 
# # 
# # # save.image(file = paste0(path.script,"input.mtx.RData"))
# # # load(paste0(path.script,"input.mtx.RData"))
# # 
# # 
# # 
# # ### check matched data
# # for(i in 1:length(SCSeq.mtx)){
# #   if(i==1|i==3){
# #     intersect(colnames(SCSeq.mtx[[i]]), metadata[[1]]$Cell) %>%length()%>% print()
# #     }else if (i==2) {
# #       intersect(colnames(SCSeq.mtx[[i]]), metadata[[2]]$id) %>%length()%>% print()
# #     }else{
# #       intersect(colnames(SCSeq.mtx[[i]]), metadata[[2]]$Cell) %>%length()%>% print()
# #   }
# # }
# # 
# # # [1] 4407
# # # [1] 4471
# # # [1] 4407
# # # [1] 8696
# # 
# # 


# # ### Note: 
# # # NormalCondition_finalaggr.mtx - metadata[[1]]$Cell
# # # StarvationCondition_finalaggr.mtx - metadata[[2]]$Cell
# # # merged
# 


###restart
file.bc<-list.files(path=path.raw.data,pattern = "*.barcodes.tsv",recursive = TRUE,full.names = T)
file.bc
file.gene<-list.files(path=path.raw.data,pattern = "*genes.tsv",recursive = TRUE,full.names = T)
file.gene
file.mtx<-list.files(path=path.raw.data,pattern = "*matrix.mtx",recursive = TRUE,full.names = T)
file.mtx
file.meta<-list.files(path=path.raw.data,pattern = "*annotation.tsv",recursive = T,full.names = T)
file.meta

samplenames <- lapply(file.mtx,
                        FUN = function(x){gsub(x = x,pattern = "_finalaggr_10X_matrix.mtx|GSE134722_FirstInstarLarvalBrain","", ignore.case = T) %>%
                            strsplit(split = "//") %>%
                            unlist() %>% .[2] }) %>% unlist()
samplenames
# [1] "NormalCondition"     "StarvationCondition"

sample.meta.names<-lapply(file.meta,
                        FUN = function(x){gsub(x = x,pattern = "_annotation.tsv","", ignore.case = T) %>%
                            strsplit(split = "Metadata/") %>%
                             unlist() %>% .[2] }) %>% unlist()
sample.meta.names
# [1] "Normal" "Merge"



###read matrix
library(Matrix)
SCSeq.mtx<-list()
SCSeq.colnames.barcode<-list()
SCSeq.rownames.gene<-list()

for(i in 1:length(file.mtx)){
    SCSeq.mtx[[i]]<- readMM(file.mtx[i])
    SCSeq.colnames.barcode[[i]]<- read.delim(file.bc[i],
                                             stringsAsFactors = F,header = F)
    SCSeq.rownames.gene[[i]]<- read.delim(file.gene[i],
                                          stringsAsFactors = F,header = F)
    
    colnames(SCSeq.mtx[[i]]) <- SCSeq.colnames.barcode[[i]]$V1
    rownames(SCSeq.mtx[[i]]) <- SCSeq.rownames.gene[[i]]$V2
    names(SCSeq.mtx)[i]<-samplenames[i]
    dim(SCSeq.mtx[[i]]) %>% print()
}

#     genes  cells
# [1] 17493  4708
# [1] 17493  4645


metadata<-list()
for(i in 1:length(file.meta)){
  metadata[[i]] <- read.delim(file.meta[i],
                              stringsAsFactors = F,header = T)
  if(i==2){
    metadata[[i]]$id<-gsub(metadata[[i]]$Cell, pattern = "_.", replacement = "")
  }
  names(metadata)[i]<-sample.meta.names[i]
  dim(metadata[[i]]) %>% print()
}
#     cells
# [1] 4407    4
# [1] 8696    5


metadata[[2]] %>% str()


intersect(metadata[[2]]$Cell, paste0(colnames(SCSeq.mtx[[1]]),"_1")) %>% length()
# [1] 4349
intersect(metadata[[2]]$Cell, paste0(colnames(SCSeq.mtx[[2]]),"_2")) %>% length()
# [1] 4347 *** using this one
intersect(metadata[[2]]$id, colnames(SCSeq.mtx[[2]])) %>% length()
# [1] 4471 *** will have duplicated barcode
intersect(metadata[[2]]$Cell, paste0(colnames(SCSeq.mtx[[2]]),"_1")) %>% length()
# [1] 469


raw.merged.SCSeq.mtx<-list()
for (i in 1:length(SCSeq.mtx)) {
  raw.merged.SCSeq.mtx[[i]]<-SCSeq.mtx[[i]]
  colnames(raw.merged.SCSeq.mtx[[i]])<-paste0(colnames(SCSeq.mtx[[i]]),"_",i)
  tmp<-intersect(metadata[[2]]$Cell, colnames(raw.merged.SCSeq.mtx[[i]]))
  raw.merged.SCSeq.mtx[[i]]<-raw.merged.SCSeq.mtx[[i]][,tmp]
  print(dim(raw.merged.SCSeq.mtx[[i]]))
}
all.equal(rownames(raw.merged.SCSeq.mtx[[1]]),
          rownames(raw.merged.SCSeq.mtx[[2]]))

raw.metadata.3<-metadata[[2]][,c("Cell", "Cluster")]
raw.metadata.3 %>% head()




###final data using
metadata[[1]]<-metadata[[1]][,c("Cell", "Cluster")]
metadata[[1]]  %>% str()


colnames(SCSeq.mtx[[2]])<-paste0(colnames(SCSeq.mtx[[2]]),"_2")
metadata[[2]]<-metadata[[2]][,c("Cell", "Cluster")]
metadata[[2]]<-subset(metadata[[2]], Cell %in% colnames(SCSeq.mtx[[2]]))
metadata[[2]]  %>% str()


metadata[[3]] <- raw.metadata.3

for (i in 1:length(SCSeq.mtx)) {
  SCSeq.mtx[[i]] <- SCSeq.mtx[[i]][ , metadata[[i]]$Cell]
  SCSeq.mtx[[i]] %>% str() %>% print()
}
SCSeq.mtx[[3]]<-cbind(raw.merged.SCSeq.mtx[[1]],raw.merged.SCSeq.mtx[[2]])
names(SCSeq.mtx)[3]<-"Merged"

for (i in 1:length(SCSeq.mtx)) {
  all.equal(metadata[[i]]$Cell, colnames(SCSeq.mtx[[i]])) %>% print()
}

# save.image(file=paste0(path.script,GSE.id,".mtx.RData"))
# load(paste0(path.script,GSE.id,".mtx.RData"))


###normallize matrix
# i<-1
for (i in 3:length(SCSeq.mtx)) {
  logNormalizeMatrix <- log1p(sweep(SCSeq.mtx[[i]], 2, Matrix::colSums(SCSeq.mtx[[i]]), FUN = "/") * 10000)
  
  mean(SCSeq.mtx[[i]][,1])
  mean(logNormalizeMatrix[,1])
  write.csv(logNormalizeMatrix,
            file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
            row.names = F)


###annotate cluster to cell
str(metadata[[i]])
cell2clusterAssignment<-metadata[[i]]
dim(cell2clusterAssignment)
names(cell2clusterAssignment)<-c("Barcode","Cluster")

write.csv(cell2clusterAssignment,
          file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
          row.names = F)

###pseudoBulkMatrix
  mat<-SCSeq.mtx[[i]]
  df<- mat %>% as.matrix() %>% as.data.frame() %>% t()
  str(df)
  
  all.equal(cell2clusterAssignment$Barcode, row.names(df))
  
  tail(cell2clusterAssignment$Barcode)
  tail(row.names(df))
  
  df <- cbind(cell2clusterAssignment[ , c("Cluster"), drop = FALSE], df)
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
            file = paste0(path.result,names(SCSeq.mtx)[i],"_", PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))
  
  
  
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








