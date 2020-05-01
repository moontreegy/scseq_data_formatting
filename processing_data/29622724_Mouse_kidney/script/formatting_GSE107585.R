###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"29622724"
Species<-"Mouse"
Tissue<-"kidney"
GSE.id<-"GSE107585"


### create directory
foldername<-paste(PMID,Species,Tissue,sep="_")
# system(paste0("mkdir ",foldername))
# system(paste0("cd ", foldername))
# system(paste0("mkdir ",foldername,"/raw_data"))
# system(paste0("mkdir ",foldername,"/results"))
# system(paste0("mkdir ",foldername,"/script"))
# system(paste0("touch ",foldername,"/README.md"))
# system(paste0("touch ",foldername,"/script/formatting_",GSE.id,".R"))


# direc<- system("pwd", intern=T)
# setwd(direc)

direc_rawdata<-paste0("./",foldername,"/raw_data/")
direc_results<-paste0("./",foldername,"/results/")
setwd("/n/groups/flyrnai/Yue/scseq_data_formatting/processing_data/")
system("tree ./")


# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir=direc_rawdata)
# system(paste0("gunzip ",direc_rawdata,"GSE*"))
filename<-list.files(path=direc_rawdata,pattern = "GSE*",recursive = TRUE,full.names = T)



###read matrix
SCSeq.mtx<-fread(filename)
class(SCSeq.mtx)
head(SCSeq.mtx[,1:3])

# save(SCSeq.mtx,file="GSE107585.mtx.RData")
# load("./GSE107585.mtx.RData")


###annotate cluster to cell
cell2clusterAssignment<-data.frame(Barcode=names(SCSeq.mtx),
                                   Cluster=t(SCSeq.mtx[1,]))
cell2clusterAssignment<-cell2clusterAssignment[-1,]

###normallize matrix
expressionMatrix<-SCSeq.mtx[-1,]
names(expressionMatrix)[1]<-"Gene"

head(expressionMatrix[,1:5])

cols<-names(expressionMatrix)[2:dim(expressionMatrix)[2]]

for (j in cols){
  num.sum<-sum(expressionMatrix[[j]])
  set(expressionMatrix, j = j, value = log1p((expressionMatrix[[j]]/num.sum)*10000))
}

mean(expressionMatrix$`AAACCTGAGATATGCA-1`)


###export 2 tables
write.csv(expressionMatrix,
		file=paste0(direc_results,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
            row.names = F)
write.csv(cell2clusterAssignment,
		file=paste0(direc_results,PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)



###pseudoBulkMatrix
mat<-SCSeq.mtx
df <- transpose(mat, make.names = "V1")
df[1:3, 1:3]

# sum gene counts by each celltype
df_group_by_celltype <- df %>%
  group_by(Cluster_Number) %>%
  summarise_all(sum) %>%
  as.data.frame()
class(df_group_by_celltype)

df_group_by_celltype <- transpose(df_group_by_celltype, 
                                  keep.names = "GeneName",
                                  make.names = "Cluster_Number")
df_group_by_celltype[1:3, ]
str(df_group_by_celltype)


### export pseudoBulkMatrix
pseudoBulkMatrix <- df_group_by_celltype
write.csv(pseudoBulkMatrix,
          file = paste0(direc_results, PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))



### dot plot
mat<-SCSeq.mtx[-1,]
mat<-as.data.frame(mat)
mat[1:3,1:3]
names(mat)
rownames(mat)<-mat$V1
mat<-mat[,-1]

library(Matrix)
sparse.Mat<-as.matrix(mat)
sparse.Mat <- Matrix(sparse.Mat, sparse = TRUE)

seuratObj <- CreateSeuratObject(counts = sparse.Mat)
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

cluster.anno<-read.delim(paste0(direc_results,"annotate_cluster_name.txt"),
                         header = F)
names(cluster.anno)<-c("celltype_id","cell_name")

clusterMetadataTable<-merge(clusterMetadataTable, cluster.anno)

# export clusterMetadataTable
write.csv(clusterMetadataTable,
          file = paste0(direc_results, PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)

# extract data matrix from DotPlot function
dot <- DotPlot(object = seuratObj, features = all_genes)
gene2clusterTable <- dot$data
colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")

# export gene2clusterTable
write.csv(gene2clusterTable,
          file = paste0(direc_results, PMID, "_gene2clusterTable_", Species, "_", Tissue, ".csv"),
          row.names = FALSE)




