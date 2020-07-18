### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(reshape2)

PMID<-"30315278"
paper.title<-"The adult human testis transcriptional cell atlas"
Species<-"Human"
Tissue<-"Testis"
GSE.id<-"GSE112013"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")
GEO.link<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id)
  
### create directory
list.files()
# path.use <- "./"
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
write(paste0("Note: barcode2cluster is from supplementary table 1 in paper","\n"),file=file.rmd.name,append=TRUE)



getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
system(paste0("gunzip ",path.raw.data,"*.gz"))
# 
list.files(path=path.raw.data,recursive = TRUE,full.names = T)


SCSeq.mtx<- read.delim(paste0(path.raw.data,
                                   "GSE112013_Combined_UMI_table.txt"),
                            stringsAsFactors = F)
rownames(SCSeq.mtx)<-SCSeq.mtx$Gene
SCSeq.mtx<-SCSeq.mtx[,-1]
dim(SCSeq.mtx)
SCSeq.mtx[1:3, 1:3]


# metadata
metadata<- read.delim(paste0(path.raw.data,
                              "barcode2cluster.txt"),sep = " ",
                       stringsAsFactors = F, skip = 1)
str(metadata)

metadata.anno<-read.delim(paste0(path.raw.data,
                                 "cluster_annotation.txt"),
                          stringsAsFactors = F)
str(metadata.anno)

unique(metadata$Final_clusters)
unique(metadata.anno$Cluster.ID)

names(metadata.anno)[1]<-"Final_clusters"
metadata<-merge(metadata, metadata.anno, by="Final_clusters", all.x=T)
str(metadata)

names(SCSeq.mtx) %>% head()
metadata$CellID %>% head()
  
metadata$CellID <- gsub(metadata$CellID, pattern = "-",replacement = ".")

setdiff(names(SCSeq.mtx), metadata$CellID)
setdiff(metadata$CellID, names(SCSeq.mtx))

str(metadata)


###normallize matrix
logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)

mean(SCSeq.mtx[,1])
mean(logNormalizeMatrix[,1])
write.csv(logNormalizeMatrix,
          file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
          row.names = F)


###annotate cluster to cell
cell2clusterAssignment<-metadata[,c("CellID", "Cell.type.name", "Final_clusters")]
names(cell2clusterAssignment)[1:2]<-c("Barcode","Cluster")
write.csv(cell2clusterAssignment,
		file=paste0(path.result,PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)



###pseudoBulkMatrix
mat<-SCSeq.mtx
mat<-mat[, cell2clusterAssignment$Barcode]
df <- mat %>% as.matrix() %>% as.data.frame() %>% t()
df[, 1:3] %>% head()

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
          file = paste0(path.result, PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))



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





