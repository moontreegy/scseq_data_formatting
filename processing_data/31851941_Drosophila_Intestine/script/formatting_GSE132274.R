###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"31851941"
Species<-"Drosophila"
Tissue<-"Intestine"
GSE.id<-"GSE132274"
paper.title<-"The Cellular Diversity and Transcription Factor Code of Drosophila Enteroendocrine Cells"
paper.link<-"https://pubmed.ncbi.nlm.nih.gov/31851941/"
GEO.link<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132274"
  
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


# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))
# 
list.files(path=path.raw.data,recursive = TRUE,full.names = T)


###read matrix
library(Matrix)
SCSeq.mtx<- readMM(paste0(path.raw.data,"GSM3855437_matrix.mtx"))
SCSeq.colnames.barcode<- read.delim(paste0(path.raw.data,"GSM3855437_barcodes.tsv"),
                                    stringsAsFactors = F,header = F)
SCSeq.rownames.gene<- read.delim(paste0(path.raw.data,"GSM3855437_genes.tsv"),
                                 stringsAsFactors = F,header = F)

dim(SCSeq.mtx)
SCSeq.mtx[1:3, 1:3]
colnames(SCSeq.mtx) <- SCSeq.colnames.barcode$V1
rownames(SCSeq.mtx) <- SCSeq.rownames.gene$V2
head(SCSeq.mtx[, 1:3])

# metadata
metadata <- read.csv(paste0(path.raw.data, "cell.metadata.csv"), header = TRUE,stringsAsFactors = F)
str(metadata)
unique(metadata$res.0.4)
metadata$id <- paste0(metadata$CellBarcode, "-1")

intersect(metadata$id, SCSeq.colnames.barcode$V1)

SCSeq.mtx <- SCSeq.mtx[ , metadata$id]
dim(SCSeq.mtx)


class(SCSeq.mtx)
head(SCSeq.mtx[,1:3])

# save(SCSeq.mtx,file=paste0(path.script,GSE.id,".mtx.RData"))
# load(paste0(path.script,GSE.id,".mtx.RData"))

###normallize matrix
logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)

mean(SCSeq.mtx[,1])
mean(logNormalizeMatrix[,1])
write.csv(logNormalizeMatrix,
          file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
          row.names = F)


###annotate cluster to cell
cell2clusterAssignment<-metadata[,c(7,6)]
names(cell2clusterAssignment)<-c("Barcode","Cluster")
write.csv(cell2clusterAssignment,
		file=paste0(path.result,PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)



###pseudoBulkMatrix
mat<-SCSeq.mtx
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





