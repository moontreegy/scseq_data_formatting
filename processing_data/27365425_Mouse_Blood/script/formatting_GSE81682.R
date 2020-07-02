### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(reshape2)

PMID<-"27365425"
paper.title<-"A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation"
Species<-"Mouse"
Tissue<-"Blood"
GSE.id<-"GSE81682"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")
GEO.link<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id)
  
### create directory
list.files()
# path.use <- "/Users/gaoyu/Downloads/"
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
# write(paste0("Filterd normalized data (1654 cells) and annotation are download from:","http://blood.stemcells.cam.ac.uk/single_cell_atlas.html",
#              "(Nestorowa et al., 2016)","\n"),file=file.rmd.name,append=TRUE)



getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
system(paste0("gunzip ",path.raw.data,"*.gz"))
# 
list.files(path=path.raw.data,recursive = TRUE,full.names = T)



###http://blood.stemcells.cam.ac.uk/single_cell_atlas.html#data
SCSeq.mtx.norm<- read.delim(paste0(path.raw.data,
                                   "coordinates_gene_counts_flow_cytometry.txt"),
                            stringsAsFactors = F)
# Download matrix of normalised gene counts and flow cytometry protein levels
# here. The first column contains the cell name, the next 3 columns contain the
# diffusion map coordinates where the cell is located; the following 10 columns
# contain flow cytometry data of 9 cell-surface markers and FSC-H; the rest of
# the columns contain gene expression data labelled with ensembl gene IDs.
# If you would like to locate a particular cell in the diffusion map, introduce
# the cell name here.

SCSeq.mtx.anno<- read.delim(paste0(path.raw.data,"all_cell_types.txt"),
                                 stringsAsFactors = F,header = T)
# To find cell type identity download the cell type matrix here. Cell names
# match the above matrices, and the cell types are given for both the broad and
# narrow cell type categories that were retrospectively identified using the
# index data. An entry of "1" in the table indicates a cell belongs to the cell
# type category, whereas a "0" means it does not. Due to the presence of both
# narrow and broad gates in the table, along with cell types that are subsets
# of other cell types (e.g. MPP1 and MPP), cells can belong to multiple categories.


SCSeq.mtx.anno$cell.Name<-gsub(rownames(SCSeq.mtx.anno),pattern = "-",replacement = ".")

setdiff(SCSeq.mtx.anno$cell.Name, SCSeq.mtx.norm$cell.Name)
setdiff(SCSeq.mtx.norm$cell.Name, SCSeq.mtx.anno$cell.Name )


SCSeq.mtx.anno<-subset(SCSeq.mtx.anno, cell.Name %in% SCSeq.mtx.norm$cell.Name)


SCSeq.mtx.anno.broad<-melt(SCSeq.mtx.anno[,-c(11:23)])
SCSeq.mtx.anno.broad<-subset(SCSeq.mtx.anno.broad, value==1)
SCSeq.mtx.anno.broad<-aggregate(variable~cell.Name,data = SCSeq.mtx.anno.broad,
                          function(x){paste(unique(x),collapse = ";")})
names(SCSeq.mtx.anno.broad)<-c("Barcode", "Annotaion_broad")

SCSeq.mtx.anno.narrow<-melt(SCSeq.mtx.anno[,-c(1:10,23)])
SCSeq.mtx.anno.narrow<-subset(SCSeq.mtx.anno.narrow, value==1)
SCSeq.mtx.anno.narrow<-aggregate(variable~cell.Name,data = SCSeq.mtx.anno.narrow,
                                function(x){paste(unique(x),collapse = ";")})
names(SCSeq.mtx.anno.narrow)<-c("Barcode", "Annotaion_narrow")



SCSeq.mtx.raw<- read.delim(paste0(path.raw.data,"GSE81682_HTSeq_counts.txt"))
setdiff(SCSeq.mtx.norm$cell.Name, colnames(SCSeq.mtx.raw))

SCSeq.mtx<-SCSeq.mtx.raw[, SCSeq.mtx.norm$cell.Name]
rownames(SCSeq.mtx)<-SCSeq.mtx.raw$ID
colnames(SCSeq.mtx)
dim(SCSeq.mtx)
SCSeq.mtx[1:3, 1:3]


# metadata
metadata <- data.frame(Barcode=colnames(SCSeq.mtx), stringsAsFactors = F)
metadata$Cluster <- gsub(metadata$Barcode, pattern = "_\\S+",replacement = "", perl = T)
metadata<-merge(metadata, SCSeq.mtx.anno.broad, by = "Barcode", all.x = T)
metadata<-merge(metadata, SCSeq.mtx.anno.narrow, by = "Barcode", all.x = T)

unique(metadata$Cluster)
str(metadata)


###normallize matrix
logNormalizeMatrix <- log1p(sweep(SCSeq.mtx, 2, Matrix::colSums(SCSeq.mtx), FUN = "/") * 10000)

mean(SCSeq.mtx[,1])
mean(logNormalizeMatrix[,1])
write.csv(logNormalizeMatrix,
          file=paste0(path.result,PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
          row.names = F)


###annotate cluster to cell
cell2clusterAssignment<-metadata
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





