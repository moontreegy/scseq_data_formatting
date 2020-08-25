### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

# --------------------------see email 0825 from Ilyas
# Please find attached the barcodes2cluster file (Cosacaketal2019_barcodes2cellcluster.tab) for Cosacak et al. 2019.
# The file has 4 columns;
# 
# CellName: This is an arbitrary name we used. It has SampleName and the index number in the column. For instance, PBS_Cells_1, is the first cell in the raw data matrix, provided as count table on GEO (GSE118577) for PBS sample.
# barcode  : the barcode as colnames in the raw data matrix, provided as count table.
# SampleName: We have three samples; (i) PBS (Control samples), (ii) AB42 (amyloid-beta-42 injected), (iii) IL4 (interleukin-4 injected).
# ClusterID : The ClusterID/Name we in our paper. The ClusterID or names are as; NN (0-7), PC (0-8), OPCOD_1, OPCOD_2, Im.
# 
# I provided all cells and the cells that have been filtered out based on thresholds (nUMI, nGene, mitochondrial genes percentage) are given as "N/A".
# 
# I also downloaded and attached the count tables from GEO  for you (GSE118577_RAW.tar).


PMID<-"31018142"
Species<-"Zebrafish"
Tissue<-"Brain"
paper.title<-"Single-Cell Transcriptomics Analyses of Neural Stem Cell Heterogeneity and Contextual Plasticity in a Zebrafish Brain Model of Amyloid Toxicity"
GSE.id<-"GSE118577"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")
GEO.link<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id)
  
### create directory
# path.use<-"C:/Users/gaoyu/Downloads/"
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
# write(paste0("Note: I generated the results for 4 conditions. PBS(609 cells), AB42(737 cells), IL4(450 cells) and a merged version(1796 cells)\n"),file=file.rmd.name,append=TRUE)


# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))
# 

list.files(path=path.raw.data,recursive = TRUE,full.names = T)
files<-list.files(path=path.raw.data, pattern = "*.tab",recursive = TRUE,full.names = T)
files

file.name<-list.files(path=path.raw.data,pattern = "*.tab",recursive = TRUE,full.names = F)
file.name
file.name <- lapply(file.name[2:4], 
                    FUN = function(x){gsub(x = x,pattern = ".tab","", ignore.case = T) %>%
                        strsplit(split = "_") %>% 
                        unlist() %>% .[3] }) %>% unlist()
file.name



SCSeq.mtx<-list()
m<-1
for (i in 2:length(files)) {
  SCSeq.mtx[[m]]<-read.delim(files[i], stringsAsFactors = F, sep = " ")
  m<-m+1
}




# metadata
metadata.raw <- read.delim(files[1], stringsAsFactors = F)
# length(metadata.raw$CellName)
# length(metadata.raw$CellName %>% unique())
# dim(SCSeq.mtx[[1]])[2]+dim(SCSeq.mtx[[2]])[2]+dim(SCSeq.mtx[[3]])[2]
metadata.raw <- metadata.raw[metadata.raw$ClusterID!="N/A",]
length(metadata.raw$CellName)
### PBS 609,   AB42 737,  IL4 450   in the paper
609+737+450


file.name
metadata<-list()
for (i in 1:length(file.name)) {
  metadata[[i]]<-subset(metadata.raw, SampleName==file.name[i])
  metadata[[i]]<- metadata[[i]][,c(2,4)]
  names( metadata[[i]])<-c("Barcode", "Cluster")
}


for (i in 1:length(SCSeq.mtx)) {
  SCSeq.mtx[[i]] <- SCSeq.mtx[[i]][ , metadata[[i]]$Barcode]
  # SCSeq.mtx[[i]] %>% str() %>% print()
  all.equal(metadata[[i]]$Barcode, colnames(SCSeq.mtx[[i]])) %>% print()
}

SCSeq.mtx[[4]]<-cbind(SCSeq.mtx[[1]],SCSeq.mtx[[2]],SCSeq.mtx[[3]])
names(SCSeq.mtx[[4]])<-c(paste0(file.name[1],"_",colnames(SCSeq.mtx[[1]])),
                         paste0(file.name[2],"_",colnames(SCSeq.mtx[[2]])),
                         paste0(file.name[3],"_",colnames(SCSeq.mtx[[3]])))
names(SCSeq.mtx[[4]])
names(SCSeq.mtx)<-c(file.name,"mergerd")

metadata[[4]]<- rbind(metadata[[1]],metadata[[2]],metadata[[3]])
metadata[[4]]$Barcode<-c(paste0(file.name[1],"_",colnames(SCSeq.mtx[[1]])),
  paste0(file.name[2],"_",colnames(SCSeq.mtx[[2]])),
  paste0(file.name[3],"_",colnames(SCSeq.mtx[[3]])))
names(metadata)<-c(file.name,"mergerd")

all.equal(metadata[[4]]$Barcode, colnames(SCSeq.mtx[[4]])) %>% print()



# i<-1
for (i in 4:length(SCSeq.mtx)) {
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
