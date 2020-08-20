### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"30479347"
Species<-"Drosophila"
Tissue<-"eye_disc"
paper.title<-"Single cell RNA-sequencing identifies a metabolic aspect of apoptosis in Rbf mutant"
GSE.id<-"GSE115476"
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
# write(paste0("Note: \n
# \tThe metadata got from: https://singlecell.broadinstitute.org/single_cell\n
# \tThere are two datasets in this paper, WT (11 samples) and RBF mutant (3 samples), GEO provides raw counts for all samples. \n
# \tThe paper analyse the WT alone and they also combined the WT and RBF to do analyse, so they have the 2 metadata tables, one for WT and another is for combined.\n
# \tAccording to metadata, WT (11,416 cells) and WT_RBF (5,589 wild type cells and 5,197 Rbf mutant cells)\n
# \tI generated 4 results, one for WT using WT metadata, other 3 for combined using combined metadata\n"),file=file.rmd.name,append=TRUE)


getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
system(paste0("gunzip ",path.raw.data,"*.gz"))


files<-list.files(path=path.raw.data, pattern = "*.txt",recursive = TRUE,full.names = T)
files
file.name<-list.files(path=path.raw.data,pattern = "*.txt",recursive = TRUE,full.names = F)
file.name <- lapply(file.name, 
                      FUN = function(x){gsub(x = x,pattern = ".expr.txt|.txt","", ignore.case = T) %>%
                          strsplit(split = "_") %>% 
                          unlist() %>% .[2] }) %>% unlist()
file.name


# metadata
metadata<-list()
metadata[[1]] <- read.delim(paste0(path.raw.data, "metaData_WT.txt"), stringsAsFactors = F,skip = 1)
metadata[[2]] <- read.delim(paste0(path.raw.data, "metaData_Combined_Mutant.txt"),stringsAsFactors = F,skip = 1)
names(metadata)<-c("WT", "Combined")

# metadata[[1]]$TYPE
# tmp<-metadata[[2]][metadata[[2]]$group=="WT",]
# tmp$bc<-gsub(tmp$TYPE, pattern = "WT_",replacement = "")
# setdiff(tmp$bc, metadata[[1]]$TYPE)
# setdiff(metadata[[1]]$TYPE, tmp$bc)
### ***two metadata's WT are not totally match, so process sepreately




###read matrix
#### 11 WT datasets
expr.wt.list<-list()
for (i in 1:11) {
  expr.wt.list[[i]]<-read.delim(files[i], header = T, stringsAsFactors = F)
  names(expr.wt.list)[i]<-file.name[i]
}
str(expr.wt.list)

###read matrix
#### 3 rbf datasets
expr.rbf.list<-list()
m<-1
for (i in 12:14) {
  expr.rbf.list[[m]]<-read.delim(files[i], header = T, stringsAsFactors = F)
  names(expr.rbf.list)[m]<-file.name[i]
  m<-m+1
}
str(expr.rbf.list)



# mat.wt
tmp.expr<-expr.wt.list[[1]]
names(tmp.expr)<-paste0(names(expr.wt.list)[1],"_",names(tmp.expr))
names(tmp.expr)[1]<-"GENE"

selected.col<-intersect(metadata[[1]]$TYPE,names(tmp.expr))
mat.wt<-tmp.expr[,c("GENE",selected.col)]
# i<-2
for (i in 2:length(expr.wt.list)) {
  tmp.expr<-expr.wt.list[[i]]
  names(tmp.expr)<-paste0(names(expr.wt.list)[i],"_",names(tmp.expr))
  names(tmp.expr)[1]<-"GENE"

  ####using metadata[[1]] for WT data
  selected.col<-intersect(metadata[[1]]$TYPE,names(tmp.expr))
  tmp<-tmp.expr[,c("GENE",selected.col)]
  mat.wt<-merge(mat.wt,tmp, by = "GENE", all=T)
}

mat.wt<-lapply(mat.wt, function(x) replace(x, is.na(x) ,0 )) %>% as.data.frame()
dim(mat.wt)




# mat.Combined.rbf
tmp.expr<-expr.rbf.list[[1]]
names(tmp.expr)<-paste0("rbf120a_",names(expr.rbf.list)[1],"_",names(expr.rbf.list[[1]]))
names(tmp.expr)[1]<-"GENE"

selected.col<-intersect(metadata[[2]]$TYPE,names(tmp.expr))
mat.Combined.rbf<-tmp.expr[,c("GENE",selected.col)]
# i<-2
for (i in 2:length(expr.rbf.list)) {
  tmp.expr<-expr.rbf.list[[i]]
  names(tmp.expr)<-paste0("rbf120a_",names(expr.rbf.list)[i],"_",names(expr.rbf.list[[i]]))
  names(tmp.expr)[1]<-"GENE"
  
  ####using metadata[[2]] for Combined.rbf data
  selected.col<-intersect(metadata[[2]]$TYPE,names(tmp.expr))
  tmp<-tmp.expr[,c("GENE",selected.col)]
  mat.Combined.rbf<-merge(mat.Combined.rbf,tmp, by = "GENE", all=T)
  }

mat.Combined.rbf<-lapply(mat.Combined.rbf, function(x) replace(x, is.na(x) ,0 )) %>% as.data.frame()
dim(mat.Combined.rbf)





# mat.Combined.wt
tmp.expr<-expr.wt.list[[1]]
names(tmp.expr)<-paste0("WT_",names(expr.wt.list)[1],"_",names(expr.wt.list[[1]]))
names(tmp.expr)[1]<-"GENE"

selected.col<-intersect(metadata[[2]]$TYPE,names(tmp.expr))
mat.Combined.wt<-tmp.expr[,c("GENE",selected.col)]
i<-8
for (i in 2:length(expr.wt.list)) {
  tmp.expr<-expr.wt.list[[i]]
  names(tmp.expr)<-paste0("WT_",names(expr.wt.list)[i],"_",names(expr.wt.list[[i]]))
  names(tmp.expr)[1]<-"GENE"
  
  ####using   metadata[[2]] for Combined.wt data
  selected.col<-intersect(metadata[[2]]$TYPE,names(tmp.expr))
  
  if(length(selected.col)!=0){
    tmp<-tmp.expr[,c("GENE",selected.col)]
    mat.Combined.wt<-merge(mat.Combined.wt,tmp, by = "GENE", all=T)
  }
}

mat.Combined.wt<-lapply(mat.Combined.wt, function(x) replace(x, is.na(x) ,0 )) %>% as.data.frame()
dim(mat.Combined.wt)

mat.Combined<-merge(mat.Combined.wt, mat.Combined.rbf, by = "GENE", all = T)
mat.Combined<-lapply(mat.Combined, function(x) replace(x, is.na(x) ,0 )) %>% as.data.frame()
dim(mat.Combined)






SCSeq.mtx<-list()
SCSeq.mtx[[1]]<- mat.wt
SCSeq.mtx[[2]]<- mat.Combined.wt
SCSeq.mtx[[3]]<- mat.Combined.rbf
SCSeq.mtx[[4]]<- mat.Combined
  
names(SCSeq.mtx)<-c("WT",
                    "Combined-WT",
                    "Combined-Rbf",
                    "Combined")

tmp<-metadata
metadata<-list()
metadata[[1]]<- subset(tmp[[1]], TYPE %in% names(SCSeq.mtx[[1]]))
metadata[[2]]<- subset(tmp[[2]][,c(1,3)], TYPE %in% names(SCSeq.mtx[[2]]))
metadata[[3]]<- subset(tmp[[2]][,c(1,3)], TYPE %in% names(SCSeq.mtx[[3]]))
metadata[[4]]<- subset(tmp[[2]][,c(1,3)], TYPE %in% names(SCSeq.mtx[[4]]))
names(metadata)<-c("WT",
                    "Combined-WT",
                    "Combined-Rbf",
                    "Combined")



i<-1
for (i in 1:length(SCSeq.mtx)) {
  names(metadata[[i]])<-c("Barcode","Cluster")
  rownames(SCSeq.mtx[[i]])<-SCSeq.mtx[[i]]$GENE
  SCSeq.mtx[[i]]<-SCSeq.mtx[[i]][,-1]
  SCSeq.mtx[[i]] <- SCSeq.mtx[[i]][ , metadata[[i]]$Barcode]
  # SCSeq.mtx[[i]] %>% str() %>% print()
  all.equal(metadata[[i]]$Barcode, colnames(SCSeq.mtx[[i]])) %>% print()
}

names(SCSeq.mtx)

str(SCSeq.mtx)
str(metadata)

# save(SCSeq.mtx,metadata, file=paste0(path.script,GSE.id,"_formatted.RData"))
# load(paste0(path.script,GSE.id,"_formatted.RData"))


###normallize matrix
# i<-1
for (i in 1:length(SCSeq.mtx)) {
  # logNormalizeMatrix <- log1p(sweep(SCSeq.mtx[[i]], 2, Matrix::colSums(SCSeq.mtx[[i]]), FUN = "/") * 10000)
  # 
  # mean(SCSeq.mtx[[i]][,1])
  # mean(logNormalizeMatrix[,1])
  # write.csv(logNormalizeMatrix,
  #           file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
  #           row.names = F)
  # 
  
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



