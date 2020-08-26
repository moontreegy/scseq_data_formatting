### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"30726734"
Species<-"Human"
Tissue<-"Testis"
paper.title<-"The Neonatal and Adult Human Testis Defined at the Single-Cell Level"
GSE.id<-"GSE124263"
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


write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("GEO link: ", GEO.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Note: \n  
      \tThere are two data in this paper, adult (18723 cells), neonate (14862 cells), metadata has low- and high- resolution, I generated results according to different resolutions for adult and fetal respectively\n
      \tAll cell name I use short name, see details in metadata excel table\n"),
      file=file.rmd.name,append=TRUE)




getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
system(paste0("gunzip ",path.raw.data,"*.gz"))
# 
list.files(path=path.raw.data,recursive = TRUE,full.names = T)





# metadata
meta.file<-list.files(path=path.raw.data, pattern = "*.xlsx",recursive = TRUE,full.names = T)
meta.file
metadata.raw<-list()

for (i in 1:5) {
  metadata.raw[[i]] <- readxl::read_excel(meta.file, sheet = i)
}


### cell number
metadata.raw[[2]]$...1 %>% unique() %>% length() #Neonatal Cells 14862
# metadata.raw[[3]]$...1 %>% unique() %>% length() #Neonatal Germ Cells 179
metadata.raw[[4]]$...1 %>% unique() %>% length() #Adult Cells 18723
# metadata.raw[[5]]$...1 %>% unique() %>% length() #Adult SPG Cells 5839

### *** metadata.raw[[2]] [[4]] have all identities
metadata.raw[[1]]$`Cluster Identity`
# [1] "SC = Sertoli Cell"                     "LC = Leydig Cell"                     
# [3] "PTM = Peri Tubular Myoid"              "GC = Germ Cell"                       
# [5] "EC = Endothelial"                      "PGCL = Primordial Germ Cell Like"     
# [7] "SPG = Spermatogonia"                   "St = Spermatid"                       
# [9] "SPC = Spermatocyte"                    "M = Macrophage"                       
# [11] "SSC = Spermatogonial Stem Cell"        "TC = Transition Cell"                 
# [13] "Diff SPG = Differentiating SPG"        "Early SPG = Early differentiating SPG"



### cell clusters
### Sudhir et.al curated Neonate clusters:
### neoLC neoPTM PGCL PreSPG1 PreSPG2
names(metadata.raw[[2]])
metadata.Neonate<-merge(metadata.raw[[2]], metadata.raw[[3]], all=T, by="...1")
metadata.Neonate[is.na(metadata.Neonate$percent),]
metadata.Neonate$percent %>% unique()
### format name using legend
metadata.Neonate$percent<-gsub(metadata.Neonate$percent, pattern = " \\(\\S+\\)", replacement = "", perl = T)
metadata.Neonate$percent<-gsub(metadata.Neonate$percent, pattern = "/Blood", replacement = "", perl = T)
metadata.Neonate$percent %>% unique()
# [1] "GC"       "SC"       "EC/Blood" "PTM"      "LC" 
metadata.Neonate$NeoNatGerm %>% unique()
# [1] ***"PreSPG-2" NA         ***"PGCL"     ***"PreSPG-1"



### generate high resolution of metadata.Neonate
for (i in 1:length(metadata.Neonate$...1)) {
  if(!is.na(metadata.Neonate$NeoNatGerm[i])){
    metadata.Neonate$high.resolution[i]<-metadata.Neonate$NeoNatGerm[i]
  }else{
    metadata.Neonate$high.resolution[i]<-metadata.Neonate$percent[i]
  }
}



### Sudhir et.al curated adult clusters:
### adultLC PTM-NL adultPTM SSC-1 SSC-1A SSC-1B SSC-1C SSC-2 TC E.Diff DiffSPG
metadata.Adult<-merge(metadata.raw[[4]], metadata.raw[[5]][,c(1,4:6)], all=T, by="...1")

metadata.Adult[is.na(metadata.Adult$AdultCellAnnotation),]
unique(metadata.Adult$AdultCellAnnotation)
table(metadata.Adult$AdultCellAnnotation)
### format name using legend
metadata.Adult$AdultCellAnnotation<-
  stringr::str_match(metadata.Adult$AdultCellAnnotation, pattern = "\\(\\D+\\)")
unique(metadata.Adult$AdultCellAnnotation)
metadata.Adult$AdultCellAnnotation<-
  gsub(metadata.Adult$AdultCellAnnotation, pattern = "\\(|\\)", replacement = "")

metadata.Adult$AdultCellAnnotation[is.na(metadata.Adult$AdultCellAnnotation)]<-"EC"
metadata.Adult$AdultCellAnnotation<-as.character(metadata.Adult$AdultCellAnnotation)
table(metadata.Adult$AdultCellAnnotation)
unique(metadata.Adult$AdultCellAnnotation)
# LC    M  PTM  SPC  SPG   St   EC 
# 2500  553 4284 1169 6913  307 2997 


### format SPG_Subsets name using legend
metadata.raw[[1]]$`Cluster Identity`
unique(metadata.Adult$SPG_Subsets)
# [1] NA                    "Early Diff SPG"      "SSC1"                "Differentiating SPG"
# [5] "SSC2"
metadata.Adult$SPG_Subsets<-
  gsub(metadata.Adult$SPG_Subsets, pattern = "Early Diff SPG", replacement = "Early SPG", perl = T)
metadata.Adult$SPG_Subsets<-
  gsub(metadata.Adult$SPG_Subsets, pattern = "Differentiating SPG", replacement = "Diff SPG", perl = T)
unique(metadata.Adult$SPG_Subsets)


### format SSC1_Subsets name using legend
unique(metadata.Adult$SSC1_Subsets)
# [1] NA               "Early diff SPG" "SSC1-B"         "Diff SPG"       "SSC1-C"        
# [6] "SSC1-A"         "SSC2"           "Early Diff SPG"
metadata.Adult$SSC1_Subsets<-
  gsub(metadata.Adult$SSC1_Subsets, pattern = "Early Diff SPG", replacement = "Early SPG", perl = T, ignore.case = T)
unique(metadata.Adult$SSC1_Subsets)


### format TransitionCellIdentity name using legend
unique(metadata.Adult$TransitionCellIdentity)
# [1] NA                "Early diff SPG"  "SSC1"            "Diff SPG"       
# [5] "Transition Cell" "SSC2"
metadata.Adult$TransitionCellIdentity<-
  gsub(metadata.Adult$TransitionCellIdentity, pattern = "Early Diff SPG", replacement = "Early SPG", perl = T, ignore.case = T)
metadata.Adult$TransitionCellIdentity<-
  gsub(metadata.Adult$TransitionCellIdentity, pattern = "Transition Cell", replacement = "TC", perl = T, ignore.case = T)

unique(metadata.Adult$TransitionCellIdentity)

#test
# tmp<-subset(metadata.Adult, SPG_Subsets!=SSC1_Subsets|
#               SSC1_Subsets!=TransitionCellIdentity|
#               TransitionCellIdentity!=SPG_Subsets)

for (i in 1:length(metadata.Adult$orig.ident)) {
  if(!is.na(metadata.Adult$TransitionCellIdentity[i])&metadata.Adult$TransitionCellIdentity[i]=="TC"){
    metadata.Adult$high.resolution[i]<-metadata.Adult$TransitionCellIdentity[i]
  }else if(!is.na(metadata.Adult$SSC1_Subsets[i])){
    metadata.Adult$high.resolution[i]<-metadata.Adult$SSC1_Subsets[i]
  }else if(!is.na(metadata.Adult$SPG_Subsets[i])){
    metadata.Adult$high.resolution[i]<-metadata.Adult$SPG_Subsets[i]
  }else{
    metadata.Adult$high.resolution[i]<-metadata.Adult$AdultCellAnnotation[i]
  }
}





metadata<-list()
metadata[[1]]<-metadata.Adult[,c(1,4)]
unique(metadata[[1]]$AdultCellAnnotation) %>% length() #7
metadata[[2]]<-metadata.Adult[,c(1,8)]
unique(metadata[[2]]$high.resolution)# %>% length() #14
metadata[[3]]<-metadata.Neonate[,c(1,4)]
unique(metadata[[3]]$percent) %>% length() #5
metadata[[4]]<-metadata.Neonate[,c(1,6)]
unique(metadata[[4]]$high.resolution) %>% length() #8

names(metadata)<-c("LowRes-Adult","HighRes-Adult","LowRes-Neonate","HighRes-Neonate")


for (i in 1:length(metadata)) {
  names( metadata[[i]])<-c("Barcode", "Cluster")
}







### read matrix
files.mtx<-list.files(path=path.raw.data, pattern = "*matrix.mtx",recursive = TRUE,full.names = T)
files.mtx

files.genes<-list.files(path=path.raw.data, pattern = "*genes.tsv",recursive = TRUE,full.names = T)
files.genes

files.bc<-list.files(path=path.raw.data, pattern = "*barcodes.tsv",recursive = TRUE,full.names = T)
files.bc

# file.name<-list.files(path=path.raw.data,pattern = "*matrix.mtx",recursive = TRUE,full.names = F)
# file.name
# file.name <- lapply(file.name, 
#                     FUN = function(x){gsub(x = x,pattern = ".tsv","", ignore.case = T) %>%
#                         strsplit(split = "_") %>% 
#                         unlist() %>% .[2] }) %>% unlist()
# file.name
# # "d2"      "D2Total" "d7"      "D7total" "A1"      "A1Total" "A2"      "A2" 
file.name<-c("D2I", "D2T", "D7I" , "D7T", "A1I", "A1T", "A2I", "A2T" )




library(Matrix)
SCSeq.mtx.raw<-list()
i<-1
for (i in 1:length(files.mtx)) {
  
  SCSeq.mtx.raw[[i]]<- readMM(files.mtx[[i]])
  SCSeq.colnames.barcode<- read.delim(files.bc[[i]],
                                      stringsAsFactors = F,header = F)
  SCSeq.colnames.barcode %>% tail()
  SCSeq.colnames.barcode$V1<-gsub(SCSeq.colnames.barcode$V1, pattern = "-\\d", replacement = "")
  SCSeq.colnames.barcode$V1<-paste0(file.name[i],"_",SCSeq.colnames.barcode$V1)
  
  SCSeq.rownames.gene<- read.delim(files.genes[[i]],
                                   stringsAsFactors = F,header = F)
  
  dim(SCSeq.mtx.raw[[i]])
  SCSeq.mtx.raw[[i]][1:3, 1:3]
  colnames(SCSeq.mtx.raw[[i]]) <- SCSeq.colnames.barcode$V1
  rownames(SCSeq.mtx.raw[[i]]) <- SCSeq.rownames.gene$V2
  head(SCSeq.mtx.raw[[i]][, 1:3])
  
}
names(SCSeq.mtx.raw)<-file.name


#### In the paper, 14862 neonatal cells, 18723 adult cells
SCSeq.mtx.raw[[1]]@Dim[2]+SCSeq.mtx.raw[[2]]@Dim[2]+
  SCSeq.mtx.raw[[3]]@Dim[2]+SCSeq.mtx.raw[[4]]@Dim[2] ###14890

SCSeq.mtx.raw[[5]]@Dim[2]+SCSeq.mtx.raw[[6]]@Dim[2]+
  SCSeq.mtx.raw[[7]]@Dim[2]+SCSeq.mtx.raw[[8]]@Dim[2] ###19839

SCSeq.mtx.Neonate<-cbind(SCSeq.mtx.raw[[1]],SCSeq.mtx.raw[[2]],
                         SCSeq.mtx.raw[[3]],SCSeq.mtx.raw[[4]])
SCSeq.mtx.Adult<-cbind(SCSeq.mtx.raw[[5]],SCSeq.mtx.raw[[6]],
                       SCSeq.mtx.raw[[7]],SCSeq.mtx.raw[[8]])



SCSeq.mtx<-list()
SCSeq.mtx[[1]]<-SCSeq.mtx.Adult[, metadata.Adult$...1]
SCSeq.mtx[[2]]<-SCSeq.mtx.Neonate[, metadata.Neonate$...1]
names(SCSeq.mtx)<-c("Adult", "Neonate")



# for (i in 1:length(SCSeq.mtx)) {
# logNormalizeMatrix <- log1p(sweep(SCSeq.mtx[[i]], 2, Matrix::colSums(SCSeq.mtx[[i]]), FUN = "/") * 10000)
# mean(SCSeq.mtx[[i]][,1])
# mean(logNormalizeMatrix[,1])
# write.csv(logNormalizeMatrix,
#           file=paste0(path.result,names(SCSeq.mtx)[i],"_",PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
#           row.names = F)
# }


i<-1
for (i in 1:length(metadata)) {
  if(i<3){
    m<-1
  }else{
    m<-2
  }
  all.equal(metadata[[i]]$Barcode, colnames(SCSeq.mtx[[m]])) %>% print()

  ###annotate cluster to cell
  str(metadata[[i]])
  cell2clusterAssignment<-metadata[[i]]
  dim(cell2clusterAssignment)
  names(cell2clusterAssignment)<-c("Barcode","Cluster")
  
  write.csv(cell2clusterAssignment,
            file=paste0(path.result,names(metadata)[i],"_",PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)
  
  ###pseudoBulkMatrix
  mat<-SCSeq.mtx[[m]]
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
            file = paste0(path.result,names(metadata)[i],"_", PMID, "_pseudoBulkMatrix_", Species, "_", Tissue, ".csv"))
  
  
  
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
            file = paste0(path.result,names(metadata)[i],"_",PMID, "_clusterMetadataTable_", Species, "_", Tissue, ".csv"),
            row.names = FALSE)
  
  
  # extract data matrix from DotPlot function
  dot <- DotPlot(object = seuratObj, features = all_genes)
  gene2clusterTable <- dot$data
  colnames(gene2clusterTable) <- c("avg_exp", "pct_exp", "gene", "celltype_id", "avg_exp_scaled")
  
  # export gene2clusterTable
  write.csv(gene2clusterTable,
            file = paste0(path.result,names(metadata)[i],"_",  PMID, "_gene2clusterTable_", Species, "_", Tissue,".csv"),
            row.names = FALSE)
}

