### R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)

PMID<-"29980650"
Species<-"Human"
Tissue<-"Kidney"
paper.title<-"Single-Cell Transcriptomics of a Human Kidney Allograft Biopsy Specimen Defines a Diverse Inflammatory Response"
GSE.id.1<-"GSE109564"
GSE.id.2<-"GSE114156"
paper.link<-paste0("https://pubmed.ncbi.nlm.nih.gov/",PMID,"/")
GEO.link.1<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id.1)
GEO.link.2<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",GSE.id.2)

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
file.script.name<-paste0(path.script, "formatting_",GSE.id.1,"_", GSE.id.2,".R")

if (!file.exists(file.rmd.name)) {file.create(file.rmd.name)}
if (!file.exists(file.script.name)) {file.create(file.script.name)}


write(paste0("PMID: ", PMID,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Species: ", Species,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Tissue: ", Tissue,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper title: ", paper.title,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Paper link: ", paper.link,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("GEO link of allograft biopsy (single cell data): ", GEO.link.1,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("GEO link of healthy kidney (single nuclei data): ", GEO.link.2,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Result location: ", path.result,"\n"),file=file.rmd.name,append=TRUE)
write(paste0("Note:\n
             \tIn this paper they have two samples, healthy adult kidney (single-nuclei Seq, 4259 cells) and a single kidney transplant biopsy core (single-cell Seq, 4487 cells)\n,
             \tBesides, they compare the epithelia between healthy- and allograft sample, so I generated 5 kinds of results:\n
             \t\thealthy kidney: 4259 nucleis with 6 clusters, results with prefix 'healthy'\n
             \t\tallograft kidney: 4487 cells with 16 clusters, results with prefix 'allograft'\n
             \t\tcombined epithelia: 4643 cells (healthy 3532 cells, allograft 1111 cells), results with prefix 'Combined-Epithelia-...' \n"),file=file.rmd.name,append=TRUE)

# Comparison of healthy kidney epithelial transcriptomes with biopsy specimen counterparts identified novel segment-specific proinflammatory responses in rejection. Endothelial cells formed three distinct subclusters: resting cells and two activated endothelial cell groups.



# getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir= path.raw.data)
# system(paste0("tar -xvf  ",path.raw.data,"GSE* -C", path.raw.data))
# system(paste0("gunzip ",path.raw.data,"*.gz"))
# 


### wget https://www.dropbox.com/s/z7m8hfqfqwtmx6d/biopsy.seurat.Robj?dl=0
### wget https://www.dropbox.com/s/seeyu0aju5xl56e/MTS_kidney2_ycut2_pc8_merged.Robj?dl=0
### wget https://www.dropbox.com/s/3oxkzddud4i9tdg/epithelia_compare.rds?dl=0

list.files(path=path.raw.data,recursive = TRUE,full.names = T)

load(paste0(path.raw.data,"biopsy.seurat.Robj"))
load(paste0(path.raw.data,"MTS_kidney2_ycut2_pc8_merged.Robj"))
epithelia.comb<-readRDS(paste0(path.raw.data,"epithelia_compare.rds"))




### read data for allograft biopsy
meta.allograft<-biopsy@ident %>% as.data.frame()
meta.allograft$Barcode<-rownames(meta.allograft)
names(meta.allograft)[1]<-"Cluster"
meta.allograft$Cluster<-as.character(meta.allograft$Cluster)
meta.allograft<-meta.allograft[,c(2,1)]
unique(meta.allograft$Cluster)
# [1] "LOH (DL)"      "Plasma2"       "T cells"       "EC"            "CD"           
# [6] "PT"            "LOH (AL)"      "Mono1"         "Pericyte"      "Mono2"        
# [11] "Myofibroblast" "Plasma1"       "B cells"       "Fibroblast"    "Cycling"      
# [16] "Mast cells" 

SCSeq.mtx.allograft<-biopsy@raw.data
setdiff(colnames(SCSeq.mtx.allograft),meta.allograft$Barcode)
setdiff(meta.allograft$Barcode,colnames(SCSeq.mtx.allograft))
SCSeq.mtx.allograft<-SCSeq.mtx.allograft[,meta.allograft$Barcode]
dim(SCSeq.mtx.allograft)
# [1] 20477  4487


###read data for healthy kidney
meta.healthy<-MTS_kidney2@ident %>% as.data.frame()
meta.healthy$Barcode<-rownames(meta.healthy)
names(meta.healthy)[1]<-"Cluster"
meta.healthy$Cluster<-as.character(meta.healthy$Cluster)
meta.healthy<-meta.healthy[,c(2,1)]
unique(meta.healthy$Cluster)
#[1] "LH"    "PT"    "CD:PC" "DT"    "P"     "CD:IC"


SCSeq.mtx.healthy<-MTS_kidney2@raw.data
setdiff(colnames(SCSeq.mtx.healthy),meta.healthy$Barcode)
setdiff(meta.healthy$Barcode,colnames(SCSeq.mtx.healthy))
SCSeq.mtx.healthy<-SCSeq.mtx.healthy[,meta.healthy$Barcode]
dim(SCSeq.mtx.healthy)
# [1] 21384  4259



###read data for epithelia_compare
meta.comb<-epithelia.comb@meta.data %>% as.data.frame()
meta.comb$orig.ident %>% unique()
meta.comb$Barcode<-rownames(meta.comb)
names(meta.comb)[1]<-"Cluster"
meta.comb$Cluster<-as.character(meta.comb$Cluster)
meta.comb<-meta.comb[,c(10,1,2)]
# meta.comb$Barcode<-gsub(meta.comb$Barcode, pattern = "\\.\\S+_",  replacement = "_", perl = T)
meta.comb$Cluster<-gsub(meta.comb$Cluster, pattern = "\\S+\\.",  replacement = "", perl = T)
unique(meta.comb$Cluster)
#[1] "LH"    "PT"    "CD:PC" "DT"    "P"     "CD:IC"


SCSeq.mtx.comb<-epithelia.comb@assays$RNA@counts
setdiff(colnames(SCSeq.mtx.comb),meta.comb$Barcode)
setdiff(meta.comb$Barcode,colnames(SCSeq.mtx.comb))
SCSeq.mtx.comb<-SCSeq.mtx.comb[,meta.comb$Barcode]
dim(SCSeq.mtx.comb)
# [1] 21384  4259

meta.comb.healthy<-subset(meta.comb, protocol=="MTS")
SCSeq.mtx.comb.healthy<-SCSeq.mtx.comb[,meta.comb.healthy$Barcode]

meta.comb.allograft<-subset(meta.comb, protocol=="Bx")
SCSeq.mtx.comb.allograft<-SCSeq.mtx.comb[,meta.comb.allograft$Barcode]




SCSeq.mtx<- list()

SCSeq.mtx[[1]]<-SCSeq.mtx.allograft
SCSeq.mtx[[2]]<-SCSeq.mtx.healthy
SCSeq.mtx[[3]]<-SCSeq.mtx.comb
SCSeq.mtx[[4]]<-SCSeq.mtx.comb.allograft
SCSeq.mtx[[5]]<-SCSeq.mtx.comb.healthy

names(SCSeq.mtx)<-c("allograft", "healthy", "Combined-Epithelia", "Combined-Epithelia-allograft","Combined-Epithelia-healthy")


# metadata
metadata <- list()

metadata[[1]]<-meta.allograft
metadata[[2]]<-meta.healthy
metadata[[3]]<-meta.comb[,1:2]
metadata[[4]]<-meta.comb.allograft[,1:2]
metadata[[5]]<-meta.comb.healthy[,1:2]

names(metadata)<-c("allograft", "healthy", "Combined-Epithelia", "Combined-Epithelia-allograft","Combined-Epithelia-healthy")


# save(SCSeq.mtx,file=paste0(path.script,GSE.id,".mtx.RData"))
# load(paste0(path.script,GSE.id,".mtx.RData"))




###normallize matrix
i<-1
for (i in 1:length(SCSeq.mtx)) {
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




