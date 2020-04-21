###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)

PMID<-"29622724"
Species<-"Mouse"
Tissue<-"kidney"
GSE.id<-"GSE107585"

foldername<-paste(PMID,Species,Tissue,sep="_")
system(paste0("mkdir ",foldername))
system(paste0("cd "foldername))
system(paste0("mkdir ",foldername,"/raw_data"))
system(paste0("mkdir ",foldername,"/results"))
system(paste0("mkdir ",foldername,"/script"))
system(paste0("touch ",foldername,"/README.md"))
system(paste0("touch ",foldername,"/script/formatting_",GSE.id,".R"))

system("tree ./")

# direc<- system("pwd", intern=T)
# setwd(direc)

direc_rawdata<-paste0("./",foldername,"/raw_data")
direc_results<-paste0("./",foldername,"/results")

getGEOSuppFiles(GSE.id,makeDirectory = F,baseDir=direc_rawdata)
system(paste0("gunzip ",direc_rawdata,"/GSE*"))
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

res.dt<-expressionMatrix
cols<-names(res.dt)[2:dim(res.dt)[2]]

for (j in cols){
  num.sum<-sum(res.dt[[j]])
  set(res.dt, j = j, value = log1p((res.dt[[j]]/num.sum)*10000))
}

mean(expressionMatrix$`AAACCTGAGATATGCA-1`)
mean(res.dt$`AAACCTGAGATATGCA-1`)


write.csv(dir=,res.dt,
		file=paste0(direc_results,"/",PMID,"_expressionMatrix_",Species,"_",Tissue,".csv"),
            row.names = F)
write.csv(cell2clusterAssignment,
		file=paste0(direc_results,"/",PMID,"_cell2clusterAssignment_",Species,"_",Tissue,".csv"),
            row.names = F)





