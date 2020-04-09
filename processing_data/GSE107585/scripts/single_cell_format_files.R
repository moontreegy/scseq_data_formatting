###R format GEO data
library(GEOquery)
library(magrittr)
library(data.table)

setwd("~/work/scseq/")

list.files()
getGEOSuppFiles('GSE107585')
system("gunzip ./GSE107585/GSE107585_Mouse_kidney_single_cell_datamatrix.txt.gz")

###read matrix
SCSeq.mtx<-fread("./GSE107585/GSE107585_Mouse_kidney_single_cell_datamatrix.txt")
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

head(expressionMatrix[,1:3])

res.dt<-expressionMatrix
cols<-names(res.dt)[2:dim(res.dt)[2]]

for (j in cols){
  num.sum<-sum(res.dt[[j]])
  set(res.dt, j = j, value = log1p((res.dt[[j]]/num.sum)*10000))
}

mean(res.dt$`AAACCTGAGATATGCA-1`)


write.csv(res.dt,
            "29622724_expressionMatrix_Mouse_kidney.csv",
            row.names = F)
write.csv(cell2clusterAssignment,
            "29622724_cell2clusterAssignment_Mouse_kidney.csv",
            row.names = F)




