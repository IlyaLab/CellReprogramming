
library(plyr)
library(stringr)
library(Rtsne)
library(scRecover)


getExp=function (cells, dem_){ 
  names(cells)=NULL
  if(length(cells)>0){
    cells=str_replace(cells, "-", ".")
    dem_data=read.csv2(dem_, sep="\t",  header=TRUE, row.names=1,stringsAsFactors = TRUE,  as.is=TRUE) 
    data_final=as.matrix(dem_data[,colnames(dem_data)%in%cells])
    colnames(data_final)=colnames(dem_data)[colnames(dem_data)%in%cells]
    return(data_final)
    }
}

cell_files<-function(path){
   setwd(path)
   files=list.files("./", recursive=TRUE)
   anno_files=files[which(str_detect(files, 'anno')==TRUE)]
   dem_files=files[which(str_detect(files, 'dem')==TRUE)]
   prob_gene_exp=list()
   mono_gene_exp=list()
  for (anno in anno_files){
    print(anno)
    start_ind=str_locate(anno, '_')
    end_ind=str_locate(anno, '\\.')
    sample_id=substr(anno,start_ind, end_ind)
    matched_dem=str_detect(dem_files, fixed(sample_id))
    dem_=dem_files[matched_dem]
    anno_data=read.csv2(anno, sep="\t")
    cells_prob=unlist(subset(anno_data, CellType=='ProB', select=c(Cell)))
    prob_gene_exp[[anno]]=getExp(cells_prob,dem_)
    cells_mono=unlist(subset(anno_data, CellType=='Mono', select=c(Cell)))
    mono_gene_exp[[anno]]=getExp(cells_mono, dem_)
  }
   setwd("../.")
   return(list(prob_exp=prob_gene_exp, mono_exp=mono_gene_exp))
}

result=cell_files("Healthy")
prob_exp_matrix=do.call(cbind, result$prob_exp)
mono_exp_matrix=do.call(cbind, result$mono_exp)

all_data=cbind(prob_exp_matrix,mono_exp_matrix)

sample_label=rbind(cbind(colnames(prob_exp_matrix), c(rep.int("ProB", ncol(prob_exp_matrix)))),cbind(colnames(mono_exp_matrix), c(rep.int("Mono", ncol(mono_exp_matrix)))))
labels=sample_label[,2]

colnames(sample_label)=c("CellID","CellType")
write.table(sample_label, file="labels.txt",  quote=FALSE, col.names=TRUE, row.names=FALSE)

library(BiocParallel)
param <- MulticoreParam(workers = 18, progressbar = TRUE)
register(param)
scRecover(all_data, labels=labels, outputDir="./imputed/", parallel = TRUE, BPPARAM = param)

