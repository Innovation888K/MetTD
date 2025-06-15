load("D:\\RSpace\\MetTD\\human_model.Rda")
load("D:\\RSpace\\MetTD\\hmdb_ms1_database0.0.3.rda")
load("D:\\RSpace\\MetTD\\Mixed_db_final.rda")
local_path = "D:\\RSpace\\MetTD\\"
source(paste(local_path, "MetTD.R", sep = ""))
library(rsbml)
library(parallel)
library(dplyr)
library(readr)
library(KEGGREST)
library(igraph)
library(ggplot2)
library(ggrepel)
library(future)
library(future.apply)
library(openxlsx)
library(tidymass)
library(progressr)
library(metpath)
filter='(M+H)+'
filepath="D:\\"
for(i in list.files(filepath)){
root=paste0(filepath,'\\',i)
aux_root=root
drug_num='result'
polar='positive'
  metabolism_annotation(
    path = root,
    polar = polar,
    database = hmdb_ms1_database0.0.3,
    filter = filter,
    min_fraction=0
  )
  if(!is.na(filter)){
    x=readxl::read_excel(paste0(root,"\\Result\\filted.xlsx"))
  }else{
    x=readxl::read_excel(paste0(root,"\\Result\\All_annotations.xlsx"))
  }
  express=read.csv(paste0(root,"\\Result\\knn_expression_data.csv"))
  rownames(express)=express$X
  express=express[x$variable_id,]
  x$Level=express
  temp=express
  model_peak = x$Level[, 2:((ncol(express)+1)/2)]
  baseline =x$Level[,(((ncol(express)+1)/2)+1):(ncol(express))]
  bl=lapply(1:nrow(baseline),function(x){
    rdt=t(baseline[x,])
    std=sd(rdt)
    mean=mean(rdt)
    ret=rdt[rdt>(mean-3*std)&rdt<mean+3*std]
    return(mean(ret))
  })
  bl=do.call(rbind,bl)
  mp=lapply(1:nrow(baseline),function(x){
    rdt=t(model_peak[x,])
    std=sd(rdt)
    mean=mean(rdt)
    ret=rdt[rdt>(mean-1.96*std)&rdt<mean+1.96*std]
    return(mean(ret))
  })
  mp=do.call(rbind,mp)
  ml=(mp-bl)
  model_level = (model_peak - baseline)
  model_level = data.frame(Level = rowMeans(model_level))
  model_level=ml
  temp = dplyr::select(x, -'Level')
  temp$Level = ml
  tryCatch({
    temp$Level=do.call(rbind,temp$Level)
  },error=function(e){
    
  })
  
  temp = temp[!is.na(temp$Level), ]
  Model$ref = get_ref_idx(temp, Model, default_annotation_database)
  Model$local_score = get_local_score(Model, do_zscore = T)
  p = get_p(Model, 500, do_zscore = T,n_cores=8)
  Model$p = p
  scatter_plot(Model, 500, first_n = 5, p_limit = 0.5)
  ggsave("result.pdf",width=17,height=10)
  setwd(root)
  LPS = data.frame(
    p = p,
    ls = Model$local_score,
    enzyme = Model$enzyme$Compound.name
  )
  write.xlsx(LPS, file = paste0(root, '\\', drug_num, '.xlsx'))
  sample_data=LPS
  sample_data=sample_data[order(t(sample_data[,"score"]), decreasing = T),]
  sample_data=sample_data[order(t(sample_data[,'p'])), ]
  
  a=lapply(sample_data$enzyme[1:30],function(x){
    ret=c()
    tryCatch({
      ret=keggGet(x)
    }, error = function(e){
      print(x)
    })
    return(ret)
  })
  gene=lapply(a,
              function(x){
                G=x[[1]]$GENES[1]
                if(!is.null(G)){
                  if(substr(G,1,3)=="HSA"){
                    left=c()
                    right=c()
                    for (i in 5:nchar(G)) {
                      if(substr(G,i,i)=="("){
                        left=c(left,i)
                      }else if(substr(G,i,i)==")"){
                        right=c(right,i)
                      }
                    }
                    ret=c()
                    for (i in 1:length(left)) {
                      ret=c(ret,substr(G,left[i]+1,right[i]-1))
                    }
                    return(ret)
                  }
                }
                return(NA)
              }
  )
  ge=c()
  for (x in gene) {
    ge=c(ge,x)
  }
  ge=ge[!is.na(ge)]
  ge=unique(ge)
  write.csv(ge,file="first30.csv")
}
  
