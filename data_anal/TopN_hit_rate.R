library(tidyverse)
library(ggplot2)
drugs=list.files("D:\\swiss")
drugs=strsplit(drugs,"\\.")
methods=list.files("D:")
g2n=read.csv('D:\\MetTD\\EC2genename\\transfor.csv',row.names = T)
rownames(g2n)=g2n$enzyme
all_result=list()
n=30
for (m in methods) {
  if(m=="swiss"){
    all_result[[m]]=list()
    for (i in 1:13) {
      temp=read.csv(paste0("D:\\",m,"\\",drugs[[i]][1],".",drugs[[i]][2]))
      temp=temp[temp$Target.Class=="Enzyme",]
      all_result[[m]][[drugs[[i]][1]]]=temp$Common.name
    }
  }else if(m=="sea"){
    all_result[[m]]=list()
    for (i in 1:13) {
      temp=read.csv(paste0("D:\\",m,"\\",drugs[[i]][1],".",'xls'))
      temp=temp[do.call(rbind,lapply(temp$Target.ID,function(x){
        return(strsplit(x,"_")[[1]][2])
      }))=="HUMAN",]
      all_result[[m]][[drugs[[i]][1]]]=temp$Name
    }
  }else if(m=="stitch"){
    all_result[[m]]=list()
    for (i in 1:13) {
      all_result[[m]][[drugs[[i]][1]]]=logical(0)
      tryCatch({
        temp=read.table(paste0("D:\\",m,"\\",drugs[[i]][1],".",'tsv'),fill = TRUE,sep = "\t")
        temp=temp[temp[,2]==drugs[[i]][1],]
        all_result[[m]][[drugs[[i]][1]]]=temp[,1]
      },error=function(e){
      })
    }
  }else if(m=="PharmMapper"){
    for (i in 1:13) {
      temp=read.csv(paste0("D:\\",m,"\\",drugs[[i]][1],".",'csv'))
      temp=temp[!is.na(do.call(rbind,lapply(temp$Uniplot,function(x){
        return(strsplit(x,"_")[[1]][2])
      }))=="HUMAN"),]
      all_result[[m]][[drugs[[i]][1]]]=do.call(rbind,lapply(temp$Uniplot,function(x){
        return(strsplit(x,"_")[[1]][1])
      }))
      if(is.null(nrow(all_result[[m]][[drugs[[i]][1]]]))){
        all_result[[m]][[drugs[[i]][1]]]=logical(0)
      }
    }
  }else if(m=='NetInfer'){
    for (i in 1:13) {
      all_result[[m]][[drugs[[i]][1]]]=logical(0)
      temp=read.csv(paste0("D:\\",m,"\\",drugs[[i]][1],".",'csv'))
      all_result[[m]][[drugs[[i]][1]]]=temp$Gene.symbol
    }
  }else if(m=="MetTD"){
    for (i in 1:13) {
      all_result[[m]][[drugs[[i]][1]]]=logical(0)
      temp=readxl::read_excel(paste0("D:\\",m,'\\',drugs[[i]][1],'\\result.xlsx'))
      temp=temp[temp$p<0.1,]
      temp=temp[order(temp$p),]
      temp=g2n[temp$enzyme,-1]
      temp=temp[rowSums(!is.na(temp))>0,]
      all_result[[m]][[drugs[[i]][1]]]=temp
    }
  }else if(m=="BindingDB"){
    for (i in 1:13) {
      all_result[[m]][[drugs[[i]][1]]]=logical(0)
      temp=read.table(paste0("D:\\",m,"\\",drugs[[i]][1],".",'tsv'),header = TRUE,fill=T,sep = "\t")
      temp=temp[temp$Target.Source.Organism.According.to.Curator.or.DataSource=="Homo sapiens",]
      all_result[[m]][[drugs[[i]][1]]]=do.call(rbind,lapply(temp$UniProt..SwissProt..Entry.Name.of.Target.Chain,function(x){
        return(strsplit(x,"_")[[1]][1])
      }))
      all_result[[m]][[drugs[[i]][1]]]=all_result[[m]][[drugs[[i]][1]]][!is.na(all_result[[m]][[drugs[[i]][1]]])]
    }
  }
}

for (m in methods) {
  if(!m=="MetTD"){
    tnrate=rep(0,n)
    for (i in 1:13) {
      drug=drugs[[i]][1]
      target=read.csv(paste0("D:\\",drug,".csv"))
      target=unique(target$genename)
      testlist=all_result[[m]][[i]]
      for(tc in 1:min(20,length(testlist))){
        if(length(testlist>0)){
          if(testlist[tc]%in%target){
            tnrate=tnrate+c(rep(0,tc-1),rep(1,n-tc+1))
            all_result[["hit"]][[m]][[drug]]=tc
            break
          }
        }
      }
    }
    all_result[["Topnhitrate"]][[m]]=tnrate
  }else{
    tnrate=rep(0,n)
    for (i in 1:13) {
      drug=drugs[[i]][1]
      all_result[["hit"]][[m]][[drug]]=0
      target=read.csv(paste0("D:\\",drug,".csv"))
      target=unique(target$genename)
      testframe=all_result[[m]][[i]]
      for(tc in 1:min(n,nrow(testframe))){
        tl=testframe[tc,]
        tl=tl[!is.na(tl)]
        sig=sum(tl%in%target)
        if(sig>0){
          tnrate=tnrate+c(rep(0,tc-1),rep(1,n-tc+1))
          all_result[["hit"]][[m]][[drug]]=tc
          break
        }
      }
    }
    all_result[["Topnhitrate"]][[m]]=tnrate
  }
}
plotdata=data.frame()
for (m in methods) {
  plotdata=rbind(plotdata,all_result[['Topnhitrate']][[m]])
}
rownames(plotdata)=methods
colnames(plotdata)=1:30
plotdata=plotdata/13
plotdata_long <- plotdata %>%
  rownames_to_column(var = "Method") %>%
  pivot_longer(
    cols = -Method,
    names_to = "Category",
    values_to = "Value"
  )
plotdata_long$Category <- factor(plotdata_long$Category, levels = unique(plotdata_long$Category))
p=ggplot(plotdata_long, aes(x = Category, y = Value, group = Method, color = Method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  #geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Performance Comparison of Different Methods",
    x = "TopN",
    y = "Hitrate",
    color = "Method"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    legend.position = "top" 
  );p
ggsave("D:\\plot.pdf",width=16,height=6)

###########################
drug=drugs[[i]][1]
temp=readxl::read_excel(paste0("D:\\",m,"\\",drugs[[i]][1],"\\result.",'xlsx'))
enzymes=ec2gen(temp$enzyme)
geneframe=do.call(plyr::rbind.fill,lapply(enzymes,function(x){
  return(data.frame(t(x)))
}))
write.csv(geneframe,file=paste0("D:\\MetTD\\EC2genename\\transfor.csv"))
data=readxl::read_excel("D:\\MetTD\\EC2genename\\remain.xlsx")
ec2gen=function(x){
  a=lapply(x,function(x){
    ret=c()
    tryCatch({
      ret=KEGGREST::keggGet(x)
    }, error = function(e){
      print(x)
    })
    return(ret)
  })
  gene_kidney=lapply(a,function(x){
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
}









