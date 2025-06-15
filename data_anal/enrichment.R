library(ggplot2)
library(reshape2)
library(dplyr)
root="G:\\enrichment"
data=data.frame()
for (file in list.files(root)) {
  temp=read.csv(paste0(root,"\\",file))
  temp=temp[temp$Raw.p<0.1,]
  data=rbind(data,data.frame((temp[,c(1,5)]),group=file))
}
data=data.frame(t(data))
data$Raw.p=-log10(data$Raw.p)
pathways=t(unique(data[1,]))
map=data.frame()
for (i in pathways) {
  tryCatch({
    pathwayid=rownames(data.frame(KEGGREST::keggFind("pathway",i)))[1]
    pathway_info=KEGGREST::keggGet(pathwayid)
    map=rbind(map,data.frame(i,strsplit(pathway_info[[1]]$CLASS,'; ')[[1]][2]))
  },error=function(e){
    print(i)
  })
}
colnames(map)=c("from","to")
map=rbind(map,data.frame(from='Alanine, aspartate and glutamate metabolism',to='Amino acid metabolism'))
map=rbind(map,data.frame(from='Valine, leucine and isoleucine biosynthesis',to='Amino acid metabolism'))
map=rbind(map,data.frame(from='Glycine, serine and threonine metabolism',to='Amino acid metabolism'))
map=rbind(map,data.frame(from='Phenylalanine, tyrosine and tryptophan biosynthesis',to='Amino acid metabolism'))
map=rbind(map,data.frame(from='Drug metabolism - other enzymes',to='Xenobiotics biodegradation and metabolism'))
map=rbind(map,data.frame(from='Glycolysis / Gluconeogenesis',to='Carbohydrate metabolism'))
map=rbind(map,data.frame(from='Neomycin, kanamycin and gentamicin biosynthesis',to='Biosynthesis of other secondary metabolites'))
map=rbind(map,data.frame(from='Valine, leucine and isoleucine degradation',to='Amino acid metabolism'))
data["supergroup",]=""
for (i in 1:ncol(data)) {
  temp=map[map$from==data[1,i],]
  data["supergroup",i]=temp$to[1]
}
data=data.frame(t(data))

plotdata=data.frame()
for (i in unique(data$supergroup)) {
  plotdata=rbind(plotdata,data[data$supergroup==i,])
}
for(i in 1:nrow(plotdata)){
  plotdata$groupx[i]=substr(as.character(plotdata$group[i]),1,nchar(as.character(plotdata$group[i]))-4)
}
plotdata$groupx <- factor(plotdata$groupx, levels = unique(plotdata$groupx))
plotdata$supergroup <- factor(plotdata$supergroup, levels = unique(plotdata$supergroup))
plotdata$X <- factor(plotdata$X, levels = unique(plotdata$X))
supergroup_colors <- setNames(rainbow(n=10),
                              levels(plotdata$supergroup))
color=c(rep(supergroup_colors[1],5),rep(supergroup_colors[2],13),rep(supergroup_colors[3],11),rep(supergroup_colors[4],4),
        rep(supergroup_colors[5],5),rep(supergroup_colors[6],2),rep(supergroup_colors[7],2),rep(supergroup_colors[8],2),
        rep(supergroup_colors[9],1),rep(supergroup_colors[10],1))
p <- ggplot(plotdata, aes(x = X, y = groupx, color = supergroup)) + 
  geom_point(
    aes(size = -log10(as.numeric(Raw.p))),
    alpha = 0.7
  ) +
  scale_size_continuous(
    name = "-log10(p)",
    range = c(1, 3)
  ) +
  labs(x = "Pathway", y = "Condition") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,color =color),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y.left = element_text(angle = 90, hjust = 1), 
    panel.spacing.y = unit(0.5, "lines") 
  ) +
  facet_grid(supergroup ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(name = "Supergroup", values = supergroup_colors);p 

ggsave(p, file = "D:\\grouped_colored_y.pdf", width = 12, height = 10)



