### Network properties
## read cor and p
#### Threshold 0.3
seed_cor_df_padj <- read.table("1_0.3_edge_seed.tsv", sep='\t', header =T)

nodeattrib_seed_combine <- data.frame(node=union(seed_cor_df_padj$Source,seed_cor_df_padj$Target))
nodeattrib_seed_combine$kingdom <- 0

for (i in as.character(nodeattrib_seed_combine$node))
{
  if (i%in%nodeattrib_seed_combine$node[(grep("^B", nodeattrib_seed_combine$node))] == TRUE)
  {nodeattrib_seed_combine[nodeattrib_seed_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_seed_combine[nodeattrib_seed_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_seed_combine) <- as.character(nodeattrib_seed_combine$node)


all_seed_net <- graph_from_data_frame(seed_cor_df_padj,direct=F, vertices=nodeattrib_seed_combine)

## Number of nodes
length(V(all_seed_net)) # 872

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_seed_net)))) #521
length(grep("^F",names(V(all_seed_net)))) #351


## Connections 
bb_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_seed) #1440

ff_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_seed) #1014

fb_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_seed) #2035



## Network properties
meta_degree <- sort(igraph::degree(all_seed_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "40"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "10.295871559633"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_seed_net, directed = FALSE)))
#"average shortest path length =   3.78939371596499"
print(paste("mean clustering coefficient = ", igraph::transitivity(all_seed_net, "global")))
#"mean clustering coefficient =  0.350295117152567"
print(paste("mean betweenness centrality = ", mean(betweenness(all_seed_net, directed = FALSE, normalized = TRUE))))
#"mean betweenness centrality = 0.00320619967352297"
print(paste("mean closeness centrality = ", mean(closeness(all_seed_net, normalized = TRUE))))
#"mean closeness centrality =  0.268180348643262"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_seed_net, V(all_seed_net))))))
#"mean number of neighbors = 10.295871559633"
##

net <- all_seed_net
seed_all_deg <- igraph::degree(net,mode="all")
seed_all_betweenness <- betweenness(net, normalized = TRUE)
seed_all_closeness <- closeness(net, normalized = TRUE)
seed_all_transitivity <- igraph::transitivity(net, "local", vids = V(net))
names(seed_all_transitivity)<- V(net)$name
seed_all_transitivity[is.na(seed_all_transitivity)] <- 0

### network properties of bacteria and fungi
df.seed.degree<-data.frame(seed_all_deg)
head(df.seed.degree)
df.seed.degree$Group <- "seed"
names(df.seed.degree)[1] <- c("Degree")

df.seed.closeness<-data.frame(seed_all_closeness)
head(df.seed.closeness)
df.seed.closeness$Group <- "seed"
names(df.seed.closeness)[1] <- c("Closeness")

df.seed.betweenness<-data.frame(seed_all_betweenness)
head(df.seed.betweenness)
df.seed.betweenness$Group <- "seed"
names(df.seed.betweenness)[1] <- c("Betweenness")

df.seed.degree$kingdom <- ifelse(grepl("^B",rownames(df.seed.degree)),'Bacteria', 'Fungi')
df.seed.closeness$kingdom <- ifelse(grepl("^B",rownames(df.seed.closeness)),'Bacteria', 'Fungi')
df.seed.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.seed.betweenness)),'Bacteria', 'Fungi')

shapiro.test(df.seed.degree$Degree) #W = 0.92132, p-value < 2.2e-16, not normal
shapiro.test(df.seed.betweenness$Betweenness) #W = 0.67267, p-value < 2.2e-16, not normal
shapiro.test(df.seed.closeness$Closeness) # = 0.98684, p-value = 4.711e-07, not normal

print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.seed.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.seed.betweenness)
p2

## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.seed.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()


#### Roles of prevalent ASVs and rare ASVs in the seed microbial network
### Ordering ASVs by prevalence values
data.b.2<- data.b[-c(71)]

data.b.sum<-data.frame(rowSums(data.b.2))
names(data.b.sum)[1] <- "sumPresence"
data.b.sum$Prevalence <- data.b.sum$sumPresence/70
data.b.sum$ID <- rownames(data.b.sum)

data.b.prev<- data.b.sum[c(2,3)]
data.b.prev <- data.b.prev%>%group_by(ID)%>% arrange(desc(Prevalence))


data.f.2<- data.f[-c(71)]

data.f.sum<-data.frame(rowSums(data.f.2))
names(data.f.sum)[1] <- "sumPresence"
data.f.sum$Prevalence <- data.f.sum$sumPresence/70
data.f.sum$ID <- rownames(data.f.sum)

data.f.prev<- data.f.sum[c(2,3)]
data.f.prev <- data.f.prev%>%group_by(ID)%>% arrange(desc(Prevalence))




##Over 80% prevalence
prev.b.list<-data.b.prev$ID[which(data.b.prev$Prevalence >= 0.8)]
prev.f.list<-data.f.prev$ID[which(data.f.prev$Prevalence >= 0.8)]
prev.ASV.list <- c(prev.b.list, prev.f.list)

### Less than 30% prevalence
rare.b.list<-data.b.prev$ID[which(data.b.prev$Prevalence <= 0.3)]
rare.f.list<-data.f.prev$ID[which(data.f.prev$Prevalence <= 0.3)]
rare.ASV.list <- c(rare.b.list, rare.f.list)

### Designate prevalent and rare ASVs in the network property results
df.seed.degree$nodeCate <- 0
for (i in rownames(df.seed.degree)){
  df.seed.degree$nodeCate[which(rownames(df.seed.degree)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
    }
  
df.seed.betweenness$nodeCate <- 0
for (i in rownames(df.seed.betweenness)){
  df.seed.betweenness$nodeCate[which(rownames(df.seed.betweenness)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

df.seed.closeness$nodeCate <- 0
for (i in rownames(df.seed.closeness)){
  df.seed.closeness$nodeCate[which(rownames(df.seed.closeness)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

ggplot(data=df.seed.degree, aes(x=nodeCate, y=Degree, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

ggplot(data=df.seed.betweenness, aes(x=nodeCate, y=Betweenness, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Betweenness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

ggplot(data=df.seed.closeness, aes(x=nodeCate, y=Closeness, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Closeness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")







#### Roles of prevalent ASVs and rare ASVs in the seed microbial network
### Ordering ASVs by prevalence values
data.b.2<- data.b[-c(71)]

data.b.sum<-data.frame(rowSums(data.b.2))
names(data.b.sum)[1] <- "sumPresence"
data.b.sum$Prevalence <- data.b.sum$sumPresence/70
data.b.sum$ID <- rownames(data.b.sum)

data.b.prev<- data.b.sum[c(2,3)]
data.b.prev <- data.b.prev%>%group_by(ID)%>% arrange(desc(Prevalence))


data.f.2<- data.f[-c(71)]

data.f.sum<-data.frame(rowSums(data.f.2))
names(data.f.sum)[1] <- "sumPresence"
data.f.sum$Prevalence <- data.f.sum$sumPresence/70
data.f.sum$ID <- rownames(data.f.sum)

data.f.prev<- data.f.sum[c(2,3)]
data.f.prev <- data.f.prev%>%group_by(ID)%>% arrange(desc(Prevalence))




##Over 80% prevalence
prev.b.list<-data.b.prev$ID[which(data.b.prev$Prevalence >= 0.8)]
prev.f.list<-data.f.prev$ID[which(data.f.prev$Prevalence >= 0.8)]
prev.ASV.list <- c(prev.b.list, prev.f.list)

### Less than 30% prevalence
rare.b.list<-data.b.prev$ID[which(data.b.prev$Prevalence <= 0.3)]
rare.f.list<-data.f.prev$ID[which(data.f.prev$Prevalence <= 0.3)]
rare.ASV.list <- c(rare.b.list, rare.f.list)

### Designate prevalent and rare ASVs in the network property results
df.seed.degree$nodeCate <- 0
for (i in rownames(df.seed.degree)){
  df.seed.degree$nodeCate[which(rownames(df.seed.degree)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

df.seed.betweenness$nodeCate <- 0
for (i in rownames(df.seed.betweenness)){
  df.seed.betweenness$nodeCate[which(rownames(df.seed.betweenness)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

df.seed.closeness$nodeCate <- 0
for (i in rownames(df.seed.closeness)){
  df.seed.closeness$nodeCate[which(rownames(df.seed.closeness)==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

ggplot(data=df.seed.degree, aes(x=nodeCate, y=Degree, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

ggplot(data=df.seed.betweenness, aes(x=nodeCate, y=Betweenness, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Betweenness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

ggplot(data=df.seed.closeness, aes(x=nodeCate, y=Closeness, fill=nodeCate))+ geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Closeness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")