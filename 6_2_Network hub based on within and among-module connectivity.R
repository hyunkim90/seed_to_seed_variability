### New hub estimation
#install.packages('brainGraph')
library(brainGraph)
library(igraph)
seed_cor_df_padj <- read.table("8_0.2247_edge_seed.tsv", sep='\t', header =T)

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
length(V(all_seed_net)) # 400

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_seed_net)))) #208
length(grep("^F",names(V(all_seed_net)))) #192


## Connections 
bb_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_seed) #358

ff_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_seed) #377

fb_occur_seed <- droplevels(seed_cor_df_padj[with(seed_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_seed) #373



## Network properties
meta_degree <- sort(igraph::degree(all_seed_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "26"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "5.54"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_seed_net, directed = FALSE)))
#"average shortest path length =   4.8224233784747"
print(paste("mean clustering coefficient = ", igraph::transitivity(all_seed_net, "global")))
#"mean clustering coefficient =  0.339733432752353"
print(paste("mean betweenness centrality = ", mean(betweenness(all_seed_net, directed = FALSE, normalized = TRUE))))
#"mean betweenness centrality = 0.00844268334151963"
print(paste("mean closeness centrality = ", mean(closeness(all_seed_net, normalized = TRUE))))
#"mean closeness centrality =  0.0318857910074282"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_seed_net, V(all_seed_net))))))
#"mean number of neighbors = 5.54"
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




network_hub_z_p<-function(net.soil.spar, keyword){####Assigning community membership (module)
  ## Perform cluster analysis using greedy clustering algorithm 
  cfg_soil <- cluster_fast_greedy(as.undirected(net.soil.spar))
  
  ### Among-module connectivity
  print(membership(cfg_soil))
  Pi<-part_coeff(net.soil.spar, membership(cfg_soil))
  z_p_table<-data.frame(Pi)
  head(data.frame(Pi))
  ### within-module connectivity
  Zi<-within_module_deg_z_score(net.soil.spar, membership(cfg_soil))
  z_p_table$zi <- Zi
  
  ###Assigning Kingdom
  z_p_table$kingdom <- 0
  z_p_table$kingdom[grep("^B",rownames(z_p_table))] <- "Bacteria"
  z_p_table$kingdom[grep("^F",rownames(z_p_table))] <- "Fungi"
  
  z_p_table$Role <- 0
  z_p_table$Role[z_p_table$zi >= 2.5 & z_p_table$Pi >= 0.62] <- "Network hub"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi < 0.62] <- "Peripheral"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi >= 0.62] <- "Connector"
  z_p_table$Role[z_p_table$zi >= 2.5 & z_p_table$Pi < 0.62] <- "Module hub"
  
  
  ### Assigning Module numbers
  cluster.raw<- membership(cfg_soil)
  
  df.cluster <- data.frame(matrix(ncol = 1, nrow = 0))
  names(df.cluster)[1] <- "Module"
  
  for(i in 1:length(cluster.raw)){
    df.add<-data.frame(cluster.raw[i])
    names(df.add)[1] <- "Module"
    df.cluster <- rbind(df.cluster, df.add)
  }
  
  z_p_table <- merge(z_p_table, df.cluster, by = "row.names")
  names(z_p_table)[1] <- "Node"
  write.csv(z_p_table, paste0("z_p_table_igraph_",keyword,".csv"))
  
  ##plotting
  p<-ggplot(z_p_table, aes(x=Pi, y=zi, color=kingdom)) +
    xlab('\n Among-module connectivity (Pi)')+
    ylab("Within-module connectivity (Zi) \n") +
    geom_point(size=3, alpha=0.7) +
    theme(aspect.ratio = 1)+
    # ggtitle("Volcano Plot \n") +
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size=FALSE) +
    #scale_x_continuous(breaks=seq(0,1,0.2))+
    #scale_y_continuous(breaks=seq(-20,0,-5))+
    scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
    geom_hline(yintercept=2.5, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', size = 0.75)
  
  return(p)
}


#### Hub identification with the SparCC-based network
network_hub_z_p(all_seed_net, "0.22_seed")

## Assigning rare and prev
z.p.table <- read.csv("z_p_table_igraph_0.22_seed.csv")
head(z.p.table)

z.p.table$NodeCate <- 0

for (i in as.character(z.p.table$Node)){
  z.p.table$NodeCate[which(z.p.table$Node==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

write.csv(z.p.table, "z_p_table_igraph_0.22_seed_prev_rare.csv")



z.p.table.SparCC <- read.csv("z_p_table_igraph_0.22_seed_prev_rare.csv")
z.p.table.SparCC$Color_dot <- factor(z.p.table.SparCC$Color_dot, levels = rev(c("Prev","Rare","Others","Peripheral")))
ggplot(z.p.table.SparCC, aes(x=Pi, y=zi, color=Color_dot, shape = Shape_dot)) +
  xlab('\n Among-module connectivity (Pi)')+
  ylab("Within-module connectivity (Zi) \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  #scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  scale_shape_manual(values=c("Keystone"=16,"Peripheral"=1)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=2.5, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', size = 0.75)

### connectivity of Prev, rare
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


#### Top network connectivity

df.seed.netProperties <- df.seed.degree 

for (i in rownames(df.seed.netProperties)){
  df.seed.netProperties$Betweenness[which(rownames(df.seed.netProperties) == i)] <- df.seed.betweenness$Betweenness[which(rownames(df.seed.betweenness) == i)]
  df.seed.netProperties$Closeness[which(rownames(df.seed.netProperties) == i)] <- df.seed.closeness$Closeness[which(rownames(df.seed.closeness) == i)]
}

write.xlsx(df.seed.netProperties, "Network properties_SparCC network_thres 0.2274.xlsx")
### 


### log normal distribution and hub

df.seed.netProperties$nodeCate <- factor(df.seed.netProperties$nodeCate, levels = c("Prev", "Rare", "Others"))
pD.B<-ggplot(df.seed.netProperties, aes(x=Betweenness, y=Degree, color=nodeCate)) +
  xlab('\n Betweenness')+
  ylab("Degree \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  #scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  #scale_shape_manual(values=c("Prev_hub"=17,"Rare_hub"=18, "Others" = 16)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=17, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.04, linetype='dashed', color='black', size = 0.75)

pD.C<-ggplot(df.seed.netProperties, aes(x=Closenness, y=Degree, color=nodeCate)) +
  xlab('\n Closeness')+
  ylab("Degree \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  #scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  #scale_shape_manual(values=c("Prev_hub"=17,"Rare_hub"=18, "Others" = 16)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=17, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.0353378797272164, linetype='dashed', color='black', size = 0.75)

library(gridExtra)
grid.arrange(pD.B, pD.C, nrow = 2)
df.seed.netProperties$Closenness[which(rownames(df.seed.netProperties)=="B242")]
min(df.seed.netProperties$Closenness)


pD.C2<-ggplot(subset(df.seed.netProperties, Closenness > 0.002518876), aes(x=Closenness, y=Degree, color=nodeCate)) +
  xlab('\n Closeness')+
  ylab("Degree \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  #scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  #scale_shape_manual(values=c("Prev_hub"=17,"Rare_hub"=18, "Others" = 16)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=17, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.0353378797272164, linetype='dashed', color='black', size = 0.75)



HighDegree<-rownames(df.seed.netProperties)[which(df.seed.netProperties$Degree > 17)]
HighBetween<-rownames(df.seed.netProperties)[which(df.seed.netProperties$Betweenness > 0.04)]
HighClose<-rownames(df.seed.netProperties)[which(df.seed.netProperties$Closenness > 0.0353378797272164)]

HighConnect<-Reduce(intersect, list(HighDegree,HighBetween, HighClose))
subset(df.seed.netProperties, rownames(df.seed.netProperties)%in% HighConnect)

install.packages("outliers")
library(outliers)
test <- dixon.test(df.seed.netProperties$Betweenness)
test

#### Hub identification with the Spearman-based network
network_hub_z_p(all_SPseed_net, "0.34_SPseed")

## Assigning rare and prev
prev.ASV.list
rare.ASV.list

z.p.table <- read.csv("z_p_table_igraph_0.34_SPseed.csv")
head(z.p.table)

z.p.table$nodeCate <- 0

for (i in as.character(z.p.table$Node)){
  z.p.table$nodeCate[which(z.p.table$Node==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

write.csv(z.p.table, "z_p_table_igraph_0.34_SPseed_prev_rare.csv")


#### Hub identification with the SparCC-based network
network_hub_z_p(all_seed_net, "0.15_seed")

## Assigning rare and prev
z.p.table <- read.csv("z_p_table_igraph_0.15_seed.csv")
head(z.p.table)

z.p.table$NodeCate <- 0

for (i in as.character(z.p.table$Node)){
  z.p.table$nodeCate[which(z.p.table$Node==i)] <- ifelse(i %in% prev.ASV.list, "Prev", ifelse(i %in% rare.ASV.list, "Rare","Others"))
}

write.csv(z.p.table, "z_p_table_igraph_0.15_seed_prev_rare.csv")





#### Plotting (different shapes for prevalent and rare nodes)
z.p.table.Spear <- read.csv("z_p_table_igraph_0.34_SPseed_prev_rare.csv")
z.p.table.Spear$Cate2 <- ifelse(z.p.table.Spear$Role %in% c("Connector","Module hub", "Network hub")& z.p.table.Spear$nodeCate == "Prev", "Prev_hub", 
                                ifelse(z.p.table.Spear$Role %in% c("Connector","Module hub", "Network hub")& z.p.table.Spear$nodeCate == "Rare", "Rare_hub","Others"))



z.p.table.SparCC <- read.csv("z_p_table_igraph_0.15_seed_prev_rare.csv")
z.p.table.SparCC$Cate2 <- ifelse(z.p.table.SparCC$Role %in% c("Connector","Module hub", "Network hub")& z.p.table.SparCC$nodeCate == "Prev", "Prev_hub", 
                                ifelse(z.p.table.SparCC$Role %in% c("Connector","Module hub", "Network hub")& z.p.table.SparCC$nodeCate == "Rare", "Rare_hub","Others"))


ggplot(z.p.table.Spear, aes(x=Pi, y=zi, color=kingdom, shape = Cate2)) +
  xlab('\n Among-module connectivity (Pi)')+
  ylab("Within-module connectivity (Zi) \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  scale_shape_manual(values=c("Prev_hub"=17,"Rare_hub"=18, "Others" = 16)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=2.5, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', size = 0.75)

ggplot(z.p.table.SparCC, aes(x=Pi, y=zi, color=kingdom, shape = Cate2)) +
  xlab('\n Among-module connectivity (Pi)')+
  ylab("Within-module connectivity (Zi) \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size="none") +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) +
  scale_shape_manual(values=c("Prev_hub"=17,"Rare_hub"=18, "Others" = 16)) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=2.5, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', size = 0.75)
