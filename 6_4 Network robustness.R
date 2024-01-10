#### Contribution of prevalence and rare nodes in the robustness of the network
#(1) Extinction scheme in the order of highest Degree to lowest
links <- read.table("8_0.2247_edge_seed.tsv", sep='\t', header =T)
web <- graph.data.frame(links, directed = F)

mytheme_2d <- theme_bw() + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "top")

# function for getting largest connected component
largest_component <- function(wild_web){
  components <- igraph::clusters(wild_web)
  biggest_cluster_id <- which.max(components$csize)
  # ids
  vert_ids <- V(wild_web)[components$membership == biggest_cluster_id]
  # subgraph
  largest_component <- igraph::induced_subgraph(wild_web, vert_ids)
  return(largest_component)
}

# largest component 
lg_web <- largest_component(web)
# vertex count
vcount(lg_web) #375

plot(lg_web, vertex.size=5, vertex.label=NA, vertex.color='lightblue', edge.arrow.width=0.3, edge.arrow.curve=0.5)

# get vertex list of the largest component and get the degree for those nodes
lg_cen <- tibble::rownames_to_column(df.seed.netProperties, var="Node") %>% filter(Node %in% V(lg_web)$name)

write.csv(lg_cen,"centrality table.csv")

# Set extinction order
lg_cen<-read.csv("centrality table.csv",header = T) ### Degree

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction, "Exinction table.csv")
# plot 
deg.p<-ggplot(df_extinction, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)

#####Betweenness
lg_cen.btw<-read.csv("centrality table.btw.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen.btw%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.btw <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.btw, "Exinction table.btw.csv")
# plot 
btw.p<-ggplot(df_extinction.btw, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)

### Closeness
lg_cen.close<-read.csv("centrality table.close.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen.close%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.close <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.close, "Exinction table.close.csv")
# plot 
close.p<-ggplot(df_extinction.close, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)


### Remove rare first
lg_cen<-read.csv("./robustness/remove rare node first/centrality table.rare.csv",header = T) ### Degree

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.rare <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.rare, "Exinction table.rare.csv")
# plot 
deg.rare.p<-ggplot(df_extinction.rare, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.03733, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.832, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)

#####Betweenness
lg_cen.btw<-read.csv("./robustness/remove rare node first/centrality table.btw.rare.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen.btw%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.btw.rare <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.btw.rare, "Exinction table.btw.rare.csv")
# plot 
btw.rare.p<-ggplot(df_extinction.btw.rare, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.03467, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.832, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)

### Closeness
lg_cen.close<-read.csv("./robustness/remove rare node first/centrality table.close.rare.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare, and remaining nodes
ext_order.test <- lg_cen.close%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.close.rare <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.close.rare, "Exinction table.close.rare.csv")
# plot 
close.p.rare<-ggplot(df_extinction.close.rare, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.03467, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.832, linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept=0.84667, linetype='dashed', color='black', size = 0.75)


### Compare
df_extinction.deg.PrevRare<-read.csv("./robustness/remove rare node first/Exinction table.degree.csv",header = T) ### Betweenness

ggplot(df_extinction.deg.PrevRare, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)

df_extinction.btw.PrevRare<-read.csv("./robustness/remove rare node first/Exinction table.betweenness.csv",header = T) ### Betweenness

ggplot(df_extinction.btw.PrevRare, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)

df_extinction.close.PrevRare<-read.csv("./robustness/remove rare node first/Exinction table.closeness.csv",header = T) ### Betweenness

ggplot(df_extinction.close.PrevRare, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.04267, linetype='dashed', color='black', size = 0.75)
  

AUC(df_extinction.close.PrevRare$Removed[which(df_extinction.close.PrevRare$Attack_order == "Rare_first")],
    df_extinction.close.PrevRare$Lg_size[which(df_extinction.close.PrevRare$Attack_order == "Rare_first")])

AUC(df_extinction.close.PrevRare$Removed[which(df_extinction.close.PrevRare$Attack_order == "Prev_first")],
    df_extinction.close.PrevRare$Lg_size[which(df_extinction.close.PrevRare$Attack_order == "Prev_first")])

AUC(df_extinction.deg.PrevRare$Removed[which(df_extinction.deg.PrevRare$Attack_order == "Rare_first")],
    df_extinction.deg.PrevRare$Lg_size[which(df_extinction.deg.PrevRare$Attack_order == "Rare_first")])

AUC(df_extinction.deg.PrevRare$Removed[which(df_extinction.deg.PrevRare$Attack_order == "Prev_first")],
    df_extinction.deg.PrevRare$Lg_size[which(df_extinction.deg.PrevRare$Attack_order == "Prev_first")])

AUC(df_extinction.btw.PrevRare$Removed[which(df_extinction.btw.PrevRare$Attack_order == "Rare_first")],
    df_extinction.btw.PrevRare$Lg_size[which(df_extinction.btw.PrevRare$Attack_order == "Rare_first")])

AUC(df_extinction.btw.PrevRare$Removed[which(df_extinction.btw.PrevRare$Attack_order == "Prev_first")],
    df_extinction.btw.PrevRare$Lg_size[which(df_extinction.btw.PrevRare$Attack_order == "Prev_first")])


#### Attack order (from largest to smallest values)

lg_cen<-read.csv("./robustness/remove prevalent node first/centrality table.csv",header = T) ### Degree
lg_cen <- lg_cen[sample(1:dim(lg_cen)[1]),] # shuffling because there are nodes that have same degree
ext_order <- lg_cen %>% arrange(desc(Degree)) %>% dplyr::select(Node)
ext_order <- ext_order$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtractedeb <- lg_web - vertices(ext_order[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtractedeb)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.deg.arrange <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
df_extinction.deg.arrange$Attack_order <- "HighToLow"

### Betweenness
ext_order <- lg_cen %>% arrange(desc(Betweenness)) %>% dplyr::select(Node)
ext_order <- ext_order$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtractedeb <- lg_web - vertices(ext_order[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtractedeb)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.btw.arrange <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
df_extinction.btw.arrange$Attack_order <- "HighToLow"


### Closeness
ext_order <- lg_cen %>% arrange(desc(Closeness)) %>% dplyr::select(Node)
ext_order <- ext_order$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtractedeb <- lg_web - vertices(ext_order[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtractedeb)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.close.arrange <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
df_extinction.close.arrange$Attack_order <- "HighToLow"


df_extinction.deg.all<-rbind(df_extinction.deg.PrevRare,df_extinction.deg.arrange)
df_extinction.btw.all<-rbind(df_extinction.btw.PrevRare,df_extinction.btw.arrange)
df_extinction.close.all<-rbind(df_extinction.close.PrevRare,df_extinction.close.arrange)


ggplot(df_extinction.deg.all, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d


ggplot(df_extinction.btw.all, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

ggplot(df_extinction.close.all, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

### Remove rare nodes last
### Remove rare first
lg_cen<-read.csv("./robustness/remove rare node first/centrality table.rare.last.csv",header = T) ### Degree

#Extinction order arranged by prevalent, rare.last, and remaining nodes
ext_order.test <- lg_cen%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.rare.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.rare.last, "Exinction table.rare.last.csv")
# plot 
deg.rare.last.p<-ggplot(df_extinction.rare.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

#####Betweenness
lg_cen.btw<-read.csv("./robustness/remove rare node first/centrality table.btw.rare.last.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare.last, and remaining nodes
ext_order.test <- lg_cen.btw%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.btw.rare.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.btw.rare.last, "Exinction table.btw.rare.last.csv")
# plot 
btw.rare.last.p<-ggplot(df_extinction.btw.rare.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

### Closeness
lg_cen.close<-read.csv("./robustness/remove rare node first/centrality table.close.rare.last.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, rare.last, and remaining nodes
ext_order.test <- lg_cen.close%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.close.rare.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.close.rare.last, "Exinction table.close.rare.last.csv")
# plot 
close.p.rare.last<-ggplot(df_extinction.close.rare.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

### Remove prevalent nodes last
lg_cen<-read.csv("./robustness/remove prevalent node first/centrality table.prev.last.csv",header = T) ### Degree

#Extinction order arranged by prevalent, prev.last, and remaining nodes
ext_order.test <- lg_cen%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.prev.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.prev.last, "Exinction table.prev.last.csv")
# plot 
deg.prev.last.p<-ggplot(df_extinction.prev.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

#####Betweenness
lg_cen.btw<-read.csv("./robustness/remove prevalent node first/centrality table.btw.prev.last.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, prev.last, and remaining nodes
ext_order.test <- lg_cen.btw%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.btw.prev.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.btw.prev.last, "Exinction table.btw.prev.last.csv")
# plot 
btw.prev.last.p<-ggplot(df_extinction.btw.prev.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d

### Closeness
lg_cen.close<-read.csv("./robustness/remove prevalent node first/centrality table.close.prev.last.csv",header = T) ### Betweenness

#Extinction order arranged by prevalent, prev.last, and remaining nodes
ext_order.test <- lg_cen.close%>%dplyr::select(Node)
ext_order.test <- ext_order.test$Node

# set vectors that will become x and y of the 2D plot
lg_size_vec <- c(vcount(lg_web))
remove_num_vec <- 0:vcount(lg_web)

# Execute extinction
for (i in 1:vcount(lg_web)){
  # print(i)
  # print(ext_order[i])
  # subtract the vertex in order one by one
  subtracted_web <- lg_web - vertices(ext_order.test[1:i])
  # calculate the largest component
  subtracted_lg <- largest_component(subtracted_web)
  # print(paste0('size of largest component: ',vcount(subtracted_lg)))
  # append the size of largest component
  # if the size of the largest component is 1, it is regarded as 0 because no links exist.
  if (vcount(subtracted_lg) == 1){
    lg_size_vec <- c(lg_size_vec,0)
  } else{
    lg_size_vec <- c(lg_size_vec,vcount(subtracted_lg))
  }
}

lg_size_vec %>% length()
remove_num_vec %>% length()

# create dataframe
df_extinction.close.prev.last <- tibble(Removed= remove_num_vec/vcount(lg_web), Lg_size = lg_size_vec/vcount(lg_web))
write.csv(df_extinction.close.prev.last, "Exinction table.close.prev.last.csv")
# plot 
close.p.prev.last<-ggplot(df_extinction.close.prev.last, aes(x=Removed, y=Lg_size))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d


df_extinction.prev.last$Attack_order <- "Prev_last"
df_extinction.btw.prev.last$Attack_order <- "Prev_last"
df_extinction.close.prev.last$Attack_order <- "Prev_last"

df_extinction.rare.last$Attack_order <- "Rare_last"
df_extinction.btw.rare.last$Attack_order <- "Rare_last"
df_extinction.close.rare.last$Attack_order <- "Rare_last"


df_extinction.deg.all.last<-rbind(df_extinction.prev.last,df_extinction.rare.last,df_extinction.deg.arrange)
df_extinction.btw.all.last<-rbind(df_extinction.btw.prev.last,df_extinction.btw.rare.last,df_extinction.btw.arrange)
df_extinction.close.all.last<-rbind(df_extinction.close.prev.last,df_extinction.close.rare.last,df_extinction.close.arrange)

ggplot(df_extinction.deg.all.last, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Degree-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.12533, linetype='dashed', color='blue', size = 0.75)+ ## prev starts to be removed 
  geom_vline(xintercept=0.168, linetype='dashed', color='blue', size = 0.75)+## rare starts to be removed
  geom_vline(xintercept=0.12533, linetype='dashed', color='green', size = 0.75)+##rare starts to be removed
  geom_vline(xintercept=0.95733, linetype='dashed', color='green', size = 0.75) ## prev starts to be removed 



ggplot(df_extinction.btw.all.last, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Betweenness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.12533, linetype='dashed', color='blue', size = 0.75)+ ## prev starts to be removed 
  geom_vline(xintercept=0.168, linetype='dashed', color='blue', size = 0.75)+## rare starts to be removed
  geom_vline(xintercept=0.12533, linetype='dashed', color='green', size = 0.75)+##rare starts to be removed
  geom_vline(xintercept=0.95733, linetype='dashed', color='green', size = 0.75) ## prev starts to be removed 


ggplot(df_extinction.close.all.last, aes(x=Removed, y=Lg_size, color = Attack_order))+ #, group=name)) +
  geom_line(size=0.5) +
  ggtitle(paste0("Closeness-based attack \n")) +
  xlab('\n Fraction of vertices removed')+
  ylab("Fractional size of largest component  \n") +
  mytheme_2d+
  geom_vline(xintercept=0.12533, linetype='dashed', color='blue', size = 0.75)+ ## prev starts to be removed 
  geom_vline(xintercept=0.168, linetype='dashed', color='blue', size = 0.75)+## rare starts to be removed
  geom_vline(xintercept=0.12533, linetype='dashed', color='green', size = 0.75)+##rare starts to be removed
  geom_vline(xintercept=0.95733, linetype='dashed', color='green', size = 0.75) ## prev starts to be removed 


write.csv(df_extinction.close.all.last, "raw data for robustness_closeness part.csv")
write.csv(df_extinction.btw.all.last, "raw data for robustness_betweenness part.csv")
write.csv(df_extinction.deg.all.last, "raw data for robustness_degree part.csv")
