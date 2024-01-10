#### Phylogenetic tree with prevalence information


# install.packages(c("tidytree", "ggstar", "ggnewscale", "TDbook"))
# BiocManager::install("treeio")
# BiocManager::install("ggtreeExtra")
# BiocManager::install("ggtree")

library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

# load data from TDbook, including tree_hmptree, 
# df_tippoint (the abundance and types of microbes),
# df_ring_heatmap (the abundance of microbes at different body sites),
# and df_barplot_attr (the abundance of microbes of greatest prevalence)

### Bacteria
sample_data(bac.clean.ss.f)
tree <- phyloseq::phy_tree(bac.clean.ss.f)

## Data 1: Tip information (taxonomic information)
bac.list.2 <- bac.list
bac.list.2$Phylum[is.na(bac.list.2$Phylum)] <- "Unidentified"
bac.list.2$Class[is.na(bac.list.2$Class)] <- "Unidentified"
bac.list.2$Phylum2 <- bac.list.2$Phylum
bac.list.2$Phylum2[which(bac.list.2$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Phylum == "Proteobacteria" & bac.list.2$Class == "Unidentified")] <- "Unidentified Proteobacteria"

bac.list.2[-c(2:10)]
data1 <- bac.list.2[-c(2:10)]
names(data1)[1]<-"ID"

## Data 2: Specificity of ASVs (Heat map)
bac.clean.preval <- bac.clean.ss.f
bac.melt.preval<-psmelt(bac.clean.preval)
bac.melt.preval$Presence <- ifelse(bac.melt.preval$Abundance > 0, 1, 0)
bac.melt.preval <- bac.melt.preval%>% group_by(OTU, Branch_number) %>% summarise(sumPresence=sum(Presence))
names(bac.melt.preval)[1] <- "ID"
data2 <- bac.melt.preval
data2$sumPresence <- as.character(data2$sumPresence)
data2$sumPresence <- factor(data2$sumPresence, levels = c("0","1","2","3","4","5","6","7","8","9","10"))

## Data 3: Prevalence data of each ASVs in all 70 seeds (bar plot)
otu.G <- otu_table(bac.clean.ss.f)
df.otu.G <- data.frame(otu.G)
head(df.otu.G)

df.otu.G[df.otu.G>0] <- 1

df.otu.sum<-data.frame(rowSums(df.otu.G))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$ID <- rownames(df.otu.sum)

data3<- df.otu.sum[c(2,3)]

bac.prev.tab<-data3
names(bac.prev.tab)[2] <- "OTU" 
bac.prev.tab<-merge(bac.list, bac.prev.tab, by = c("OTU" = "OTU"))
write.csv(bac.prev.tab, "bac.prev.tab.csv")

#write.csv(bac.melt.preval, "bac.melt.preval.csv")

tree <- phyloseq::phy_tree(bac.clean.ss.f)

nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
#nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
# poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
#              1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
#              0.3, 0.4, 0.3)
#labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)
#labdf <- data.frame(node=nodeids, label=nodelab)
# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=0,branch.length='none')
p <- p %<+% data1+geom_tippoint(mapping=aes(color=Phylum2), 
                size=0.2,
                show.legend=T)+scale_color_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                                                                 "Deltaproteobacteria" = "darkolivegreen4", "Unidentified Proteobacteria" = "#003333","Actinobacteriota"="indianred2",
                                                                 "Firmicutes" ="tan1","Bacteroidota"="steelblue1","Cyanobacteria"="#557882","Planctomycetota"="#dc00f4",
                                                                 "Acidobacteriota"="#F5DC50","Patescibacteria"="#64508C","Bdellovibrionota"="#9BD2D2", "Deinococcota"="#ffaade",
                                                                 "Unidentified" = "black"))
p <- p + new_scale_fill() +
  geom_fruit(data=data2, geom=geom_tile,
             mapping=aes(y=ID, x=Branch_number, fill=sumPresence),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_fill_manual(values=c("0" = "#FFFFFF","1"="#7CB5D2","2"="#96CEB7","3"="#F8DFC8",
                             "4"="#BB94B7", "5"="#BDD74B","6"="#F68D51","7"="#FCE05E",
                             "8"="#18A68D","9"="#AC2324","10"="#104577"))+
  geom_fruit(data=data3, geom=geom_bar,
             mapping=aes(y=ID, x=Prevalence), fill="#4E362F",
             pwidth=0.38, 
             orientation="y", 
             stat="identity")+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p

dev.off()

### Fungi
tree <- phyloseq::phy_tree(fun.clean.ss.f)

## Data 1: Tip information (taxonomic information)
fun.list.2 <- fun.list
fun.list.2$Phylum[is.na(fun.list.2$Phylum)] <- "Unidentified"
fun.list.2$Class[is.na(fun.list.2$Class)] <- "Unidentified"
fun.list.2$Class2 <- fun.list.2$Class
fun.list.2$Class2[which(fun.list.2$Class == "Unidentified"& fun.list.2$Phylum == "Ascomycota")] <- "Unidentified Ascomycota"
fun.list.2$Class2[which(fun.list.2$Class == "Unidentified"& fun.list.2$Phylum == "Basidiomycota")] <- "Unidentified Basidiomycota"

fun.list.2[-c(3:10)]
data1 <- fun.list.2[-c(3:10)]
names(data1)[1]<-"ID"

#write.csv(data1,"data1_fungi.csv")
## Data 2: Specificity of ASVs (Heat map)
fun.clean.preval <- fun.clean.ss.f
fun.melt.preval<-psmelt(fun.clean.preval)
fun.melt.preval$Presence <- ifelse(fun.melt.preval$Abundance > 0, 1, 0)
fun.melt.preval <- fun.melt.preval%>% group_by(OTU, Branch_number) %>% summarise(sumPresence=sum(Presence))
names(fun.melt.preval)[1] <- "ID"
data2 <- fun.melt.preval
data2$sumPresence <- as.character(data2$sumPresence)
data2$sumPresence <- factor(data2$sumPresence, levels = c("0","1","2","3","4","5","6","7","8","9","10"))


## Data 3: Prevalence data of each ASVs in all 70 seeds (bar plot)
otu.G <- otu_table(fun.clean.ss.f)
df.otu.G <- data.frame(otu.G)
head(df.otu.G)

df.otu.G[df.otu.G>0] <- 1

df.otu.sum<-data.frame(rowSums(df.otu.G))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$ID <- rownames(df.otu.sum)

data3<- df.otu.sum[c(2,3)]
fun.prev.tab<-data3
names(fun.prev.tab)[2] <- "OTU" 
fun.prev.tab<-merge(fun.list, fun.prev.tab, by = c("OTU" = "OTU"))
write.csv(fun.prev.tab, "fun.prev.tab.csv")

tree <- phyloseq::phy_tree(fun.clean.ss.f)
unique(data1$Class2)
# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=0,branch.length='none')
p <- p %<+% data1+geom_tippoint(mapping=aes(color=Class2), 
                                size=0.2,
                                show.legend=T)+scale_color_manual(values = c("Dothideomycetes" = "#5195D1","Sordariomycetes"= "#1E63AF", "Eurotiomycetes"= "#6DA9DC",
                                                                             "Leotiomycetes"= "#11335F","Saccharomycetes"="#758C32","Pezizomycetes"="#ECE9B5",
                                                                             "Pezizomycotina_cls_Incertae_sedis"="#60848E","Unidentified Ascomycota" = "#13008A",
                                                                             "Ustilaginomycetes"="#ffcc33","Tremellomycetes"="#BE4146","Cystobasidiomycetes" = "#A871AE",
                                                                             "Microbotryomycetes" ="#DC9A9E","Agaricomycetes" = "#CC6C71","Malasseziomycetes" = "#99cc99",
                                                                             "Agaricostilbomycetes"="#003333",
                                                                             "Unidentified Basidiomycota"="#6A0014","Spizellomycetes" = "#D2C5BC",
                                                                             "Unidentified" ="black"))
p <- p + new_scale_fill() +
  geom_fruit(data=data2, geom=geom_tile,
             mapping=aes(y=ID, x=Branch_number, fill=sumPresence),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_fill_manual(values=c("0" = "#FFFFFF","1"="#7CB5D2","2"="#96CEB7","3"="#F8DFC8",
                             "4"="#BB94B7", "5"="#BDD74B","6"="#F68D51","7"="#FCE05E",
                             "8"="#18A68D","9"="#AC2324","10"="#104577"))+
  geom_fruit(data=data3, geom=geom_bar,
             mapping=aes(y=ID, x=Prevalence), fill="#4E362F",
             pwidth=0.38, 
             orientation="y", 
             stat="identity")+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p
