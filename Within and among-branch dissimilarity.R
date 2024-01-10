### Within-branch dissimilarity and Among-branch dissimilarity
##Correct metadata
b.meta <- read.table("./Bacteria/sample_metadata.tsv", sep = '\t', header =T)
f.meta <- read.table("./Fungi/sample_metadata.tsv", sep = '\t', header =T)

rownames(b.meta) <- b.meta$SampleID 
rownames(f.meta) <- f.meta$SampleID 

sample_data(bac.clean.ss) <- sample_data(b.meta)
sample_data(bac.clean.ss.f) <- sample_data(b.meta)
sample_data(bac.clean.nolog) <- sample_data(b.meta)
sample_data(bac.clean.log) <- sample_data(b.meta)
sample_data(bac.rarefied) <- sample_data(b.meta)

sample_data(fun.clean.ss) <- sample_data(f.meta)
sample_data(fun.clean.ss.f) <- sample_data(f.meta)
sample_data(fun.clean.nolog) <- sample_data(f.meta)
sample_data(fun.clean.log) <- sample_data(f.meta)
sample_data(fun.rarefied) <- sample_data(f.meta)

#### Jaccard dissimilarity - Bacteria
b.otu.lognorm <- otu_table(bac.clean.log)
jaccard.dist.bac<-vegdist(t(b.otu.lognorm), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.bac <-as.matrix(jaccard.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

jaccard.dist.bac.lower<-get_lower_tri(jaccard.dist.bac)
jaccard.dist.bac.melt <- melt(as.matrix(jaccard.dist.bac.lower), na.rm = T)
head(jaccard.dist.bac.melt)
names(jaccard.dist.bac.melt)[3] <- "Jaccard_distance"
jaccard.dist.bac.melt <- subset(jaccard.dist.bac.melt, Jaccard_distance != 0)

b.meta
#### Distance within each branch
Br2<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_2")]
Br3<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_3")]
Br4<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_4")]
Br5<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_5")]
Br6<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_6")]
Br7<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_7")]
Br8<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_8")]

## Within branch
jaccard.dist.bac.br2 <- subset(jaccard.dist.bac.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(jaccard.dist.bac.br2$Jaccard_distance)
jaccard.dist.bac.br2$Comparison <- "Branch_2"
mean(jaccard.dist.bac.br2$Jaccard_distance) #0.7905298

jaccard.dist.bac.br3 <- subset(jaccard.dist.bac.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(jaccard.dist.bac.br3$Jaccard_distance)
jaccard.dist.bac.br3$Comparison <- "Branch_3"
mean(jaccard.dist.bac.br3$Jaccard_distance) #0.7754927

jaccard.dist.bac.br4 <- subset(jaccard.dist.bac.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(jaccard.dist.bac.br4$Jaccard_distance)
jaccard.dist.bac.br4$Comparison <- "Branch_4"
mean(jaccard.dist.bac.br4$Jaccard_distance) #0.7748675

jaccard.dist.bac.br5 <- subset(jaccard.dist.bac.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(jaccard.dist.bac.br5$Jaccard_distance)
jaccard.dist.bac.br5$Comparison <- "Branch_5"
mean(jaccard.dist.bac.br5$Jaccard_distance) #0.7843951

jaccard.dist.bac.br6 <- subset(jaccard.dist.bac.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(jaccard.dist.bac.br6$Jaccard_distance)
jaccard.dist.bac.br6$Comparison <- "Branch_6"
mean(jaccard.dist.bac.br6$Jaccard_distance) #0.7457121

jaccard.dist.bac.br7 <- subset(jaccard.dist.bac.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(jaccard.dist.bac.br7$Jaccard_distance)
jaccard.dist.bac.br7$Comparison <- "Branch_7"
mean(jaccard.dist.bac.br7$Jaccard_distance) #0.7044251

jaccard.dist.bac.br8 <- subset(jaccard.dist.bac.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(jaccard.dist.bac.br8$Jaccard_distance)
jaccard.dist.bac.br8$Comparison <- "Branch_8"
mean(jaccard.dist.bac.br8$Jaccard_distance) #0.8543361

jaccard.dist.within.br <- rbind(jaccard.dist.bac.br2, jaccard.dist.bac.br3,jaccard.dist.bac.br4,
                                jaccard.dist.bac.br5,jaccard.dist.bac.br6,jaccard.dist.bac.br7,
                                jaccard.dist.bac.br8)
jaccard.dist.within.br$Group <- "Within_branch"

##Among branch
jaccard.dist.bac.br2.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
jaccard.dist.bac.br2.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
jaccard.dist.bac.br2.other <- rbind(jaccard.dist.bac.br2.other.1, jaccard.dist.bac.br2.other.2)
length(jaccard.dist.bac.br2.other$Jaccard_distance)
jaccard.dist.bac.br2.other$Comparison <- "Br2_Others"
mean(jaccard.dist.bac.br2.other$Jaccard_distance) #0.7944298

jaccard.dist.bac.br3.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
jaccard.dist.bac.br3.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
jaccard.dist.bac.br3.other <- rbind(jaccard.dist.bac.br3.other.1, jaccard.dist.bac.br3.other.2)
length(jaccard.dist.bac.br3.other$Jaccard_distance)
jaccard.dist.bac.br3.other$Comparison <- "Br3_Others"
mean(jaccard.dist.bac.br3.other$Jaccard_distance) #0.7950433

jaccard.dist.bac.br4.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
jaccard.dist.bac.br4.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
jaccard.dist.bac.br4.other <- rbind(jaccard.dist.bac.br4.other.1, jaccard.dist.bac.br4.other.2)
length(jaccard.dist.bac.br4.other$Jaccard_distance)
jaccard.dist.bac.br4.other$Comparison <- "Br4_Others"
mean(jaccard.dist.bac.br4.other$Jaccard_distance) #0.795574

jaccard.dist.bac.br5.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
jaccard.dist.bac.br5.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
jaccard.dist.bac.br5.other <- rbind(jaccard.dist.bac.br5.other.1, jaccard.dist.bac.br5.other.2)
length(jaccard.dist.bac.br5.other$Jaccard_distance)
jaccard.dist.bac.br5.other$Comparison <- "Br5_Others"
mean(jaccard.dist.bac.br5.other$Jaccard_distance) #0.8299794

jaccard.dist.bac.br6.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
jaccard.dist.bac.br6.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
jaccard.dist.bac.br6.other <- rbind(jaccard.dist.bac.br6.other.1, jaccard.dist.bac.br6.other.2)
length(jaccard.dist.bac.br6.other$Jaccard_distance)
jaccard.dist.bac.br6.other$Comparison <- "Br6_Others"
mean(jaccard.dist.bac.br6.other$Jaccard_distance) #0.8045498

jaccard.dist.bac.br7.other.1 <- subset(jaccard.dist.bac.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
jaccard.dist.bac.br7.other.2 <- subset(jaccard.dist.bac.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
jaccard.dist.bac.br7.other <- rbind(jaccard.dist.bac.br7.other.1, jaccard.dist.bac.br7.other.2)
length(jaccard.dist.bac.br7.other$Jaccard_distance)
jaccard.dist.bac.br7.other$Comparison <- "Br7_Others"
mean(jaccard.dist.bac.br7.other$Jaccard_distance) #0.7981205

jaccard.dist.among.br <- rbind(jaccard.dist.bac.br2.other, jaccard.dist.bac.br3.other,jaccard.dist.bac.br4.other,
                                jaccard.dist.bac.br5.other,jaccard.dist.bac.br6.other,jaccard.dist.bac.br7.other)
jaccard.dist.among.br$Group <- "Among_branch"


jaccard.dist.bac.br <- rbind(jaccard.dist.within.br,jaccard.dist.among.br)

### Statistical analysis
shapiro.test(jaccard.dist.bac.br$Jaccard_distance) #p-value = 6.218e-16
var.test(jaccard.dist.bac.br$Jaccard_distance~ jaccard.dist.bac.br$Group)

#wilcox.test(jaccard.dist.within.br$Jaccard_distance, jaccard.dist.among.br$Jaccard_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(jaccard.dist.within.br$Jaccard_distance) #0.7756797
mean(jaccard.dist.among.br$Jaccard_distance) #0.8010119

jaccard.dist.bac.br$Group <- factor(jaccard.dist.bac.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(jaccard.dist.bac.br, aes(x=Group, y= Jaccard_distance, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#### Jaccard dissimilarity - Fungi
f.otu.lognorm <- otu_table(fun.clean.log)
jaccard.dist.fun<-vegdist(t(f.otu.lognorm), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.fun <-as.matrix(jaccard.dist.fun)

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.fun)
jaccard.dist.fun.melt <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt)
names(jaccard.dist.fun.melt)[3] <- "Jaccard_distance"
jaccard.dist.fun.melt <- subset(jaccard.dist.fun.melt, Jaccard_distance != 0)

#### Distance within each branch
Br2<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_2")]
Br3<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_3")]
Br4<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_4")]
Br5<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_5")]
Br6<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_6")]
Br7<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_7")]
Br8<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_8")]

## Within branch
jaccard.dist.fun.br2 <- subset(jaccard.dist.fun.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(jaccard.dist.fun.br2$Jaccard_distance)
jaccard.dist.fun.br2$Comparison <- "Branch_2"
mean(jaccard.dist.fun.br2$Jaccard_distance) #0.6915578

jaccard.dist.fun.br3 <- subset(jaccard.dist.fun.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(jaccard.dist.fun.br3$Jaccard_distance)
jaccard.dist.fun.br3$Comparison <- "Branch_3"
mean(jaccard.dist.fun.br3$Jaccard_distance) #0.5773316

jaccard.dist.fun.br4 <- subset(jaccard.dist.fun.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(jaccard.dist.fun.br4$Jaccard_distance)
jaccard.dist.fun.br4$Comparison <- "Branch_4"
mean(jaccard.dist.fun.br4$Jaccard_distance) #0.5645032

jaccard.dist.fun.br5 <- subset(jaccard.dist.fun.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(jaccard.dist.fun.br5$Jaccard_distance)
jaccard.dist.fun.br5$Comparison <- "Branch_5"
mean(jaccard.dist.fun.br5$Jaccard_distance) #0.5947695

jaccard.dist.fun.br6 <- subset(jaccard.dist.fun.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(jaccard.dist.fun.br6$Jaccard_distance)
jaccard.dist.fun.br6$Comparison <- "Branch_6"
mean(jaccard.dist.fun.br6$Jaccard_distance) #0.6601957

jaccard.dist.fun.br7 <- subset(jaccard.dist.fun.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(jaccard.dist.fun.br7$Jaccard_distance)
jaccard.dist.fun.br7$Comparison <- "Branch_7"
mean(jaccard.dist.fun.br7$Jaccard_distance) # 0.6406678

jaccard.dist.fun.br8 <- subset(jaccard.dist.fun.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(jaccard.dist.fun.br8$Jaccard_distance)
jaccard.dist.fun.br8$Comparison <- "Branch_8"
mean(jaccard.dist.fun.br8$Jaccard_distance) #0.5475861

jaccard.dist.within.br.f <- rbind(jaccard.dist.fun.br2, jaccard.dist.fun.br3,jaccard.dist.fun.br4,
                                jaccard.dist.fun.br5,jaccard.dist.fun.br6,jaccard.dist.fun.br7,
                                jaccard.dist.fun.br8)
jaccard.dist.within.br.f$Group <- "Within_branch"

##Among branch
jaccard.dist.fun.br2.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
jaccard.dist.fun.br2.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
jaccard.dist.fun.br2.other <- rbind(jaccard.dist.fun.br2.other.1, jaccard.dist.fun.br2.other.2)
length(jaccard.dist.fun.br2.other$Jaccard_distance)
jaccard.dist.fun.br2.other$Comparison <- "Br2_Others"
mean(jaccard.dist.fun.br2.other$Jaccard_distance) #0.6776587

jaccard.dist.fun.br3.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
jaccard.dist.fun.br3.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
jaccard.dist.fun.br3.other <- rbind(jaccard.dist.fun.br3.other.1, jaccard.dist.fun.br3.other.2)
length(jaccard.dist.fun.br3.other$Jaccard_distance)
jaccard.dist.fun.br3.other$Comparison <- "Br3_Others"
mean(jaccard.dist.fun.br3.other$Jaccard_distance) #0.6188844

jaccard.dist.fun.br4.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
jaccard.dist.fun.br4.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
jaccard.dist.fun.br4.other <- rbind(jaccard.dist.fun.br4.other.1, jaccard.dist.fun.br4.other.2)
length(jaccard.dist.fun.br4.other$Jaccard_distance)
jaccard.dist.fun.br4.other$Comparison <- "Br4_Others"
mean(jaccard.dist.fun.br4.other$Jaccard_distance) #0.6222914

jaccard.dist.fun.br5.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
jaccard.dist.fun.br5.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
jaccard.dist.fun.br5.other <- rbind(jaccard.dist.fun.br5.other.1, jaccard.dist.fun.br5.other.2)
length(jaccard.dist.fun.br5.other$Jaccard_distance)
jaccard.dist.fun.br5.other$Comparison <- "Br5_Others"
mean(jaccard.dist.fun.br5.other$Jaccard_distance) #0.637669

jaccard.dist.fun.br6.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
jaccard.dist.fun.br6.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
jaccard.dist.fun.br6.other <- rbind(jaccard.dist.fun.br6.other.1, jaccard.dist.fun.br6.other.2)
length(jaccard.dist.fun.br6.other$Jaccard_distance)
jaccard.dist.fun.br6.other$Comparison <- "Br6_Others"
mean(jaccard.dist.fun.br6.other$Jaccard_distance) #0.6567838

jaccard.dist.fun.br7.other.1 <- subset(jaccard.dist.fun.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
jaccard.dist.fun.br7.other.2 <- subset(jaccard.dist.fun.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
jaccard.dist.fun.br7.other <- rbind(jaccard.dist.fun.br7.other.1, jaccard.dist.fun.br7.other.2)
length(jaccard.dist.fun.br7.other$Jaccard_distance)
jaccard.dist.fun.br7.other$Comparison <- "Br7_Others"
mean(jaccard.dist.fun.br7.other$Jaccard_distance) #0.6193156

jaccard.dist.among.br.f <- rbind(jaccard.dist.fun.br2.other, jaccard.dist.fun.br3.other,jaccard.dist.fun.br4.other,
                               jaccard.dist.fun.br5.other,jaccard.dist.fun.br6.other,jaccard.dist.fun.br7.other)
jaccard.dist.among.br.f$Group <- "Among_branch"


jaccard.dist.fun.br <- rbind(jaccard.dist.within.br.f,jaccard.dist.among.br.f)

### Statistical analysis
shapiro.test(jaccard.dist.fun.br$Jaccard_distance) #p-value = 6.539e-10
var.test(jaccard.dist.fun.br$Jaccard_distance~ jaccard.dist.fun.br$Group)

#wilcox.test(jaccard.dist.within.br.f$Jaccard_distance, jaccard.dist.among.br.f$Jaccard_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(jaccard.dist.within.br.f$Jaccard_distance) #0.6109445
mean(jaccard.dist.among.br.f$Jaccard_distance) #0.6426395

jaccard.dist.fun.br$Group <- factor(jaccard.dist.fun.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(jaccard.dist.fun.br, aes(x=Group, y= Jaccard_distance, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Bray-Curtis
b.otu.lognorm <- otu_table(bac.clean.log)
bray.dist.bac<-vegdist(t(b.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)


bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "BC_distance"
bray.dist.bac.melt <- subset(bray.dist.bac.melt, BC_distance != 0)

#### Distance within each branch
Br2<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_2")]
Br3<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_3")]
Br4<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_4")]
Br5<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_5")]
Br6<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_6")]
Br7<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_7")]
Br8<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_8")]

## Within branch
bray.dist.bac.br2 <- subset(bray.dist.bac.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(bray.dist.bac.br2$BC_distance)
bray.dist.bac.br2$Comparison <- "Branch_2"
mean(bray.dist.bac.br2$BC_distance) #0.6629088

bray.dist.bac.br3 <- subset(bray.dist.bac.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(bray.dist.bac.br3$BC_distance)
bray.dist.bac.br3$Comparison <- "Branch_3"
mean(bray.dist.bac.br3$BC_distance) #0.6403801

bray.dist.bac.br4 <- subset(bray.dist.bac.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(bray.dist.bac.br4$BC_distance)
bray.dist.bac.br4$Comparison <- "Branch_4"
mean(bray.dist.bac.br4$BC_distance) # 0.6385357

bray.dist.bac.br5 <- subset(bray.dist.bac.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(bray.dist.bac.br5$BC_distance)
bray.dist.bac.br5$Comparison <- "Branch_5"
mean(bray.dist.bac.br5$BC_distance) #0.6498912

bray.dist.bac.br6 <- subset(bray.dist.bac.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(bray.dist.bac.br6$BC_distance)
bray.dist.bac.br6$Comparison <- "Branch_6"
mean(bray.dist.bac.br6$BC_distance) #0.6012618

bray.dist.bac.br7 <- subset(bray.dist.bac.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(bray.dist.bac.br7$BC_distance)
bray.dist.bac.br7$Comparison <- "Branch_7"
mean(bray.dist.bac.br7$BC_distance) #0.5497782

bray.dist.bac.br8 <- subset(bray.dist.bac.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(bray.dist.bac.br8$BC_distance)
bray.dist.bac.br8$Comparison <- "Branch_8"
mean(bray.dist.bac.br8$BC_distance) #0.7647707

bray.dist.within.br <- rbind(bray.dist.bac.br2, bray.dist.bac.br3,bray.dist.bac.br4,
                                bray.dist.bac.br5,bray.dist.bac.br6,bray.dist.bac.br7,
                                bray.dist.bac.br8)
bray.dist.within.br$Group <- "Within_branch"

##Among branch
bray.dist.bac.br2.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
bray.dist.bac.br2.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
bray.dist.bac.br2.other <- rbind(bray.dist.bac.br2.other.1, bray.dist.bac.br2.other.2)
length(bray.dist.bac.br2.other$BC_distance)
bray.dist.bac.br2.other$Comparison <- "Br2_Others"
mean(bray.dist.bac.br2.other$BC_distance) #0.6685275

bray.dist.bac.br3.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
bray.dist.bac.br3.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
bray.dist.bac.br3.other <- rbind(bray.dist.bac.br3.other.1, bray.dist.bac.br3.other.2)
length(bray.dist.bac.br3.other$BC_distance)
bray.dist.bac.br3.other$Comparison <- "Br3_Others"
mean(bray.dist.bac.br3.other$BC_distance) #0.6669901

bray.dist.bac.br4.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
bray.dist.bac.br4.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
bray.dist.bac.br4.other <- rbind(bray.dist.bac.br4.other.1, bray.dist.bac.br4.other.2)
length(bray.dist.bac.br4.other$BC_distance)
bray.dist.bac.br4.other$Comparison <- "Br4_Others"
mean(bray.dist.bac.br4.other$BC_distance) #0.6695561

bray.dist.bac.br5.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
bray.dist.bac.br5.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
bray.dist.bac.br5.other <- rbind(bray.dist.bac.br5.other.1, bray.dist.bac.br5.other.2)
length(bray.dist.bac.br5.other$BC_distance)
bray.dist.bac.br5.other$Comparison <- "Br5_Others"
mean(bray.dist.bac.br5.other$BC_distance) #0.716338

bray.dist.bac.br6.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
bray.dist.bac.br6.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
bray.dist.bac.br6.other <- rbind(bray.dist.bac.br6.other.1, bray.dist.bac.br6.other.2)
length(bray.dist.bac.br6.other$BC_distance)
bray.dist.bac.br6.other$Comparison <- "Br6_Others"
mean(bray.dist.bac.br6.other$BC_distance) #0.6837204

bray.dist.bac.br7.other.1 <- subset(bray.dist.bac.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
bray.dist.bac.br7.other.2 <- subset(bray.dist.bac.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
bray.dist.bac.br7.other <- rbind(bray.dist.bac.br7.other.1, bray.dist.bac.br7.other.2)
length(bray.dist.bac.br7.other$BC_distance)
bray.dist.bac.br7.other$Comparison <- "Br7_Others"
mean(bray.dist.bac.br7.other$BC_distance) #0.681077

bray.dist.among.br <- rbind(bray.dist.bac.br2.other, bray.dist.bac.br3.other,bray.dist.bac.br4.other,
                               bray.dist.bac.br5.other,bray.dist.bac.br6.other,bray.dist.bac.br7.other)
bray.dist.among.br$Group <- "Among_branch"


bray.dist.bac.br <- rbind(bray.dist.within.br,bray.dist.among.br)

### Statistical analysis
shapiro.test(bray.dist.bac.br$BC_distance) #p-value < 2.2e-16
var.test(bray.dist.bac.br$BC_distance~ bray.dist.bac.br$Group)

#wilcox.test(bray.dist.within.br$BC_distance, bray.dist.among.br$BC_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(bray.dist.within.br$BC_distance) #0.6439324
mean(bray.dist.among.br$BC_distance) #0.677232

bray.dist.bac.br$Group <- factor(bray.dist.bac.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(bray.dist.bac.br, aes(x=Group, y= BC_distance, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#### Jaccard dissimilarity - Fungi
f.otu.lognorm <- otu_table(fun.clean.log)
bray.dist.fun<-vegdist(t(f.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "BC_distance"
bray.dist.fun.melt <- subset(bray.dist.fun.melt, BC_distance != 0)

#### Distance within each branch
Br2<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_2")]
Br3<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_3")]
Br4<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_4")]
Br5<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_5")]
Br6<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_6")]
Br7<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_7")]
Br8<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_8")]

## Within branch
bray.dist.fun.br2 <- subset(bray.dist.fun.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(bray.dist.fun.br2$BC_distance)
bray.dist.fun.br2$Comparison <- "Branch_2"
mean(bray.dist.fun.br2$BC_distance) #0.5300283

bray.dist.fun.br3 <- subset(bray.dist.fun.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(bray.dist.fun.br3$BC_distance)
bray.dist.fun.br3$Comparison <- "Branch_3"
mean(bray.dist.fun.br3$BC_distance) #0.4102179

bray.dist.fun.br4 <- subset(bray.dist.fun.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(bray.dist.fun.br4$BC_distance)
bray.dist.fun.br4$Comparison <- "Branch_4"
mean(bray.dist.fun.br4$BC_distance) #0.3958264

bray.dist.fun.br5 <- subset(bray.dist.fun.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(bray.dist.fun.br5$BC_distance)
bray.dist.fun.br5$Comparison <- "Branch_5"
mean(bray.dist.fun.br5$BC_distance) #0.4254783

bray.dist.fun.br6 <- subset(bray.dist.fun.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(bray.dist.fun.br6$BC_distance)
bray.dist.fun.br6$Comparison <- "Branch_6"
mean(bray.dist.fun.br6$BC_distance) #0.4954319

bray.dist.fun.br7 <- subset(bray.dist.fun.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(bray.dist.fun.br7$BC_distance)
bray.dist.fun.br7$Comparison <- "Branch_7"
mean(bray.dist.fun.br7$BC_distance) # 0.4737623

bray.dist.fun.br8 <- subset(bray.dist.fun.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(bray.dist.fun.br8$BC_distance)
bray.dist.fun.br8$Comparison <- "Branch_8"
mean(bray.dist.fun.br8$BC_distance) #0.3784891

bray.dist.within.br.f <- rbind(bray.dist.fun.br2, bray.dist.fun.br3,bray.dist.fun.br4,
                                  bray.dist.fun.br5,bray.dist.fun.br6,bray.dist.fun.br7,
                                  bray.dist.fun.br8)
bray.dist.within.br.f$Group <- "Within_branch"

##Among branch
bray.dist.fun.br2.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
bray.dist.fun.br2.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
bray.dist.fun.br2.other <- rbind(bray.dist.fun.br2.other.1, bray.dist.fun.br2.other.2)
length(bray.dist.fun.br2.other$BC_distance)
bray.dist.fun.br2.other$Comparison <- "Br2_Others"
mean(bray.dist.fun.br2.other$BC_distance) #0.514143

bray.dist.fun.br3.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
bray.dist.fun.br3.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
bray.dist.fun.br3.other <- rbind(bray.dist.fun.br3.other.1, bray.dist.fun.br3.other.2)
length(bray.dist.fun.br3.other$BC_distance)
bray.dist.fun.br3.other$Comparison <- "Br3_Others"
mean(bray.dist.fun.br3.other$BC_distance) #0.4509968

bray.dist.fun.br4.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
bray.dist.fun.br4.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
bray.dist.fun.br4.other <- rbind(bray.dist.fun.br4.other.1, bray.dist.fun.br4.other.2)
length(bray.dist.fun.br4.other$BC_distance)
bray.dist.fun.br4.other$Comparison <- "Br4_Others"
mean(bray.dist.fun.br4.other$BC_distance) #0.4538591

bray.dist.fun.br5.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
bray.dist.fun.br5.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
bray.dist.fun.br5.other <- rbind(bray.dist.fun.br5.other.1, bray.dist.fun.br5.other.2)
length(bray.dist.fun.br5.other$BC_distance)
bray.dist.fun.br5.other$Comparison <- "Br5_Others"
mean(bray.dist.fun.br5.other$BC_distance) #0.4705515

bray.dist.fun.br6.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
bray.dist.fun.br6.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
bray.dist.fun.br6.other <- rbind(bray.dist.fun.br6.other.1, bray.dist.fun.br6.other.2)
length(bray.dist.fun.br6.other$BC_distance)
bray.dist.fun.br6.other$Comparison <- "Br6_Others"
mean(bray.dist.fun.br6.other$BC_distance) #0.4910826

bray.dist.fun.br7.other.1 <- subset(bray.dist.fun.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
bray.dist.fun.br7.other.2 <- subset(bray.dist.fun.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
bray.dist.fun.br7.other <- rbind(bray.dist.fun.br7.other.1, bray.dist.fun.br7.other.2)
length(bray.dist.fun.br7.other$BC_distance)
bray.dist.fun.br7.other$Comparison <- "Br7_Others"
mean(bray.dist.fun.br7.other$BC_distance) #0.4512072

bray.dist.among.br.f <- rbind(bray.dist.fun.br2.other, bray.dist.fun.br3.other,bray.dist.fun.br4.other,
                                 bray.dist.fun.br5.other,bray.dist.fun.br6.other,bray.dist.fun.br7.other)
bray.dist.among.br.f$Group <- "Among_branch"


bray.dist.fun.br <- rbind(bray.dist.within.br.f,bray.dist.among.br.f)

### Statistical analysis
shapiro.test(bray.dist.fun.br$BC_distance) # p-value = 0.01291
var.test(bray.dist.fun.br$BC_distance~ bray.dist.fun.br$Group)
wilcox.test(bray.dist.within.br.f$BC_distance, bray.dist.among.br.f$BC_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(bray.dist.within.br.f$BC_distance) # 0.4441763
mean(bray.dist.among.br.f$BC_distance) #0.476205

bray.dist.fun.br$Group <- factor(bray.dist.fun.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(bray.dist.fun.br, aes(x=Group, y= BC_distance, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



