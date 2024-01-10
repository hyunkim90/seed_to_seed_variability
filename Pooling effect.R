### Pooled samples vs. individual seeds
### Dissimilarity
bac.time.clean.log
fun.time.clean.log

bac.clean.log
fun.clean.log

##Jaccard
### Pooled bacteria
sample_data(bac.time.clean.log)
bac.time.clean.log.18 <- subset_samples(bac.time.clean.log, Replication == "G_141")
bac.time.clean.log.18 <- phyloseq::filter_taxa(bac.time.clean.log.18, function(x) sum(x) != 0, TRUE)


b.otu.lognorm <- otu_table(bac.time.clean.log.18)

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
jaccard.dist.bac.melt.pooled <- subset(jaccard.dist.bac.melt, Jaccard_distance != 0)
jaccard.dist.bac.melt.pooled$Condition <- "Pooled"
mean(jaccard.dist.bac.melt.pooled$Jaccard_distance)
sd(jaccard.dist.bac.melt.pooled$Jaccard_distance)

### Individual bacteria
### Select 27 samples


b.otu.lognorm <- otu_table(bac.clean.log)

set.seed(1234)

simple.random.sampling <- sample(1:nrow(t(b.otu.lognorm)), 27)

b.otu.lognorm <- t(b.otu.lognorm)[simple.random.sampling, ]
dim(b.otu.lognorm)

rownames(b.otu.lognorm)


jaccard.dist.bac<-vegdist(b.otu.lognorm, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.bac <-as.matrix(jaccard.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

jaccard.dist.bac.lower<-get_lower_tri(jaccard.dist.bac)
jaccard.dist.bac.melt <- melt(as.matrix(jaccard.dist.bac.lower), na.rm = T)
head(jaccard.dist.bac.melt)
names(jaccard.dist.bac.melt)[3] <- "Jaccard_distance"
jaccard.dist.bac.melt.ind <- subset(jaccard.dist.bac.melt, Jaccard_distance != 0)
jaccard.dist.bac.melt.ind$Condition <- "Individual"


Pool.Ind.bac<-rbind(jaccard.dist.bac.melt.pooled,jaccard.dist.bac.melt.ind)

shapiro.test(Pool.Ind.bac$Jaccard_distance) #p-value = 4.614e-06

wilcox.test(subset(Pool.Ind.bac, Condition == "Pooled")$Jaccard_distance, subset(Pool.Ind.bac, Condition == "Individual")$Jaccard_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.bac, Condition == "Pooled")$Jaccard_distance) #0.5466935
mean(subset(Pool.Ind.bac, Condition == "Individual")$Jaccard_distance) #0.7986365

Pool.Ind.bac$Condition <- factor(Pool.Ind.bac$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.bac, aes(x=Condition, y= Jaccard_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


write.csv(Pool.Ind.bac, "230119_Jaccard_27 samples.csv")
#### All samples
### Individual bacteria
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
jaccard.dist.bac.melt.ind <- subset(jaccard.dist.bac.melt, Jaccard_distance != 0)
jaccard.dist.bac.melt.ind$Condition <- "Individual"
mean(jaccard.dist.bac.melt.ind$Jaccard_distance)
sd(jaccard.dist.bac.melt.ind$Jaccard_distance)

Pool.Ind.bac<-rbind(jaccard.dist.bac.melt.pooled,jaccard.dist.bac.melt.ind)

shapiro.test(Pool.Ind.bac$Jaccard_distance) #p-value < 2.2e-16

wilcox.test(subset(Pool.Ind.bac, Condition == "Pooled")$Jaccard_distance, subset(Pool.Ind.bac, Condition == "Individual")$Jaccard_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.bac, Condition == "Pooled")$Jaccard_distance) #0.5466935
mean(subset(Pool.Ind.bac, Condition == "Individual")$Jaccard_distance) #0.7977077


write.csv(Pool.Ind.bac,"20220606_pooling effect_jaccard distance.csv")


###Bray
### Pooled bacteria
sample_data(bac.time.clean.log)
bac.time.clean.log.18 <- subset_samples(bac.time.clean.log, Replication == "G_141")
bac.time.clean.log.18 <- phyloseq::filter_taxa(bac.time.clean.log.18, function(x) sum(x) != 0, TRUE)


b.otu.lognorm <- otu_table(bac.time.clean.log.18)

bray.dist.bac<-vegdist(t(b.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "Bray_distance"
bray.dist.bac.melt.pooled <- subset(bray.dist.bac.melt, Bray_distance != 0)
bray.dist.bac.melt.pooled$Condition <- "Pooled"

mean(bray.dist.bac.melt.pooled$Bray_distance)
sd(bray.dist.bac.melt.pooled$Bray_distance)


### Individual bacteria
### Select 27 samples


b.otu.lognorm <- otu_table(bac.clean.log)

set.seed(1234)

simple.random.sampling

b.otu.lognorm <- t(b.otu.lognorm)[simple.random.sampling, ]
dim(b.otu.lognorm)




bray.dist.bac<-vegdist(b.otu.lognorm, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "Bray_distance"
bray.dist.bac.melt.ind <- subset(bray.dist.bac.melt, Bray_distance != 0)
bray.dist.bac.melt.ind$Condition <- "Individual"


Pool.Ind.bac.bray<-rbind(bray.dist.bac.melt.pooled,bray.dist.bac.melt.ind)

shapiro.test(Pool.Ind.bac.bray$Bray_distance) #p-value = 1.828e-08

wilcox.test(subset(Pool.Ind.bac.bray, Condition == "Pooled")$Bray_distance, subset(Pool.Ind.bac.bray, Condition == "Individual")$Bray_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.bac.bray, Condition == "Pooled")$Bray_distance) # 0.3843621
mean(subset(Pool.Ind.bac.bray, Condition == "Individual")$Bray_distance) #0.6755359

Pool.Ind.bac.bray$Condition <- factor(Pool.Ind.bac.bray$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.bac.bray, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

write.csv(Pool.Ind.bac.bray, "230119_Bray_27 samples.csv")

#### All samples
### Individual bacteria
b.otu.lognorm <- otu_table(bac.clean.log)


bray.dist.bac<-vegdist(t(b.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "Bray_distance"
bray.dist.bac.melt.ind <- subset(bray.dist.bac.melt, Bray_distance != 0)
bray.dist.bac.melt.ind$Condition <- "Individual"
mean(bray.dist.bac.melt.ind$Bray_distance)
sd(bray.dist.bac.melt.ind$Bray_distance)
Pool.Ind.bac.bray<-rbind(bray.dist.bac.melt.pooled,bray.dist.bac.melt.ind)
write.csv(Pool.Ind.bac.bray,"20220606_pooling effect_bray distance.csv")

shapiro.test(Pool.Ind.bac.bray$Bray_distance) #p-value = 4.614e-06

wilcox.test(subset(Pool.Ind.bac.bray, Condition == "Pooled")$Bray_distance, subset(Pool.Ind.bac.bray, Condition == "Individual")$Bray_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.bac.bray, Condition == "Pooled")$Bray_distance) #0.3843621
mean(subset(Pool.Ind.bac.bray, Condition == "Individual")$Bray_distance) #0.6728886

Pool.Ind.bac.bray$Condition <- factor(Pool.Ind.bac.bray$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.bac.bray, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())





##Jaccard
### Pooled fungi
fun.time.clean.log.18 <- subset_samples(fun.time.clean.log, Replication == "G_141")
fun.time.clean.log.18 <- phyloseq::filter_taxa(fun.time.clean.log.18, function(x) sum(x) != 0, TRUE)


f.otu.lognorm <- otu_table(fun.time.clean.log.18)

jaccard.dist.bac<-vegdist(t(f.otu.lognorm), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.bac <-as.matrix(jaccard.dist.bac)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.bac)
jaccard.dist.fun.melt <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt)
names(jaccard.dist.fun.melt)[3] <- "Jaccard_distance"
jaccard.dist.fun.melt.pooled <- subset(jaccard.dist.fun.melt, Jaccard_distance != 0)
jaccard.dist.fun.melt.pooled$Condition <- "Pooled"

### Individual fungi
### Select 27 samples


f.otu.lognorm <- otu_table(fun.clean.log)

set.seed(1234)

simple.random.sampling <- sample(1:nrow(t(f.otu.lognorm)), 27)

f.otu.lognorm <- t(f.otu.lognorm)[simple.random.sampling, ]
dim(f.otu.lognorm)




jaccard.dist.fun<-vegdist(f.otu.lognorm, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.fun <-as.matrix(jaccard.dist.fun)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.fun)
jaccard.dist.fun.melt <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt)
names(jaccard.dist.fun.melt)[3] <- "Jaccard_distance"
jaccard.dist.fun.melt.ind <- subset(jaccard.dist.fun.melt, Jaccard_distance != 0)
jaccard.dist.fun.melt.ind$Condition <- "Individual"


Pool.Ind.fun<-rbind(jaccard.dist.fun.melt.pooled,jaccard.dist.fun.melt.ind)

shapiro.test(Pool.Ind.fun$Jaccard_distance) #p-value = 0.4789

wilcox.test(subset(Pool.Ind.fun, Condition == "Pooled")$Jaccard_distance, subset(Pool.Ind.fun, Condition == "Individual")$Jaccard_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.fun, Condition == "Pooled")$Jaccard_distance) # 0.7126094
mean(subset(Pool.Ind.fun, Condition == "Individual")$Jaccard_distance) #0.6362991

Pool.Ind.fun$Condition <- factor(Pool.Ind.fun$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.fun, aes(x=Condition, y= Jaccard_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



#### All samples
### Individual fungi

f.otu.lognorm <- otu_table(fun.clean.log)


jaccard.dist.fun<-vegdist(t(f.otu.lognorm), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.fun <-as.matrix(jaccard.dist.fun)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.fun)
jaccard.dist.fun.melt <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt)
names(jaccard.dist.fun.melt)[3] <- "Jaccard_distance"
jaccard.dist.fun.melt.ind <- subset(jaccard.dist.fun.melt, Jaccard_distance != 0)
jaccard.dist.fun.melt.ind$Condition <- "Individual"


Pool.Ind.fun<-rbind(jaccard.dist.fun.melt.pooled,jaccard.dist.fun.melt.ind)

shapiro.test(Pool.Ind.fun$Jaccard_distance) #p-value = 2.2e-05

wilcox.test(subset(Pool.Ind.fun, Condition == "Pooled")$Jaccard_distance, subset(Pool.Ind.fun, Condition == "Individual")$Jaccard_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.fun, Condition == "Pooled")$Jaccard_distance) #0.7126094
mean(subset(Pool.Ind.fun, Condition == "Individual")$Jaccard_distance) #0.6385054

Pool.Ind.fun$Condition <- factor(Pool.Ind.fun$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.fun, aes(x=Condition, y= Jaccard_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

###Bray
### Pooled fungi
fun.time.clean.log.18 <- subset_samples(fun.time.clean.log, Replication == "G_141")
fun.time.clean.log.18 <- phyloseq::filter_taxa(fun.time.clean.log.18, function(x) sum(x) != 0, TRUE)


f.otu.lognorm <- otu_table(fun.time.clean.log.18)

bray.dist.fun<-vegdist(t(f.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "Bray_distance"
bray.dist.fun.melt.pooled <- subset(bray.dist.fun.melt, Bray_distance != 0)
bray.dist.fun.melt.pooled$Condition <- "Pooled"

### Individual fungi
### Select 27 samples


f.otu.lognorm <- otu_table(fun.clean.log)

set.seed(1234)

simple.random.sampling <- sample(1:nrow(t(f.otu.lognorm)), 27)

f.otu.lognorm <- t(f.otu.lognorm)[simple.random.sampling, ]
dim(f.otu.lognorm)




bray.dist.fun<-vegdist(f.otu.lognorm, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "Bray_distance"
bray.dist.fun.melt.ind <- subset(bray.dist.fun.melt, Bray_distance != 0)
bray.dist.fun.melt.ind$Condition <- "Individual"


Pool.Ind.fun<-rbind(bray.dist.fun.melt.pooled,bray.dist.fun.melt.ind)

shapiro.test(Pool.Ind.fun$Bray_distance) #p-value = 0.2449

wilcox.test(subset(Pool.Ind.fun, Condition == "Pooled")$Bray_distance, subset(Pool.Ind.fun, Condition == "Individual")$Bray_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.fun, Condition == "Pooled")$Bray_distance) # 0.5566014
mean(subset(Pool.Ind.fun, Condition == "Individual")$Bray_distance) #0.4690041

Pool.Ind.fun$Condition <- factor(Pool.Ind.fun$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.fun, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



#### All samples
### Individual fungi
f.otu.lognorm <- otu_table(fun.clean.log)


bray.dist.fun<-vegdist(t(f.otu.lognorm), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "Bray_distance"
bray.dist.fun.melt.ind <- subset(bray.dist.fun.melt, Bray_distance != 0)
bray.dist.fun.melt.ind$Condition <- "Individual"


Pool.Ind.fun<-rbind(bray.dist.fun.melt.pooled,bray.dist.fun.melt.ind)

shapiro.test(Pool.Ind.fun$Bray_distance) #p-value = 4.614e-06

wilcox.test(subset(Pool.Ind.fun, Condition == "Pooled")$Bray_distance, subset(Pool.Ind.fun, Condition == "Individual")$Bray_distance) #p-value = 2.05e-06
mean(subset(Pool.Ind.fun, Condition == "Pooled")$Bray_distance) #0.5466935
mean(subset(Pool.Ind.fun, Condition == "Individual")$Bray_distance) #0.6728886

Pool.Ind.fun$Condition <- factor(Pool.Ind.fun$Condition, levels = c("Pooled","Individual"))

ggplot(Pool.Ind.fun, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())




## Comparison of alpha diversity
bac.time.clean.ss.18 <-subset_samples(bac.time.clean.ss.f, Year == "year_2018" & Age == "141days")
bac.time.clean.ss.18 <- phyloseq::filter_taxa(bac.time.clean.ss.18, function(x) sum(x) != 0, TRUE)

min(colSums(otu_table(bac.time.clean.ss.18))) #4255
min(colSums(otu_table(bac.clean.ss.f))) #3442

bac.time.rarefied.18<-rarefy_even_depth(bac.time.clean.ss.18 , sample.size = 3000,
                                     rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
bac.rarefied.18<-rarefy_even_depth(bac.clean.ss.f , sample.size = 3000,
                                        rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


##Alpha diversity

##Richness
## Chao1
Richtab.pooled <-microbiome::richness(bac.time.rarefied.18, index = "chao1")
head(Richtab.pooled)
mean(Richtab.pooled$chao1) #31.11111
sd(Richtab.pooled$chao1) # 8.750092
Richtab.pooled$Cond <- "Pooled"
Richtab.ind <-microbiome::richness(bac.rarefied.18, index = "chao1")
head(Richtab.ind)
mean(Richtab.ind$chao1) #26.93488
sd(Richtab.ind$chao1) #11.47618
Richtab.ind$Cond <- "Individual"

RichTab <- rbind(Richtab.pooled,Richtab.ind)
write.csv(RichTab,"220603_Chao1_Pooled vs Ind.csv")

##Observed ASVs
Observedtab.pooled <-microbiome::richness(bac.time.rarefied.18, index = "observed")
head(Observedtab.pooled)
mean(Observedtab.pooled$observed) #31.11111
sd(Observedtab.pooled$observed) # 8.750092
Observedtab.pooled$Cond <- "Pooled"
Observedtab.ind <-microbiome::richness(bac.rarefied.18, index = "observed")
head(Observedtab.ind)
mean(Observedtab.ind$observed) #26.93488
sd(Observedtab.ind$observed) #11.47618
Observedtab.ind$Cond <- "Individual"

ObservedTab <- rbind(Observedtab.pooled,Observedtab.ind)
write.csv(ObservedTab,"220606_Observed ASVs_Pooled vs Ind.csv")
##27 samples
random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

Observedtab.ind.sub<-subset(Observedtab.ind, rownames(Observedtab.ind) %in% random.sample)

Observedtab.2 <- rbind(Observedtab.pooled,Observedtab.ind.sub)
write.csv(Observedtab.2,"230119_Observed ASVs_Pooled vs Ind_27samples.csv")



##Diversity
Divertab.pooled <-microbiome::diversity(bac.time.rarefied.18, index = "shannon")
head(Divertab.pooled)
mean(Divertab.pooled$shannon) #3.227032
sd(Divertab.pooled$shannon) # 0.2778303
Divertab.pooled$Cond <- "Pooled"
Divertab.ind <-microbiome::diversity(bac.rarefied.18, index = "shannon")
head(Divertab.ind)
mean(Divertab.ind$shannon) # 2.223052
sd(Divertab.ind$shannon) #0.553643
Divertab.ind$Cond <- "Individual"

DiverTab <- rbind(Divertab.pooled,Divertab.ind)
write.csv(DiverTab,"220603_Shannon_Pooled vs Ind.csv")

##27 samples
random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

Divertab.ind.sub<-subset(Divertab.ind, rownames(Divertab.ind) %in% random.sample)

Divertab.2 <- rbind(Divertab.pooled,Divertab.ind.sub)
write.csv(Divertab.2,"230119_Shannon_Pooled vs Ind_27samples.csv")

##Evenness
Eventab.pooled <-microbiome::evenness(bac.time.rarefied.18, index = "simpson")
head(Eventab.pooled)
mean(Eventab.pooled$simpson) #0.7466689
sd(Eventab.pooled$simpson) # 0.07240954
Eventab.pooled$Cond <- "Pooled"

Eventab.ind <-microbiome::evenness(bac.rarefied.18, index = "simpson")
head(Eventab.ind)
mean(Eventab.ind$simpson) # 0.2890836
sd(Eventab.ind$simpson) #0.09609175
Eventab.ind$Cond <- "Individual"

EvenTab <- rbind(Eventab.pooled,Eventab.ind)
write.csv(EvenTab,"220603_Simpson_Pooled vs Ind.csv")

##27 samples
random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

Eventab.ind.sub<-subset(Eventab.ind, rownames(Eventab.ind) %in% random.sample)

EvenTab.2 <- rbind(Eventab.pooled,Eventab.ind.sub)
write.csv(EvenTab.2,"230119_Simpson_Pooled vs Ind_27samples.csv")


##Dominance
Domintab.pooled <-microbiome::dominance(bac.time.rarefied.18, index = "gini")
head(Domintab.pooled)
mean(Domintab.pooled$gini) #0.8447217
sd(Domintab.pooled$gini) # 0.04475872
Domintab.pooled$Cond <- "Pooled"

Domintab.ind <-microbiome::dominance(bac.rarefied.18, index = "gini")
head(Domintab.ind)
mean(Domintab.ind$gini) # 0.9841488
sd(Domintab.ind$gini) #0.008513602
Domintab.ind$Cond <- "Individual"

DominTab <- rbind(Domintab.pooled,Domintab.ind)
write.csv(DominTab,"220603_Gini_Pooled vs Ind.csv")

##27 samples
random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

DominTab.sub<-subset(Domintab.ind, rownames(Domintab.ind) %in% random.sample)

DominTab.2 <- rbind(Domintab.pooled,DominTab.sub)
write.csv(DominTab.2,"230119_Dominance_Pooled vs Ind_27samples.csv")

##Rarity
Raritab.pooled <-microbiome::rarity(bac.time.rarefied.18, index = "log_modulo_skewness")
head(Raritab.pooled)
mean(Raritab.pooled$log_modulo_skewness) #0.996321
sd(Raritab.pooled$log_modulo_skewness) # 0.1918643
Raritab.pooled$Cond <- "Pooled"

Raritab.ind <-microbiome::rarity(bac.rarefied.18, index = "log_modulo_skewness")
head(Raritab.ind)
mean(Raritab.ind$log_modulo_skewness) # 1.950529
sd(Raritab.ind$log_modulo_skewness) #0.119363
Raritab.ind$Cond <- "Individual"

RariTab <- rbind(Raritab.pooled,Raritab.ind)
write.csv(RariTab,"220603_Rarity_Pooled vs Ind.csv")


###27 samples
Raritab.ind <-microbiome::rarity(bac.rarefied.18, index = "log_modulo_skewness")
head(Raritab.ind)
mean(Raritab.ind$log_modulo_skewness) # 1.950529
sd(Raritab.ind$log_modulo_skewness) #0.119363
Raritab.ind$Cond <- "Individual"

random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

Raritab.ind.sub<-subset(Raritab.ind, rownames(Raritab.ind) %in% random.sample)

RariTab.2 <- rbind(Raritab.pooled,Raritab.ind.sub)
write.csv(RariTab.2,"230119_Rarity_Pooled vs Ind_27samples.csv")




### Plotting
Pool.Ind.bac$Condition <- factor(Pool.Ind.bac$Condition, levels = c("Pooled","Individual"))

p1<-ggplot(Pool.Ind.bac, aes(x=Condition, y= Jaccard_distance, fill = Condition))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

Pool.Ind.bac.bray$Condition <- factor(Pool.Ind.bac.bray$Condition, levels = c("Pooled","Individual"))

p2<-ggplot(Pool.Ind.bac.bray, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


RariTab$Cond <- factor(RariTab$Cond, levels = c("Pooled","Individual"))

p3<-ggplot(RariTab.2, aes(x=Cond, y= log_modulo_skewness, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


EvenTab$Cond <- factor(EvenTab$Cond, levels = c("Pooled","Individual"))

p4<-ggplot(EvenTab, aes(x=Cond, y= simpson, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


DominTab$Cond <- factor(DominTab$Cond, levels = c("Pooled","Individual"))

p5<-ggplot(DominTab, aes(x=Cond, y= gini, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

ObservedTab$Cond <- factor(ObservedTab$Cond, levels = c("Pooled","Individual"))

p6<-ggplot(ObservedTab, aes(x=Cond, y= observed, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red")+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

### Composite figure-Main

pdf("220607_Pooling effect_main figure_with mean value.pdf", width = 10, height =10)
ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,
                  labels = c("a","b","c","d","e","f"), ncol = 3, nrow=2)

dev.off()


##Supple Fig. 27 samples

Pool.Ind.bac$Condition <- factor(Pool.Ind.bac$Condition, levels = c("Pooled","Individual"))

p1<-ggplot(Pool.Ind.bac, aes(x=Condition, y= Jaccard_distance, fill = Condition))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

Pool.Ind.bac.bray$Condition <- factor(Pool.Ind.bac.bray$Condition, levels = c("Pooled","Individual"))

p2<-ggplot(Pool.Ind.bac.bray, aes(x=Condition, y= Bray_distance, fill = Condition))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


RariTab.2$Cond <- factor(RariTab.2$Cond, levels = c("Pooled","Individual"))

p3<-ggplot(RariTab.2, aes(x=Cond, y= log_modulo_skewness, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


EvenTab.2$Cond <- factor(EvenTab.2$Cond, levels = c("Pooled","Individual"))

p4<-ggplot(EvenTab.2, aes(x=Cond, y= simpson, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")


DominTab.2$Cond <- factor(DominTab.2$Cond, levels = c("Pooled","Individual"))

p5<-ggplot(DominTab.2, aes(x=Cond, y= gini, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red") +
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

Observedtab.2$Cond <- factor(Observedtab.2$Cond, levels = c("Pooled","Individual"))

p6<-ggplot(Observedtab.2, aes(x=Cond, y= observed, fill = Cond))+
  geom_boxplot()+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed", colour = "red")+
  ggsignif::stat_signif(comparisons = list(c("Pooled","Individual")), test = "wilcox.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")

### Composite figure-Main

pdf("230119_Pooling effect_supple figure_with mean value.pdf", width = 10, height =10)
ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,
                  labels = c("A","B","C","D","E","F"), ncol = 3, nrow=2)

dev.off()

