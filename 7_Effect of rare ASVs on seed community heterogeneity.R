### Effect of rare ASVs on seed community heterogeneity
## Analysis 1
#Prevalence data of each ASVs in all 70 seeds
otu.G <- otu_table(bac.clean.ss.f)
df.otu.G <- data.frame(otu.G)
head(df.otu.G)

df.otu.G[df.otu.G>0] <- 1

df.otu.sum<-data.frame(rowSums(df.otu.G))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$OTU <- rownames(df.otu.sum)

bac.prev.tab<- df.otu.sum[c(2,3)]

otu.G <- otu_table(fun.clean.ss.f)
df.otu.G <- data.frame(otu.G)
head(df.otu.G)

df.otu.G[df.otu.G>0] <- 1

df.otu.sum<-data.frame(rowSums(df.otu.G))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$OTU <- rownames(df.otu.sum)

fun.prev.tab<- df.otu.sum[c(2,3)]


#### Rare ASVs
rare.bac<-bac.prev.tab$OTU[which(bac.prev.tab$Prevalence <= 0.2)] #584/619
rare.fun<-fun.prev.tab$OTU[which(fun.prev.tab$Prevalence <= 0.2)] #339/387

### unweighted Unifrac distance (all ASVs)
## Bacteria
## All ASVs
unifrac.dist.bac<-UniFrac(bac.clean.log, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.bac <-as.matrix(unifrac.dist.bac)

unifrac.dist.bac.melt <- melt(as.matrix(unifrac.dist.bac), na.rm = T)
names(unifrac.dist.bac.melt)[3] <- "All_ASVs"
unifrac.dist.bac.melt <- subset(unifrac.dist.bac.melt, All_ASVs != 0)
head(unifrac.dist.bac.melt)
mean(unifrac.dist.bac.melt$All_ASVs) #0.5995192

## without rare ASVs
bac.log.wo.rare <- pop_taxa(bac.clean.log, rare.bac)

(filt.sample <- sample_sums(bac.log.wo.rare) > 0)
sum(sample_sums(bac.log.wo.rare) <= 0)  ## 2 sample discarded
bac.log.wo.rare <- prune_samples(filt.sample, bac.log.wo.rare)
bac.log.wo.rare 

unifrac.dist.bac<-UniFrac(bac.log.wo.rare, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.bac <-as.matrix(unifrac.dist.bac)

unifrac.dist.bac.melt.wo.rare <- melt(as.matrix(unifrac.dist.bac), na.rm = T)
names(unifrac.dist.bac.melt.wo.rare)[3] <- "Wo_rare"
unifrac.dist.bac.melt.wo.rare <- subset(unifrac.dist.bac.melt.wo.rare, Wo_rare != 0)
head(unifrac.dist.bac.melt.wo.rare)
mean(unifrac.dist.bac.melt.wo.rare$Wo_rare) #0.4743252

## Combine two data frames
bac.dist.uw.unifrac <- merge(unifrac.dist.bac.melt, unifrac.dist.bac.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
bac.dist.uw.unifrac <- melt(bac.dist.uw.unifrac)
head(bac.dist.uw.unifrac)

## Normality test (Anderson-Darling normality test)
ad.test(bac.dist.uw.unifrac$value)

ggplot(bac.dist.uw.unifrac, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## Fungi
## All ASVs
unifrac.dist.fun<-UniFrac(fun.clean.log, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.fun <-as.matrix(unifrac.dist.fun)

unifrac.dist.fun.melt <- melt(as.matrix(unifrac.dist.fun), na.rm = T)
names(unifrac.dist.fun.melt)[3] <- "All_ASVs"
unifrac.dist.fun.melt <- subset(unifrac.dist.fun.melt, All_ASVs != 0)
head(unifrac.dist.fun.melt)
mean(unifrac.dist.fun.melt$All_ASVs) #0.6903096

## without rare ASVs
fun.log.wo.rare <- pop_taxa(fun.clean.log, rare.fun)

(filt.sample <- sample_sums(fun.log.wo.rare) > 0)
sum(sample_sums(fun.log.wo.rare) <= 0)  ## 0 sample discarded
fun.log.wo.rare <- prune_samples(filt.sample, fun.log.wo.rare)
fun.log.wo.rare 

unifrac.dist.fun<-UniFrac(fun.log.wo.rare, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.fun <-as.matrix(unifrac.dist.fun)

unifrac.dist.fun.melt.wo.rare <- melt(as.matrix(unifrac.dist.fun), na.rm = T)
names(unifrac.dist.fun.melt.wo.rare)[3] <- "Wo_rare"
unifrac.dist.fun.melt.wo.rare <- subset(unifrac.dist.fun.melt.wo.rare, Wo_rare != 0)
head(unifrac.dist.fun.melt.wo.rare)
mean(unifrac.dist.fun.melt.wo.rare$Wo_rare) #0.3142531

## Combine two data frames
fun.dist.uw.unifrac <- merge(unifrac.dist.fun.melt, unifrac.dist.fun.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
fun.dist.uw.unifrac <- melt(fun.dist.uw.unifrac)
head(fun.dist.uw.unifrac)

## Normality test (Anderson-Darling normality test)
ad.test(fun.dist.uw.unifrac$value)

ggplot(fun.dist.uw.unifrac, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Jaccard distance
## Bacteria
## All ASVs
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
names(jaccard.dist.bac.melt)[3] <- "All_ASVs"
jaccard.dist.bac.melt <- subset(jaccard.dist.bac.melt, All_ASVs != 0)


## wo Rare
b.otu.lognorm.wo.rare <- otu_table(bac.log.wo.rare)
b.otu.lognorm.wo.rare <- data.frame(b.otu.lognorm.wo.rare)
colSums(b.otu.lognorm.wo.rare)
jaccard.dist.bac<-vegdist(t(b.otu.lognorm.wo.rare), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
jaccard.dist.bac <-as.matrix(jaccard.dist.bac)

jaccard.dist.bac.lower<-get_lower_tri(jaccard.dist.bac)
jaccard.dist.bac.melt.wo.rare <- melt(as.matrix(jaccard.dist.bac.lower), na.rm = T)
head(jaccard.dist.bac.melt.wo.rare)
names(jaccard.dist.bac.melt.wo.rare)[3] <- "Wo_rare"
jaccard.dist.bac.melt.wo.rare <- subset(jaccard.dist.bac.melt.wo.rare, Wo_rare != 0)

## Combine two data frames
bac.dist.jaccard <- merge(jaccard.dist.bac.melt, jaccard.dist.bac.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
bac.dist.jaccard <- melt(bac.dist.jaccard)
head(bac.dist.jaccard)

## Normality test (Anderson-Darling normality test)
ad.test(bac.dist.jaccard$value)

ggplot(bac.dist.jaccard, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## Fungi
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
names(jaccard.dist.fun.melt)[3] <- "All_ASVs"
jaccard.dist.fun.melt <- subset(jaccard.dist.fun.melt, All_ASVs != 0)


## wo Rare
f.otu.lognorm.wo.rare <- otu_table(fun.log.wo.rare)
f.otu.lognorm.wo.rare <- data.frame(f.otu.lognorm.wo.rare)
colSums(f.otu.lognorm.wo.rare)
jaccard.dist.fun<-vegdist(t(f.otu.lognorm.wo.rare), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
jaccard.dist.fun <-as.matrix(jaccard.dist.fun)

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.fun)
jaccard.dist.fun.melt.wo.rare <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt.wo.rare)
names(jaccard.dist.fun.melt.wo.rare)[3] <- "Wo_rare"
jaccard.dist.fun.melt.wo.rare <- subset(jaccard.dist.fun.melt.wo.rare, Wo_rare != 0)

## Combine two data frames
fun.dist.jaccard <- merge(jaccard.dist.fun.melt, jaccard.dist.fun.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
fun.dist.jaccard <- melt(fun.dist.jaccard)
head(fun.dist.jaccard)

## Normality test (Anderson-Darling normality test)
ad.test(fun.dist.jaccard$value)

ggplot(fun.dist.jaccard, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



##Analysis 2: non-significant differences between among and within branch by removing rare ASVs?
#### Distance within each branch
Br2<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_2")]
Br3<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_3")]
Br4<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_4")]
Br5<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_5")]
Br6<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_6")]
Br7<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_7")]
Br8<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_8")]

## Within branch
jaccard.dist.bac.rare2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(jaccard.dist.bac.rare2$Wo_rare)
jaccard.dist.bac.rare2$Comparison <- "Branch_2"
mean(jaccard.dist.bac.rare2$Wo_rare) #0.6835808

jaccard.dist.bac.rare3 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(jaccard.dist.bac.rare3$Wo_rare)
jaccard.dist.bac.rare3$Comparison <- "Branch_3"
mean(jaccard.dist.bac.rare3$Wo_rare) #0.6174074

jaccard.dist.bac.rare4 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(jaccard.dist.bac.rare4$Wo_rare)
jaccard.dist.bac.rare4$Comparison <- "Branch_4"
mean(jaccard.dist.bac.rare4$Wo_rare) #0.6371094

jaccard.dist.bac.rare5 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(jaccard.dist.bac.rare5$Wo_rare)
jaccard.dist.bac.rare5$Comparison <- "Branch_5"
mean(jaccard.dist.bac.rare5$Wo_rare) #0.6169486

jaccard.dist.bac.rare6 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(jaccard.dist.bac.rare6$Wo_rare)
jaccard.dist.bac.rare6$Comparison <- "Branch_6"
mean(jaccard.dist.bac.rare6$Wo_rare) #0.59317

jaccard.dist.bac.rare7 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(jaccard.dist.bac.rare7$Wo_rare)
jaccard.dist.bac.rare7$Comparison <- "Branch_7"
mean(jaccard.dist.bac.rare7$Wo_rare) # 0.5444338

jaccard.dist.bac.rare8 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(jaccard.dist.bac.rare8$Wo_rare)
jaccard.dist.bac.rare8$Comparison <- "Branch_8"
mean(jaccard.dist.bac.rare8$Wo_rare) #0.6330365

jaccard.dist.within.rare <- rbind(jaccard.dist.bac.rare2, jaccard.dist.bac.rare3,jaccard.dist.bac.rare4,
                                jaccard.dist.bac.rare5,jaccard.dist.bac.rare6,jaccard.dist.bac.rare7,
                                jaccard.dist.bac.rare8)
jaccard.dist.within.rare$Group <- "Within_branch"

##Among branch
jaccard.dist.bac.rare2.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
jaccard.dist.bac.rare2.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
jaccard.dist.bac.rare2.other <- rbind(jaccard.dist.bac.rare2.other.1, jaccard.dist.bac.rare2.other.2)
length(jaccard.dist.bac.rare2.other$Wo_rare)
jaccard.dist.bac.rare2.other$Comparison <- "Br2_Others"
mean(jaccard.dist.bac.rare2.other$Wo_rare) #0.7944298

jaccard.dist.bac.rare3.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
jaccard.dist.bac.rare3.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
jaccard.dist.bac.rare3.other <- rbind(jaccard.dist.bac.rare3.other.1, jaccard.dist.bac.rare3.other.2)
length(jaccard.dist.bac.rare3.other$Wo_rare)
jaccard.dist.bac.rare3.other$Comparison <- "Br3_Others"
mean(jaccard.dist.bac.rare3.other$Wo_rare) #0.7950433

jaccard.dist.bac.rare4.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
jaccard.dist.bac.rare4.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
jaccard.dist.bac.rare4.other <- rbind(jaccard.dist.bac.rare4.other.1, jaccard.dist.bac.rare4.other.2)
length(jaccard.dist.bac.rare4.other$Wo_rare)
jaccard.dist.bac.rare4.other$Comparison <- "Br4_Others"
mean(jaccard.dist.bac.rare4.other$Wo_rare) #0.795574

jaccard.dist.bac.rare5.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
jaccard.dist.bac.rare5.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
jaccard.dist.bac.rare5.other <- rbind(jaccard.dist.bac.rare5.other.1, jaccard.dist.bac.rare5.other.2)
length(jaccard.dist.bac.rare5.other$Wo_rare)
jaccard.dist.bac.rare5.other$Comparison <- "Br5_Others"
mean(jaccard.dist.bac.rare5.other$Wo_rare) #0.8299794

jaccard.dist.bac.rare6.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
jaccard.dist.bac.rare6.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
jaccard.dist.bac.rare6.other <- rbind(jaccard.dist.bac.rare6.other.1, jaccard.dist.bac.rare6.other.2)
length(jaccard.dist.bac.rare6.other$Wo_rare)
jaccard.dist.bac.rare6.other$Comparison <- "Br6_Others"
mean(jaccard.dist.bac.rare6.other$Wo_rare) #0.8045498

jaccard.dist.bac.rare7.other.1 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
jaccard.dist.bac.rare7.other.2 <- subset(jaccard.dist.bac.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
jaccard.dist.bac.rare7.other <- rbind(jaccard.dist.bac.rare7.other.1, jaccard.dist.bac.rare7.other.2)
length(jaccard.dist.bac.rare7.other$Wo_rare)
jaccard.dist.bac.rare7.other$Comparison <- "Br7_Others"
mean(jaccard.dist.bac.rare7.other$Wo_rare) #0.7981205

jaccard.dist.among.rare <- rbind(jaccard.dist.bac.rare2.other, jaccard.dist.bac.rare3.other,jaccard.dist.bac.rare4.other,
                               jaccard.dist.bac.rare5.other,jaccard.dist.bac.rare6.other,jaccard.dist.bac.rare7.other)
jaccard.dist.among.rare$Group <- "Among_branch"


jaccard.dist.bac.rare <- rbind(jaccard.dist.within.rare,jaccard.dist.among.rare)

### Statistical analysis
shapiro.test(jaccard.dist.bac.rare$Wo_rare) #p-value = 6.218e-16
var.test(jaccard.dist.bac.rare$Wo_rare~ jaccard.dist.bac.rare$Group)

#wilcox.test(jaccard.dist.within.rare$Wo_rare, jaccard.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(jaccard.dist.within.rare$Wo_rare) #0.6170949
mean(jaccard.dist.among.rare$Wo_rare) #0.6605121

jaccard.dist.bac.rare$Group <- factor(jaccard.dist.bac.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(jaccard.dist.bac.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

###Fungi
#### Distance within each branch
Br2<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_2")]
Br3<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_3")]
Br4<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_4")]
Br5<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_5")]
Br6<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_6")]
Br7<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_7")]
Br8<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_8")]

## Within branch
jaccard.dist.fun.rare2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(jaccard.dist.fun.rare2$Wo_rare)
jaccard.dist.fun.rare2$Comparison <- "Branch_2"
mean(jaccard.dist.fun.rare2$Wo_rare) #0.574067

jaccard.dist.fun.rare3 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(jaccard.dist.fun.rare3$Wo_rare)
jaccard.dist.fun.rare3$Comparison <- "Branch_3"
mean(jaccard.dist.fun.rare3$Wo_rare) #0.4370992

jaccard.dist.fun.rare4 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(jaccard.dist.fun.rare4$Wo_rare)
jaccard.dist.fun.rare4$Comparison <- "Branch_4"
mean(jaccard.dist.fun.rare4$Wo_rare) #0.432106

jaccard.dist.fun.rare5 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(jaccard.dist.fun.rare5$Wo_rare)
jaccard.dist.fun.rare5$Comparison <- "Branch_5"
mean(jaccard.dist.fun.rare5$Wo_rare) #0.4373999

jaccard.dist.fun.rare6 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(jaccard.dist.fun.rare6$Wo_rare)
jaccard.dist.fun.rare6$Comparison <- "Branch_6"
mean(jaccard.dist.fun.rare6$Wo_rare) #0.5349464

jaccard.dist.fun.rare7 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(jaccard.dist.fun.rare7$Wo_rare)
jaccard.dist.fun.rare7$Comparison <- "Branch_7"
mean(jaccard.dist.fun.rare7$Wo_rare) # 0.5174208

jaccard.dist.fun.rare8 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(jaccard.dist.fun.rare8$Wo_rare)
jaccard.dist.fun.rare8$Comparison <- "Branch_8"
mean(jaccard.dist.fun.rare8$Wo_rare) #0.4159909

jaccard.dist.within.rare <- rbind(jaccard.dist.fun.rare2, jaccard.dist.fun.rare3,jaccard.dist.fun.rare4,
                                  jaccard.dist.fun.rare5,jaccard.dist.fun.rare6,jaccard.dist.fun.rare7,
                                  jaccard.dist.fun.rare8)
jaccard.dist.within.rare$Group <- "Within_branch"

##Among branch
jaccard.dist.fun.rare2.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
jaccard.dist.fun.rare2.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
jaccard.dist.fun.rare2.other <- rbind(jaccard.dist.fun.rare2.other.1, jaccard.dist.fun.rare2.other.2)
length(jaccard.dist.fun.rare2.other$Wo_rare)
jaccard.dist.fun.rare2.other$Comparison <- "Br2_Others"
mean(jaccard.dist.fun.rare2.other$Wo_rare) #0.5606695

jaccard.dist.fun.rare3.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
jaccard.dist.fun.rare3.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
jaccard.dist.fun.rare3.other <- rbind(jaccard.dist.fun.rare3.other.1, jaccard.dist.fun.rare3.other.2)
length(jaccard.dist.fun.rare3.other$Wo_rare)
jaccard.dist.fun.rare3.other$Comparison <- "Br3_Others"
mean(jaccard.dist.fun.rare3.other$Wo_rare) #0.7950433

jaccard.dist.fun.rare4.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
jaccard.dist.fun.rare4.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
jaccard.dist.fun.rare4.other <- rbind(jaccard.dist.fun.rare4.other.1, jaccard.dist.fun.rare4.other.2)
length(jaccard.dist.fun.rare4.other$Wo_rare)
jaccard.dist.fun.rare4.other$Comparison <- "Br4_Others"
mean(jaccard.dist.fun.rare4.other$Wo_rare) #0.795574

jaccard.dist.fun.rare5.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
jaccard.dist.fun.rare5.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
jaccard.dist.fun.rare5.other <- rbind(jaccard.dist.fun.rare5.other.1, jaccard.dist.fun.rare5.other.2)
length(jaccard.dist.fun.rare5.other$Wo_rare)
jaccard.dist.fun.rare5.other$Comparison <- "Br5_Others"
mean(jaccard.dist.fun.rare5.other$Wo_rare) #0.8299794

jaccard.dist.fun.rare6.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
jaccard.dist.fun.rare6.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
jaccard.dist.fun.rare6.other <- rbind(jaccard.dist.fun.rare6.other.1, jaccard.dist.fun.rare6.other.2)
length(jaccard.dist.fun.rare6.other$Wo_rare)
jaccard.dist.fun.rare6.other$Comparison <- "Br6_Others"
mean(jaccard.dist.fun.rare6.other$Wo_rare) #0.8045498

jaccard.dist.fun.rare7.other.1 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
jaccard.dist.fun.rare7.other.2 <- subset(jaccard.dist.fun.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
jaccard.dist.fun.rare7.other <- rbind(jaccard.dist.fun.rare7.other.1, jaccard.dist.fun.rare7.other.2)
length(jaccard.dist.fun.rare7.other$Wo_rare)
jaccard.dist.fun.rare7.other$Comparison <- "Br7_Others"
mean(jaccard.dist.fun.rare7.other$Wo_rare) #0.7981205

jaccard.dist.among.rare <- rbind(jaccard.dist.fun.rare2.other, jaccard.dist.fun.rare3.other,jaccard.dist.fun.rare4.other,
                                 jaccard.dist.fun.rare5.other,jaccard.dist.fun.rare6.other,jaccard.dist.fun.rare7.other)
jaccard.dist.among.rare$Group <- "Among_branch"


jaccard.dist.fun.rare <- rbind(jaccard.dist.within.rare,jaccard.dist.among.rare)

### Statistical analysis
shapiro.test(jaccard.dist.fun.rare$Wo_rare) #p-value =  3.822e-07
var.test(jaccard.dist.fun.rare$Wo_rare~ jaccard.dist.fun.rare$Group)

#wilcox.test(jaccard.dist.within.rare$Wo_rare, jaccard.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(jaccard.dist.within.rare$Wo_rare) #0.4784329
mean(jaccard.dist.among.rare$Wo_rare) #0.5150946

jaccard.dist.fun.rare$Group <- factor(jaccard.dist.fun.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(jaccard.dist.fun.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


### Unweighted unifrac
## Within branch
unifrac.dist.bac.rare2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(unifrac.dist.bac.rare2$Wo_rare)
unifrac.dist.bac.rare2$Comparison <- "Branch_2"
mean(unifrac.dist.bac.rare2$Wo_rare) #0.5264953

unifrac.dist.bac.rare3 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(unifrac.dist.bac.rare3$Wo_rare)
unifrac.dist.bac.rare3$Comparison <- "Branch_3"
mean(unifrac.dist.bac.rare3$Wo_rare) #0.6174074

unifrac.dist.bac.rare4 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(unifrac.dist.bac.rare4$Wo_rare)
unifrac.dist.bac.rare4$Comparison <- "Branch_4"
mean(unifrac.dist.bac.rare4$Wo_rare) #0.6371094

unifrac.dist.bac.rare5 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(unifrac.dist.bac.rare5$Wo_rare)
unifrac.dist.bac.rare5$Comparison <- "Branch_5"
mean(unifrac.dist.bac.rare5$Wo_rare) #0.6169486

unifrac.dist.bac.rare6 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(unifrac.dist.bac.rare6$Wo_rare)
unifrac.dist.bac.rare6$Comparison <- "Branch_6"
mean(unifrac.dist.bac.rare6$Wo_rare) #0.59317

unifrac.dist.bac.rare7 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(unifrac.dist.bac.rare7$Wo_rare)
unifrac.dist.bac.rare7$Comparison <- "Branch_7"
mean(unifrac.dist.bac.rare7$Wo_rare) # 0.5444338

unifrac.dist.bac.rare8 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(unifrac.dist.bac.rare8$Wo_rare)
unifrac.dist.bac.rare8$Comparison <- "Branch_8"
mean(unifrac.dist.bac.rare8$Wo_rare) #0.6330365

unifrac.dist.within.rare <- rbind(unifrac.dist.bac.rare2, unifrac.dist.bac.rare3,unifrac.dist.bac.rare4,
                                  unifrac.dist.bac.rare5,unifrac.dist.bac.rare6,unifrac.dist.bac.rare7,
                                  unifrac.dist.bac.rare8)
unifrac.dist.within.rare$Group <- "Within_branch"

##Among branch
unifrac.dist.bac.rare2.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
unifrac.dist.bac.rare2.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
unifrac.dist.bac.rare2.other <- rbind(unifrac.dist.bac.rare2.other.1, unifrac.dist.bac.rare2.other.2)
length(unifrac.dist.bac.rare2.other$Wo_rare)
unifrac.dist.bac.rare2.other$Comparison <- "Br2_Others"
mean(unifrac.dist.bac.rare2.other$Wo_rare) #0.7944298

unifrac.dist.bac.rare3.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
unifrac.dist.bac.rare3.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
unifrac.dist.bac.rare3.other <- rbind(unifrac.dist.bac.rare3.other.1, unifrac.dist.bac.rare3.other.2)
length(unifrac.dist.bac.rare3.other$Wo_rare)
unifrac.dist.bac.rare3.other$Comparison <- "Br3_Others"
mean(unifrac.dist.bac.rare3.other$Wo_rare) #0.7950433

unifrac.dist.bac.rare4.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
unifrac.dist.bac.rare4.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
unifrac.dist.bac.rare4.other <- rbind(unifrac.dist.bac.rare4.other.1, unifrac.dist.bac.rare4.other.2)
length(unifrac.dist.bac.rare4.other$Wo_rare)
unifrac.dist.bac.rare4.other$Comparison <- "Br4_Others"
mean(unifrac.dist.bac.rare4.other$Wo_rare) #0.795574

unifrac.dist.bac.rare5.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
unifrac.dist.bac.rare5.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
unifrac.dist.bac.rare5.other <- rbind(unifrac.dist.bac.rare5.other.1, unifrac.dist.bac.rare5.other.2)
length(unifrac.dist.bac.rare5.other$Wo_rare)
unifrac.dist.bac.rare5.other$Comparison <- "Br5_Others"
mean(unifrac.dist.bac.rare5.other$Wo_rare) #0.8299794

unifrac.dist.bac.rare6.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
unifrac.dist.bac.rare6.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
unifrac.dist.bac.rare6.other <- rbind(unifrac.dist.bac.rare6.other.1, unifrac.dist.bac.rare6.other.2)
length(unifrac.dist.bac.rare6.other$Wo_rare)
unifrac.dist.bac.rare6.other$Comparison <- "Br6_Others"
mean(unifrac.dist.bac.rare6.other$Wo_rare) #0.8045498

unifrac.dist.bac.rare7.other.1 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
unifrac.dist.bac.rare7.other.2 <- subset(unifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
unifrac.dist.bac.rare7.other <- rbind(unifrac.dist.bac.rare7.other.1, unifrac.dist.bac.rare7.other.2)
length(unifrac.dist.bac.rare7.other$Wo_rare)
unifrac.dist.bac.rare7.other$Comparison <- "Br7_Others"
mean(unifrac.dist.bac.rare7.other$Wo_rare) #0.7981205

unifrac.dist.among.rare <- rbind(unifrac.dist.bac.rare2.other, unifrac.dist.bac.rare3.other,unifrac.dist.bac.rare4.other,
                                 unifrac.dist.bac.rare5.other,unifrac.dist.bac.rare6.other,unifrac.dist.bac.rare7.other)
unifrac.dist.among.rare$Group <- "Among_branch"


unifrac.dist.bac.rare <- rbind(unifrac.dist.within.rare,unifrac.dist.among.rare)

### Statistical analysis
shapiro.test(unifrac.dist.bac.rare$Wo_rare) #p-value < 2.2e-16
var.test(unifrac.dist.bac.rare$Wo_rare~ unifrac.dist.bac.rare$Group)

#wilcox.test(unifrac.dist.within.rare$Wo_rare, unifrac.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(unifrac.dist.within.rare$Wo_rare) #0.4315543
mean(unifrac.dist.among.rare$Wo_rare) #0.4446555

unifrac.dist.bac.rare$Group <- factor(unifrac.dist.bac.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(unifrac.dist.bac.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

###Fungi
#### Distance within each branch
## Within branch
unifrac.dist.fun.rare2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(unifrac.dist.fun.rare2$Wo_rare)
unifrac.dist.fun.rare2$Comparison <- "Branch_2"
mean(unifrac.dist.fun.rare2$Wo_rare) #0.574067

unifrac.dist.fun.rare3 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(unifrac.dist.fun.rare3$Wo_rare)
unifrac.dist.fun.rare3$Comparison <- "Branch_3"
mean(unifrac.dist.fun.rare3$Wo_rare) #0.4370992

unifrac.dist.fun.rare4 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(unifrac.dist.fun.rare4$Wo_rare)
unifrac.dist.fun.rare4$Comparison <- "Branch_4"
mean(unifrac.dist.fun.rare4$Wo_rare) #0.432106

unifrac.dist.fun.rare5 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(unifrac.dist.fun.rare5$Wo_rare)
unifrac.dist.fun.rare5$Comparison <- "Branch_5"
mean(unifrac.dist.fun.rare5$Wo_rare) #0.4373999

unifrac.dist.fun.rare6 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(unifrac.dist.fun.rare6$Wo_rare)
unifrac.dist.fun.rare6$Comparison <- "Branch_6"
mean(unifrac.dist.fun.rare6$Wo_rare) #0.5349464

unifrac.dist.fun.rare7 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(unifrac.dist.fun.rare7$Wo_rare)
unifrac.dist.fun.rare7$Comparison <- "Branch_7"
mean(unifrac.dist.fun.rare7$Wo_rare) # 0.5174208

unifrac.dist.fun.rare8 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(unifrac.dist.fun.rare8$Wo_rare)
unifrac.dist.fun.rare8$Comparison <- "Branch_8"
mean(unifrac.dist.fun.rare8$Wo_rare) #0.4159909

unifrac.dist.within.rare <- rbind(unifrac.dist.fun.rare2, unifrac.dist.fun.rare3,unifrac.dist.fun.rare4,
                                  unifrac.dist.fun.rare5,unifrac.dist.fun.rare6,unifrac.dist.fun.rare7,
                                  unifrac.dist.fun.rare8)
unifrac.dist.within.rare$Group <- "Within_branch"

##Among branch
unifrac.dist.fun.rare2.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
unifrac.dist.fun.rare2.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
unifrac.dist.fun.rare2.other <- rbind(unifrac.dist.fun.rare2.other.1, unifrac.dist.fun.rare2.other.2)
length(unifrac.dist.fun.rare2.other$Wo_rare)
unifrac.dist.fun.rare2.other$Comparison <- "Br2_Others"
mean(unifrac.dist.fun.rare2.other$Wo_rare) #0.5606695

unifrac.dist.fun.rare3.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
unifrac.dist.fun.rare3.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
unifrac.dist.fun.rare3.other <- rbind(unifrac.dist.fun.rare3.other.1, unifrac.dist.fun.rare3.other.2)
length(unifrac.dist.fun.rare3.other$Wo_rare)
unifrac.dist.fun.rare3.other$Comparison <- "Br3_Others"
mean(unifrac.dist.fun.rare3.other$Wo_rare) #0.7950433

unifrac.dist.fun.rare4.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
unifrac.dist.fun.rare4.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
unifrac.dist.fun.rare4.other <- rbind(unifrac.dist.fun.rare4.other.1, unifrac.dist.fun.rare4.other.2)
length(unifrac.dist.fun.rare4.other$Wo_rare)
unifrac.dist.fun.rare4.other$Comparison <- "Br4_Others"
mean(unifrac.dist.fun.rare4.other$Wo_rare) #0.795574

unifrac.dist.fun.rare5.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
unifrac.dist.fun.rare5.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
unifrac.dist.fun.rare5.other <- rbind(unifrac.dist.fun.rare5.other.1, unifrac.dist.fun.rare5.other.2)
length(unifrac.dist.fun.rare5.other$Wo_rare)
unifrac.dist.fun.rare5.other$Comparison <- "Br5_Others"
mean(unifrac.dist.fun.rare5.other$Wo_rare) #0.8299794

unifrac.dist.fun.rare6.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
unifrac.dist.fun.rare6.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
unifrac.dist.fun.rare6.other <- rbind(unifrac.dist.fun.rare6.other.1, unifrac.dist.fun.rare6.other.2)
length(unifrac.dist.fun.rare6.other$Wo_rare)
unifrac.dist.fun.rare6.other$Comparison <- "Br6_Others"
mean(unifrac.dist.fun.rare6.other$Wo_rare) #0.8045498

unifrac.dist.fun.rare7.other.1 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
unifrac.dist.fun.rare7.other.2 <- subset(unifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
unifrac.dist.fun.rare7.other <- rbind(unifrac.dist.fun.rare7.other.1, unifrac.dist.fun.rare7.other.2)
length(unifrac.dist.fun.rare7.other$Wo_rare)
unifrac.dist.fun.rare7.other$Comparison <- "Br7_Others"
mean(unifrac.dist.fun.rare7.other$Wo_rare) #0.7981205

unifrac.dist.among.rare <- rbind(unifrac.dist.fun.rare2.other, unifrac.dist.fun.rare3.other,unifrac.dist.fun.rare4.other,
                                 unifrac.dist.fun.rare5.other,unifrac.dist.fun.rare6.other,unifrac.dist.fun.rare7.other)
unifrac.dist.among.rare$Group <- "Among_branch"


unifrac.dist.fun.rare <- rbind(unifrac.dist.within.rare,unifrac.dist.among.rare)

### Statistical analysis
shapiro.test(unifrac.dist.fun.rare$Wo_rare) #p-value =  3.822e-07
var.test(unifrac.dist.fun.rare$Wo_rare~ unifrac.dist.fun.rare$Group)

#wilcox.test(unifrac.dist.within.rare$Wo_rare, unifrac.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(unifrac.dist.within.rare$Wo_rare) #0.4784329
mean(unifrac.dist.among.rare$Wo_rare) #0.5150946

unifrac.dist.fun.rare$Group <- factor(unifrac.dist.fun.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(unifrac.dist.fun.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


#### Dissimilarity index considering ASV abundances (Bray-Curtis/ Weighted UniFrac)
### Weighted Unifrac distance (all ASVs)
## Bacteria
## All ASVs
wunifrac.dist.bac<-UniFrac(bac.clean.log, weighted = T, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.bac <-as.matrix(wunifrac.dist.bac)

wunifrac.dist.bac.melt <- melt(as.matrix(wunifrac.dist.bac), na.rm = T)
names(wunifrac.dist.bac.melt)[3] <- "All_ASVs"
wunifrac.dist.bac.melt <- subset(wunifrac.dist.bac.melt, All_ASVs != 0)
head(wunifrac.dist.bac.melt)
mean(wunifrac.dist.bac.melt$All_ASVs) #0.1850065

## without rare ASVs
bac.log.wo.rare <- pop_taxa(bac.clean.log, rare.bac)

(filt.sample <- sample_sums(bac.log.wo.rare) > 0)
sum(sample_sums(bac.log.wo.rare) <= 0)  ## 2 sample discarded
bac.log.wo.rare <- prune_samples(filt.sample, bac.log.wo.rare)
bac.log.wo.rare

wunifrac.dist.bac<-UniFrac(bac.log.wo.rare, weighted = T, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.bac <-as.matrix(wunifrac.dist.bac)

wunifrac.dist.bac.melt.wo.rare <- melt(as.matrix(wunifrac.dist.bac), na.rm = T)
names(wunifrac.dist.bac.melt.wo.rare)[3] <- "Wo_rare"
wunifrac.dist.bac.melt.wo.rare <- subset(wunifrac.dist.bac.melt.wo.rare, Wo_rare != 0)
head(wunifrac.dist.bac.melt.wo.rare)
mean(wunifrac.dist.bac.melt.wo.rare$Wo_rare) #0.3336056

## Combine two data frames
bac.dist.uw.unifrac <- merge(wunifrac.dist.bac.melt, wunifrac.dist.bac.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
bac.dist.uw.unifrac <- melt(bac.dist.uw.unifrac)
head(bac.dist.uw.unifrac)

## Normality test (Anderson-Darling normality test)
ad.test(bac.dist.uw.unifrac$value)

ggplot(bac.dist.uw.unifrac, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## Fungi
## All ASVs
wunifrac.dist.fun<-UniFrac(fun.clean.log, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.fun <-as.matrix(wunifrac.dist.fun)

wunifrac.dist.fun.melt <- melt(as.matrix(wunifrac.dist.fun), na.rm = T)
names(wunifrac.dist.fun.melt)[3] <- "All_ASVs"
wunifrac.dist.fun.melt <- subset(wunifrac.dist.fun.melt, All_ASVs != 0)
head(wunifrac.dist.fun.melt)
mean(wunifrac.dist.fun.melt$All_ASVs) #0.6903096

## without rare ASVs
fun.log.wo.rare <- pop_taxa(fun.clean.log, rare.fun)

(filt.sample <- sample_sums(fun.log.wo.rare) > 0)
sum(sample_sums(fun.log.wo.rare) <= 0)  ## 0 sample discarded
fun.log.wo.rare <- prune_samples(filt.sample, fun.log.wo.rare)
fun.log.wo.rare 

wunifrac.dist.fun<-UniFrac(fun.log.wo.rare, weighted = T, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.fun <-as.matrix(wunifrac.dist.fun)

wunifrac.dist.fun.melt.wo.rare <- melt(as.matrix(wunifrac.dist.fun), na.rm = T)
names(wunifrac.dist.fun.melt.wo.rare)[3] <- "Wo_rare"
wunifrac.dist.fun.melt.wo.rare <- subset(wunifrac.dist.fun.melt.wo.rare, Wo_rare != 0)
head(wunifrac.dist.fun.melt.wo.rare)
mean(wunifrac.dist.fun.melt.wo.rare$Wo_rare) #0.3142531

## Combine two data frames
fun.dist.uw.unifrac <- merge(wunifrac.dist.fun.melt, wunifrac.dist.fun.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
fun.dist.uw.unifrac <- melt(fun.dist.uw.unifrac)
head(fun.dist.uw.unifrac)

## Normality test (Anderson-Darling normality test)
ad.test(fun.dist.uw.unifrac$value)

ggplot(fun.dist.uw.unifrac, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



## Bray
## Bacteria
## All ASVs
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
names(bray.dist.bac.melt)[3] <- "All_ASVs"
bray.dist.bac.melt <- subset(bray.dist.bac.melt, All_ASVs != 0)


## wo Rare
b.otu.lognorm.wo.rare <- otu_table(bac.log.wo.rare)
b.otu.lognorm.wo.rare <- data.frame(b.otu.lognorm.wo.rare)
colSums(b.otu.lognorm.wo.rare)
bray.dist.bac<-vegdist(t(b.otu.lognorm.wo.rare), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
bray.dist.bac <-as.matrix(bray.dist.bac)

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt.wo.rare <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt.wo.rare)
names(bray.dist.bac.melt.wo.rare)[3] <- "Wo_rare"
bray.dist.bac.melt.wo.rare <- subset(bray.dist.bac.melt.wo.rare, Wo_rare != 0)

## Combine two data frames
bac.dist.bray <- merge(bray.dist.bac.melt, bray.dist.bac.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
bac.dist.bray <- melt(bac.dist.bray)
head(bac.dist.bray)

## Normality test (Anderson-Darling normality test)
ad.test(bac.dist.bray$value)

ggplot(bac.dist.bray, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## Fungi
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
names(bray.dist.fun.melt)[3] <- "All_ASVs"
bray.dist.fun.melt <- subset(bray.dist.fun.melt, All_ASVs != 0)


## wo Rare
f.otu.lognorm.wo.rare <- otu_table(fun.log.wo.rare)
f.otu.lognorm.wo.rare <- data.frame(f.otu.lognorm.wo.rare)
colSums(f.otu.lognorm.wo.rare)
bray.dist.fun<-vegdist(t(f.otu.lognorm.wo.rare), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
bray.dist.fun <-as.matrix(bray.dist.fun)

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt.wo.rare <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt.wo.rare)
names(bray.dist.fun.melt.wo.rare)[3] <- "Wo_rare"
bray.dist.fun.melt.wo.rare <- subset(bray.dist.fun.melt.wo.rare, Wo_rare != 0)

## Combine two data frames
fun.dist.bray <- merge(bray.dist.fun.melt, bray.dist.fun.melt.wo.rare, by = c('Var1'="Var1", "Var2"="Var2"))
fun.dist.bray <- melt(fun.dist.bray)
head(fun.dist.bray)

## Normality test (Anderson-Darling normality test)
ad.test(fun.dist.bray$value)

ggplot(fun.dist.bray, aes(x=variable, y= value, fill = variable))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("All_ASVs","Wo_rare")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Among and within branch
### weighted unifrac
## Within branch
wunifrac.dist.bac.rare2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(wunifrac.dist.bac.rare2$Wo_rare)
wunifrac.dist.bac.rare2$Comparison <- "Branch_2"
mean(wunifrac.dist.bac.rare2$Wo_rare) #0.5264953

wunifrac.dist.bac.rare3 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(wunifrac.dist.bac.rare3$Wo_rare)
wunifrac.dist.bac.rare3$Comparison <- "Branch_3"
mean(wunifrac.dist.bac.rare3$Wo_rare) #0.6174074

wunifrac.dist.bac.rare4 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(wunifrac.dist.bac.rare4$Wo_rare)
wunifrac.dist.bac.rare4$Comparison <- "Branch_4"
mean(wunifrac.dist.bac.rare4$Wo_rare) #0.6371094

wunifrac.dist.bac.rare5 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(wunifrac.dist.bac.rare5$Wo_rare)
wunifrac.dist.bac.rare5$Comparison <- "Branch_5"
mean(wunifrac.dist.bac.rare5$Wo_rare) #0.6169486

wunifrac.dist.bac.rare6 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(wunifrac.dist.bac.rare6$Wo_rare)
wunifrac.dist.bac.rare6$Comparison <- "Branch_6"
mean(wunifrac.dist.bac.rare6$Wo_rare) #0.59317

wunifrac.dist.bac.rare7 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(wunifrac.dist.bac.rare7$Wo_rare)
wunifrac.dist.bac.rare7$Comparison <- "Branch_7"
mean(wunifrac.dist.bac.rare7$Wo_rare) # 0.5444338

wunifrac.dist.bac.rare8 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(wunifrac.dist.bac.rare8$Wo_rare)
wunifrac.dist.bac.rare8$Comparison <- "Branch_8"
mean(wunifrac.dist.bac.rare8$Wo_rare) #0.6330365

wunifrac.dist.within.rare <- rbind(wunifrac.dist.bac.rare2, wunifrac.dist.bac.rare3,wunifrac.dist.bac.rare4,
                                  wunifrac.dist.bac.rare5,wunifrac.dist.bac.rare6,wunifrac.dist.bac.rare7,
                                  wunifrac.dist.bac.rare8)
wunifrac.dist.within.rare$Group <- "Within_branch"

##Among branch
wunifrac.dist.bac.rare2.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.bac.rare2.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
wunifrac.dist.bac.rare2.other <- rbind(wunifrac.dist.bac.rare2.other.1, wunifrac.dist.bac.rare2.other.2)
length(wunifrac.dist.bac.rare2.other$Wo_rare)
wunifrac.dist.bac.rare2.other$Comparison <- "Br2_Others"
mean(wunifrac.dist.bac.rare2.other$Wo_rare) #0.7944298

wunifrac.dist.bac.rare3.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.bac.rare3.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
wunifrac.dist.bac.rare3.other <- rbind(wunifrac.dist.bac.rare3.other.1, wunifrac.dist.bac.rare3.other.2)
length(wunifrac.dist.bac.rare3.other$Wo_rare)
wunifrac.dist.bac.rare3.other$Comparison <- "Br3_Others"
mean(wunifrac.dist.bac.rare3.other$Wo_rare) #0.7950433

wunifrac.dist.bac.rare4.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
wunifrac.dist.bac.rare4.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
wunifrac.dist.bac.rare4.other <- rbind(wunifrac.dist.bac.rare4.other.1, wunifrac.dist.bac.rare4.other.2)
length(wunifrac.dist.bac.rare4.other$Wo_rare)
wunifrac.dist.bac.rare4.other$Comparison <- "Br4_Others"
mean(wunifrac.dist.bac.rare4.other$Wo_rare) #0.795574

wunifrac.dist.bac.rare5.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
wunifrac.dist.bac.rare5.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
wunifrac.dist.bac.rare5.other <- rbind(wunifrac.dist.bac.rare5.other.1, wunifrac.dist.bac.rare5.other.2)
length(wunifrac.dist.bac.rare5.other$Wo_rare)
wunifrac.dist.bac.rare5.other$Comparison <- "Br5_Others"
mean(wunifrac.dist.bac.rare5.other$Wo_rare) #0.8299794

wunifrac.dist.bac.rare6.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
wunifrac.dist.bac.rare6.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
wunifrac.dist.bac.rare6.other <- rbind(wunifrac.dist.bac.rare6.other.1, wunifrac.dist.bac.rare6.other.2)
length(wunifrac.dist.bac.rare6.other$Wo_rare)
wunifrac.dist.bac.rare6.other$Comparison <- "Br6_Others"
mean(wunifrac.dist.bac.rare6.other$Wo_rare) #0.8045498

wunifrac.dist.bac.rare7.other.1 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
wunifrac.dist.bac.rare7.other.2 <- subset(wunifrac.dist.bac.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
wunifrac.dist.bac.rare7.other <- rbind(wunifrac.dist.bac.rare7.other.1, wunifrac.dist.bac.rare7.other.2)
length(wunifrac.dist.bac.rare7.other$Wo_rare)
wunifrac.dist.bac.rare7.other$Comparison <- "Br7_Others"
mean(wunifrac.dist.bac.rare7.other$Wo_rare) #0.7981205

wunifrac.dist.among.rare <- rbind(wunifrac.dist.bac.rare2.other, wunifrac.dist.bac.rare3.other,wunifrac.dist.bac.rare4.other,
                                 wunifrac.dist.bac.rare5.other,wunifrac.dist.bac.rare6.other,wunifrac.dist.bac.rare7.other)
wunifrac.dist.among.rare$Group <- "Among_branch"


wunifrac.dist.bac.rare <- rbind(wunifrac.dist.within.rare,wunifrac.dist.among.rare)

### Statistical analysis
shapiro.test(wunifrac.dist.bac.rare$Wo_rare) #p-value < 2.2e-16
var.test(wunifrac.dist.bac.rare$Wo_rare~ wunifrac.dist.bac.rare$Group)

#wilcox.test(wunifrac.dist.within.rare$Wo_rare, wunifrac.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(wunifrac.dist.within.rare$Wo_rare) #0.4315543
mean(wunifrac.dist.among.rare$Wo_rare) #0.4446555

wunifrac.dist.bac.rare$Group <- factor(wunifrac.dist.bac.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(wunifrac.dist.bac.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

###Fungi
#### Distance within each branch
## Within branch
wunifrac.dist.fun.rare2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(wunifrac.dist.fun.rare2$Wo_rare)
wunifrac.dist.fun.rare2$Comparison <- "Branch_2"
mean(wunifrac.dist.fun.rare2$Wo_rare) #0.574067

wunifrac.dist.fun.rare3 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(wunifrac.dist.fun.rare3$Wo_rare)
wunifrac.dist.fun.rare3$Comparison <- "Branch_3"
mean(wunifrac.dist.fun.rare3$Wo_rare) #0.4370992

wunifrac.dist.fun.rare4 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(wunifrac.dist.fun.rare4$Wo_rare)
wunifrac.dist.fun.rare4$Comparison <- "Branch_4"
mean(wunifrac.dist.fun.rare4$Wo_rare) #0.432106

wunifrac.dist.fun.rare5 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(wunifrac.dist.fun.rare5$Wo_rare)
wunifrac.dist.fun.rare5$Comparison <- "Branch_5"
mean(wunifrac.dist.fun.rare5$Wo_rare) #0.4373999

wunifrac.dist.fun.rare6 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(wunifrac.dist.fun.rare6$Wo_rare)
wunifrac.dist.fun.rare6$Comparison <- "Branch_6"
mean(wunifrac.dist.fun.rare6$Wo_rare) #0.5349464

wunifrac.dist.fun.rare7 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(wunifrac.dist.fun.rare7$Wo_rare)
wunifrac.dist.fun.rare7$Comparison <- "Branch_7"
mean(wunifrac.dist.fun.rare7$Wo_rare) # 0.5174208

wunifrac.dist.fun.rare8 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(wunifrac.dist.fun.rare8$Wo_rare)
wunifrac.dist.fun.rare8$Comparison <- "Branch_8"
mean(wunifrac.dist.fun.rare8$Wo_rare) #0.4159909

wunifrac.dist.within.rare <- rbind(wunifrac.dist.fun.rare2, wunifrac.dist.fun.rare3,wunifrac.dist.fun.rare4,
                                  wunifrac.dist.fun.rare5,wunifrac.dist.fun.rare6,wunifrac.dist.fun.rare7,
                                  wunifrac.dist.fun.rare8)
wunifrac.dist.within.rare$Group <- "Within_branch"

##Among branch
wunifrac.dist.fun.rare2.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.fun.rare2.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
wunifrac.dist.fun.rare2.other <- rbind(wunifrac.dist.fun.rare2.other.1, wunifrac.dist.fun.rare2.other.2)
length(wunifrac.dist.fun.rare2.other$Wo_rare)
wunifrac.dist.fun.rare2.other$Comparison <- "Br2_Others"
mean(wunifrac.dist.fun.rare2.other$Wo_rare) #0.5606695

wunifrac.dist.fun.rare3.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.fun.rare3.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
wunifrac.dist.fun.rare3.other <- rbind(wunifrac.dist.fun.rare3.other.1, wunifrac.dist.fun.rare3.other.2)
length(wunifrac.dist.fun.rare3.other$Wo_rare)
wunifrac.dist.fun.rare3.other$Comparison <- "Br3_Others"
mean(wunifrac.dist.fun.rare3.other$Wo_rare) #0.7950433

wunifrac.dist.fun.rare4.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
wunifrac.dist.fun.rare4.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
wunifrac.dist.fun.rare4.other <- rbind(wunifrac.dist.fun.rare4.other.1, wunifrac.dist.fun.rare4.other.2)
length(wunifrac.dist.fun.rare4.other$Wo_rare)
wunifrac.dist.fun.rare4.other$Comparison <- "Br4_Others"
mean(wunifrac.dist.fun.rare4.other$Wo_rare) #0.795574

wunifrac.dist.fun.rare5.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
wunifrac.dist.fun.rare5.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
wunifrac.dist.fun.rare5.other <- rbind(wunifrac.dist.fun.rare5.other.1, wunifrac.dist.fun.rare5.other.2)
length(wunifrac.dist.fun.rare5.other$Wo_rare)
wunifrac.dist.fun.rare5.other$Comparison <- "Br5_Others"
mean(wunifrac.dist.fun.rare5.other$Wo_rare) #0.8299794

wunifrac.dist.fun.rare6.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
wunifrac.dist.fun.rare6.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
wunifrac.dist.fun.rare6.other <- rbind(wunifrac.dist.fun.rare6.other.1, wunifrac.dist.fun.rare6.other.2)
length(wunifrac.dist.fun.rare6.other$Wo_rare)
wunifrac.dist.fun.rare6.other$Comparison <- "Br6_Others"
mean(wunifrac.dist.fun.rare6.other$Wo_rare) #0.8045498

wunifrac.dist.fun.rare7.other.1 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
wunifrac.dist.fun.rare7.other.2 <- subset(wunifrac.dist.fun.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
wunifrac.dist.fun.rare7.other <- rbind(wunifrac.dist.fun.rare7.other.1, wunifrac.dist.fun.rare7.other.2)
length(wunifrac.dist.fun.rare7.other$Wo_rare)
wunifrac.dist.fun.rare7.other$Comparison <- "Br7_Others"
mean(wunifrac.dist.fun.rare7.other$Wo_rare) #0.7981205

wunifrac.dist.among.rare <- rbind(wunifrac.dist.fun.rare2.other, wunifrac.dist.fun.rare3.other,wunifrac.dist.fun.rare4.other,
                                 wunifrac.dist.fun.rare5.other,wunifrac.dist.fun.rare6.other,wunifrac.dist.fun.rare7.other)
wunifrac.dist.among.rare$Group <- "Among_branch"


wunifrac.dist.fun.rare <- rbind(wunifrac.dist.within.rare,wunifrac.dist.among.rare)

### Statistical analysis
shapiro.test(wunifrac.dist.fun.rare$Wo_rare) #p-value =  3.822e-07
var.test(wunifrac.dist.fun.rare$Wo_rare~ wunifrac.dist.fun.rare$Group)

#wilcox.test(wunifrac.dist.within.rare$Wo_rare, wunifrac.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(wunifrac.dist.within.rare$Wo_rare) #0.4784329
mean(wunifrac.dist.among.rare$Wo_rare) #0.5150946

wunifrac.dist.fun.rare$Group <- factor(wunifrac.dist.fun.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(wunifrac.dist.fun.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##Bray
## Within branch
bray.dist.bac.rare2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(bray.dist.bac.rare2$Wo_rare)
bray.dist.bac.rare2$Comparison <- "Branch_2"
mean(bray.dist.bac.rare2$Wo_rare) #0.5264953

bray.dist.bac.rare3 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(bray.dist.bac.rare3$Wo_rare)
bray.dist.bac.rare3$Comparison <- "Branch_3"
mean(bray.dist.bac.rare3$Wo_rare) #0.6174074

bray.dist.bac.rare4 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(bray.dist.bac.rare4$Wo_rare)
bray.dist.bac.rare4$Comparison <- "Branch_4"
mean(bray.dist.bac.rare4$Wo_rare) #0.6371094

bray.dist.bac.rare5 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(bray.dist.bac.rare5$Wo_rare)
bray.dist.bac.rare5$Comparison <- "Branch_5"
mean(bray.dist.bac.rare5$Wo_rare) #0.6169486

bray.dist.bac.rare6 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(bray.dist.bac.rare6$Wo_rare)
bray.dist.bac.rare6$Comparison <- "Branch_6"
mean(bray.dist.bac.rare6$Wo_rare) #0.59317

bray.dist.bac.rare7 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(bray.dist.bac.rare7$Wo_rare)
bray.dist.bac.rare7$Comparison <- "Branch_7"
mean(bray.dist.bac.rare7$Wo_rare) # 0.5444338

bray.dist.bac.rare8 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(bray.dist.bac.rare8$Wo_rare)
bray.dist.bac.rare8$Comparison <- "Branch_8"
mean(bray.dist.bac.rare8$Wo_rare) #0.6330365

bray.dist.within.rare <- rbind(bray.dist.bac.rare2, bray.dist.bac.rare3,bray.dist.bac.rare4,
                                   bray.dist.bac.rare5,bray.dist.bac.rare6,bray.dist.bac.rare7,
                                   bray.dist.bac.rare8)
bray.dist.within.rare$Group <- "Within_branch"

##Among branch
bray.dist.bac.rare2.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
bray.dist.bac.rare2.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
bray.dist.bac.rare2.other <- rbind(bray.dist.bac.rare2.other.1, bray.dist.bac.rare2.other.2)
length(bray.dist.bac.rare2.other$Wo_rare)
bray.dist.bac.rare2.other$Comparison <- "Br2_Others"
mean(bray.dist.bac.rare2.other$Wo_rare) #0.7944298

bray.dist.bac.rare3.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
bray.dist.bac.rare3.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
bray.dist.bac.rare3.other <- rbind(bray.dist.bac.rare3.other.1, bray.dist.bac.rare3.other.2)
length(bray.dist.bac.rare3.other$Wo_rare)
bray.dist.bac.rare3.other$Comparison <- "Br3_Others"
mean(bray.dist.bac.rare3.other$Wo_rare) #0.7950433

bray.dist.bac.rare4.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
bray.dist.bac.rare4.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
bray.dist.bac.rare4.other <- rbind(bray.dist.bac.rare4.other.1, bray.dist.bac.rare4.other.2)
length(bray.dist.bac.rare4.other$Wo_rare)
bray.dist.bac.rare4.other$Comparison <- "Br4_Others"
mean(bray.dist.bac.rare4.other$Wo_rare) #0.795574

bray.dist.bac.rare5.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
bray.dist.bac.rare5.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
bray.dist.bac.rare5.other <- rbind(bray.dist.bac.rare5.other.1, bray.dist.bac.rare5.other.2)
length(bray.dist.bac.rare5.other$Wo_rare)
bray.dist.bac.rare5.other$Comparison <- "Br5_Others"
mean(bray.dist.bac.rare5.other$Wo_rare) #0.8299794

bray.dist.bac.rare6.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
bray.dist.bac.rare6.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
bray.dist.bac.rare6.other <- rbind(bray.dist.bac.rare6.other.1, bray.dist.bac.rare6.other.2)
length(bray.dist.bac.rare6.other$Wo_rare)
bray.dist.bac.rare6.other$Comparison <- "Br6_Others"
mean(bray.dist.bac.rare6.other$Wo_rare) #0.8045498

bray.dist.bac.rare7.other.1 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
bray.dist.bac.rare7.other.2 <- subset(bray.dist.bac.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
bray.dist.bac.rare7.other <- rbind(bray.dist.bac.rare7.other.1, bray.dist.bac.rare7.other.2)
length(bray.dist.bac.rare7.other$Wo_rare)
bray.dist.bac.rare7.other$Comparison <- "Br7_Others"
mean(bray.dist.bac.rare7.other$Wo_rare) #0.7981205

bray.dist.among.rare <- rbind(bray.dist.bac.rare2.other, bray.dist.bac.rare3.other,bray.dist.bac.rare4.other,
                                  bray.dist.bac.rare5.other,bray.dist.bac.rare6.other,bray.dist.bac.rare7.other)
bray.dist.among.rare$Group <- "Among_branch"


bray.dist.bac.rare <- rbind(bray.dist.within.rare,bray.dist.among.rare)

### Statistical analysis
shapiro.test(bray.dist.bac.rare$Wo_rare) #p-value < 2.2e-16
var.test(bray.dist.bac.rare$Wo_rare~ bray.dist.bac.rare$Group)

#wilcox.test(bray.dist.within.rare$Wo_rare, bray.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(bray.dist.within.rare$Wo_rare) #0.4315543
mean(bray.dist.among.rare$Wo_rare) #0.4446555

bray.dist.bac.rare$Group <- factor(bray.dist.bac.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(bray.dist.bac.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

###Fungi
#### Distance within each branch
## Within branch
bray.dist.fun.rare2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% Br2)
length(bray.dist.fun.rare2$Wo_rare)
bray.dist.fun.rare2$Comparison <- "Branch_2"
mean(bray.dist.fun.rare2$Wo_rare) #0.574067

bray.dist.fun.rare3 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% Br3)
length(bray.dist.fun.rare3$Wo_rare)
bray.dist.fun.rare3$Comparison <- "Branch_3"
mean(bray.dist.fun.rare3$Wo_rare) #0.4370992

bray.dist.fun.rare4 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% Br4)
length(bray.dist.fun.rare4$Wo_rare)
bray.dist.fun.rare4$Comparison <- "Branch_4"
mean(bray.dist.fun.rare4$Wo_rare) #0.432106

bray.dist.fun.rare5 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% Br5)
length(bray.dist.fun.rare5$Wo_rare)
bray.dist.fun.rare5$Comparison <- "Branch_5"
mean(bray.dist.fun.rare5$Wo_rare) #0.4373999

bray.dist.fun.rare6 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% Br6)
length(bray.dist.fun.rare6$Wo_rare)
bray.dist.fun.rare6$Comparison <- "Branch_6"
mean(bray.dist.fun.rare6$Wo_rare) #0.5349464

bray.dist.fun.rare7 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% Br7)
length(bray.dist.fun.rare7$Wo_rare)
bray.dist.fun.rare7$Comparison <- "Branch_7"
mean(bray.dist.fun.rare7$Wo_rare) # 0.5174208

bray.dist.fun.rare8 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br8 & Var2 %in% Br8)
length(bray.dist.fun.rare8$Wo_rare)
bray.dist.fun.rare8$Comparison <- "Branch_8"
mean(bray.dist.fun.rare8$Wo_rare) #0.4159909

bray.dist.within.rare <- rbind(bray.dist.fun.rare2, bray.dist.fun.rare3,bray.dist.fun.rare4,
                                   bray.dist.fun.rare5,bray.dist.fun.rare6,bray.dist.fun.rare7,
                                   bray.dist.fun.rare8)
bray.dist.within.rare$Group <- "Within_branch"

##Among branch
bray.dist.fun.rare2.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
bray.dist.fun.rare2.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
bray.dist.fun.rare2.other <- rbind(bray.dist.fun.rare2.other.1, bray.dist.fun.rare2.other.2)
length(bray.dist.fun.rare2.other$Wo_rare)
bray.dist.fun.rare2.other$Comparison <- "Br2_Others"
mean(bray.dist.fun.rare2.other$Wo_rare) #0.5606695

bray.dist.fun.rare3.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
bray.dist.fun.rare3.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
bray.dist.fun.rare3.other <- rbind(bray.dist.fun.rare3.other.1, bray.dist.fun.rare3.other.2)
length(bray.dist.fun.rare3.other$Wo_rare)
bray.dist.fun.rare3.other$Comparison <- "Br3_Others"
mean(bray.dist.fun.rare3.other$Wo_rare) #0.7950433

bray.dist.fun.rare4.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
bray.dist.fun.rare4.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
bray.dist.fun.rare4.other <- rbind(bray.dist.fun.rare4.other.1, bray.dist.fun.rare4.other.2)
length(bray.dist.fun.rare4.other$Wo_rare)
bray.dist.fun.rare4.other$Comparison <- "Br4_Others"
mean(bray.dist.fun.rare4.other$Wo_rare) #0.795574

bray.dist.fun.rare5.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
bray.dist.fun.rare5.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
bray.dist.fun.rare5.other <- rbind(bray.dist.fun.rare5.other.1, bray.dist.fun.rare5.other.2)
length(bray.dist.fun.rare5.other$Wo_rare)
bray.dist.fun.rare5.other$Comparison <- "Br5_Others"
mean(bray.dist.fun.rare5.other$Wo_rare) #0.8299794

bray.dist.fun.rare6.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
bray.dist.fun.rare6.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
bray.dist.fun.rare6.other <- rbind(bray.dist.fun.rare6.other.1, bray.dist.fun.rare6.other.2)
length(bray.dist.fun.rare6.other$Wo_rare)
bray.dist.fun.rare6.other$Comparison <- "Br6_Others"
mean(bray.dist.fun.rare6.other$Wo_rare) #0.8045498

bray.dist.fun.rare7.other.1 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% Br7 & Var2 %in% c(Br8))
bray.dist.fun.rare7.other.2 <- subset(bray.dist.fun.melt.wo.rare, Var1 %in% c(Br8) & Var2 %in% Br7)
bray.dist.fun.rare7.other <- rbind(bray.dist.fun.rare7.other.1, bray.dist.fun.rare7.other.2)
length(bray.dist.fun.rare7.other$Wo_rare)
bray.dist.fun.rare7.other$Comparison <- "Br7_Others"
mean(bray.dist.fun.rare7.other$Wo_rare) #0.7981205

bray.dist.among.rare <- rbind(bray.dist.fun.rare2.other, bray.dist.fun.rare3.other,bray.dist.fun.rare4.other,
                                  bray.dist.fun.rare5.other,bray.dist.fun.rare6.other,bray.dist.fun.rare7.other)
bray.dist.among.rare$Group <- "Among_branch"


bray.dist.fun.rare <- rbind(bray.dist.within.rare,bray.dist.among.rare)

### Statistical analysis
shapiro.test(bray.dist.fun.rare$Wo_rare) #p-value =  3.822e-07
var.test(bray.dist.fun.rare$Wo_rare~ bray.dist.fun.rare$Group)

#wilcox.test(bray.dist.within.rare$Wo_rare, bray.dist.among.rare$Wo_rare, alternative = "two.sided") #p-value = 2.05e-06
mean(bray.dist.within.rare$Wo_rare) #0.4784329
mean(bray.dist.among.rare$Wo_rare) #0.5150946

bray.dist.fun.rare$Group <- factor(bray.dist.fun.rare$Group, levels = c("Within_branch","Among_branch"))

ggplot(bray.dist.fun.rare, aes(x=Group, y= Wo_rare, fill = Group))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within_branch","Among_branch")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



### Relative abundance without rare
bac.clean.rel
fun.clean.rel <- transform(fun.clean.ss.f, transform = "compositional")

tax_table(fun.clean.rel.wo.rare)

bac.clean.rel.wo.rare <- pop_taxa(bac.clean.rel, rare.bac)
fun.clean.rel.wo.rare <- pop_taxa(fun.clean.rel, rare.fun)

df.phylum <- bac.clean.rel.wo.rare %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
#df.phylum$Phylum2[which(df.phylum$Class=="Betaproteobacteria")] <- "Betaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

head(df.phylum)
df.phylum$SampleID <- factor(df.phylum$SampleID, levels = order.sample)

df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))
unique(df.phylum$Phylum2)

levels(df.phylum$Phylum2)
levels(df.phylum$Phylum2) = c(levels(df.phylum$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(SampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100)  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.5,]$Phylum2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.proteobacteria <- c("Deltaproteobacteria","Gammaproteobacteria","Alphaproteobacteria")
vec.reorder <- append(vec.uniden.Low, vec.order)
vec.reorder <- append(vec.reorder, vec.proteobacteria)

df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=SampleID, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Deltaproteobacteria" = "darkolivegreen4", "Proteobacteria" = "#003333","Actinobacteriota"="indianred2",
                               "Firmicutes" ="tan1","Bacteroidota"="steelblue1","unidentified" = "black", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = F))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

dev.off()



df.class.fun <- fun.clean.rel.wo.rare %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100)  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 0.01,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 
df.class.fun.rel$SampleID <- factor(df.class.fun.rel$SampleID, levels = order.sample) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=SampleID, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Dothideomycetes" = "#5195D1","Ustilaginomycetes"="#ffcc33","Sordariomycetes"= "#1E63AF",
                               "Tremellomycetes"="#BE4146","Cystobasidiomycetes" = "#A871AE","Microbotryomycetes" ="#DC9A9E",
                               "Agaricomycetes" = "#CC6C71","Eurotiomycetes"= "#6DA9DC","Malasseziomycetes" = "#99cc99",
                               "Leotiomycetes"= "#11335F", "Agaricostilbomycetes"="#003333",
                               "unidentified" ="#000000", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = F))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.class.fun.rel.p1

dev.off()
