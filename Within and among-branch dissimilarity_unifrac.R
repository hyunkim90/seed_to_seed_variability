### Within-branch dissimilarity and Among-branch dissimilarity
##Correct metadata
b.meta <- read.table("./Bacteria/sample_metadata.tsv", sep = '\t', header =T)
f.meta <- read.table("./Fungi/sample_metadata.tsv", sep = '\t', header =T)

rownames(b.meta) <- b.meta$SampleID 
rownames(f.meta) <- f.meta$SampleID 

sample_data(bac.clean.ss) <- sample_data(b.meta)
sample_data(bac.clean.nolog) <- sample_data(b.meta)
sample_data(bac.clean.log) <- sample_data(b.meta)
sample_data(bac.rarefied) <- sample_data(b.meta)

sample_data(fun.clean.ss) <- sample_data(f.meta)
sample_data(fun.clean.nolog) <- sample_data(f.meta)
sample_data(fun.clean.log) <- sample_data(f.meta)
sample_data(fun.rarefied) <- sample_data(f.meta)

#### unweighted unifrac dissimilarity - Bacteria
unifrac.dist.bac<-UniFrac(bac.clean.log, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.bac <-as.matrix(unifrac.dist.bac)

unifrac.dist.bac.melt <- melt(as.matrix(unifrac.dist.bac), na.rm = T)
head(unifrac.dist.bac.melt)
names(unifrac.dist.bac.melt)[3] <- "UniFrac_distance"
unifrac.dist.bac.melt <- subset(unifrac.dist.bac.melt, UniFrac_distance != 0)

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
unifrac.dist.bac.br2 <- subset(unifrac.dist.bac.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(unifrac.dist.bac.br2$UniFrac_distance)
unifrac.dist.bac.br2$Comparison <- "Branch_2"
mean(unifrac.dist.bac.br2$UniFrac_distance) #0.641191

unifrac.dist.bac.br3 <- subset(unifrac.dist.bac.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(unifrac.dist.bac.br3$UniFrac_distance)
unifrac.dist.bac.br3$Comparison <- "Branch_3"
mean(unifrac.dist.bac.br3$UniFrac_distance) #0.5424351

unifrac.dist.bac.br4 <- subset(unifrac.dist.bac.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(unifrac.dist.bac.br4$UniFrac_distance)
unifrac.dist.bac.br4$Comparison <- "Branch_4"
mean(unifrac.dist.bac.br4$UniFrac_distance) #0.5685465

unifrac.dist.bac.br5 <- subset(unifrac.dist.bac.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(unifrac.dist.bac.br5$UniFrac_distance)
unifrac.dist.bac.br5$Comparison <- "Branch_5"
mean(unifrac.dist.bac.br5$UniFrac_distance) #0.5556779

unifrac.dist.bac.br6 <- subset(unifrac.dist.bac.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(unifrac.dist.bac.br6$UniFrac_distance)
unifrac.dist.bac.br6$Comparison <- "Branch_6"
mean(unifrac.dist.bac.br6$UniFrac_distance) #0.5653352

unifrac.dist.bac.br7 <- subset(unifrac.dist.bac.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(unifrac.dist.bac.br7$UniFrac_distance)
unifrac.dist.bac.br7$Comparison <- "Branch_7"
mean(unifrac.dist.bac.br7$UniFrac_distance) #0.5303863

unifrac.dist.bac.br8 <- subset(unifrac.dist.bac.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(unifrac.dist.bac.br8$UniFrac_distance)
unifrac.dist.bac.br8$Comparison <- "Branch_8"
mean(unifrac.dist.bac.br8$UniFrac_distance) #0.6769637

unifrac.dist.within.br <- rbind(unifrac.dist.bac.br2, unifrac.dist.bac.br3,unifrac.dist.bac.br4,
                                unifrac.dist.bac.br5,unifrac.dist.bac.br6,unifrac.dist.bac.br7,
                                unifrac.dist.bac.br8)
unifrac.dist.within.br$Group <- "Within_branch"

##Among branch
unifrac.dist.bac.br2.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
unifrac.dist.bac.br2.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
unifrac.dist.bac.br2.other <- rbind(unifrac.dist.bac.br2.other.1, unifrac.dist.bac.br2.other.2)
length(unifrac.dist.bac.br2.other$UniFrac_distance)
unifrac.dist.bac.br2.other$Comparison <- "Br2_Others"
mean(unifrac.dist.bac.br2.other$UniFrac_distance) #0.6351267

unifrac.dist.bac.br3.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
unifrac.dist.bac.br3.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
unifrac.dist.bac.br3.other <- rbind(unifrac.dist.bac.br3.other.1, unifrac.dist.bac.br3.other.2)
length(unifrac.dist.bac.br3.other$UniFrac_distance)
unifrac.dist.bac.br3.other$Comparison <- "Br3_Others"
mean(unifrac.dist.bac.br3.other$UniFrac_distance) #0.5796139

unifrac.dist.bac.br4.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
unifrac.dist.bac.br4.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
unifrac.dist.bac.br4.other <- rbind(unifrac.dist.bac.br4.other.1, unifrac.dist.bac.br4.other.2)
length(unifrac.dist.bac.br4.other$UniFrac_distance)
unifrac.dist.bac.br4.other$Comparison <- "Br4_Others"
mean(unifrac.dist.bac.br4.other$UniFrac_distance) #0.5828701

unifrac.dist.bac.br5.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
unifrac.dist.bac.br5.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
unifrac.dist.bac.br5.other <- rbind(unifrac.dist.bac.br5.other.1, unifrac.dist.bac.br5.other.2)
length(unifrac.dist.bac.br5.other$UniFrac_distance)
unifrac.dist.bac.br5.other$Comparison <- "Br5_Others"
mean(unifrac.dist.bac.br5.other$UniFrac_distance) #0.5938304

unifrac.dist.bac.br6.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
unifrac.dist.bac.br6.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
unifrac.dist.bac.br6.other <- rbind(unifrac.dist.bac.br6.other.1, unifrac.dist.bac.br6.other.2)
length(unifrac.dist.bac.br6.other$UniFrac_distance)
unifrac.dist.bac.br6.other$Comparison <- "Br6_Others"
mean(unifrac.dist.bac.br6.other$UniFrac_distance) #0.5991597

unifrac.dist.bac.br7.other.1 <- subset(unifrac.dist.bac.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
unifrac.dist.bac.br7.other.2 <- subset(unifrac.dist.bac.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
unifrac.dist.bac.br7.other <- rbind(unifrac.dist.bac.br7.other.1, unifrac.dist.bac.br7.other.2)
length(unifrac.dist.bac.br7.other$UniFrac_distance)
unifrac.dist.bac.br7.other$Comparison <- "Br7_Others"
mean(unifrac.dist.bac.br7.other$UniFrac_distance) #0.6220263

unifrac.dist.among.br <- rbind(unifrac.dist.bac.br2.other, unifrac.dist.bac.br3.other,unifrac.dist.bac.br4.other,
                               unifrac.dist.bac.br5.other,unifrac.dist.bac.br6.other,unifrac.dist.bac.br7.other)
unifrac.dist.among.br$Group <- "Among_branch"


unifrac.dist.bac.br <- rbind(unifrac.dist.within.br,unifrac.dist.among.br)

### Statistical analysis
shapiro.test(unifrac.dist.bac.br$UniFrac_distance) #p-value < 2.2e-16
var.test(unifrac.dist.bac.br$UniFrac_distance~ unifrac.dist.bac.br$Group)

#wilcox.test(unifrac.dist.within.br$UniFrac_distance, unifrac.dist.among.br$UniFrac_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(unifrac.dist.within.br$UniFrac_distance) #0.5829337
mean(unifrac.dist.among.br$UniFrac_distance) #0.602007

unifrac.dist.bac.br$Group <- factor(unifrac.dist.bac.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(unifrac.dist.bac.br, aes(x=Group, y= UniFrac_distance, fill = Group))+
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

#### UniFrac dissimilarity - Fungi
unifrac.dist.fun<-UniFrac(fun.clean.log, weighted = F, normalized = T, parallel = FALSE, fast = TRUE)
unifrac.dist.fun <-as.matrix(unifrac.dist.fun)

unifrac.dist.fun.melt <- melt(as.matrix(unifrac.dist.fun), na.rm = T)
head(unifrac.dist.fun.melt)
names(unifrac.dist.fun.melt)[3] <- "UniFrac_distance"
unifrac.dist.fun.melt <- subset(unifrac.dist.fun.melt, UniFrac_distance != 0)

#### Distance within each branch
Br2<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_2")]
Br3<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_3")]
Br4<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_4")]
Br5<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_5")]
Br6<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_6")]
Br7<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_7")]
Br8<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_8")]

## Within branch
unifrac.dist.fun.br2 <- subset(unifrac.dist.fun.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(unifrac.dist.fun.br2$UniFrac_distance)
unifrac.dist.fun.br2$Comparison <- "Branch_2"
mean(unifrac.dist.fun.br2$UniFrac_distance) # 0.7087973

unifrac.dist.fun.br3 <- subset(unifrac.dist.fun.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(unifrac.dist.fun.br3$UniFrac_distance)
unifrac.dist.fun.br3$Comparison <- "Branch_3"
mean(unifrac.dist.fun.br3$UniFrac_distance) #0.6108535

unifrac.dist.fun.br4 <- subset(unifrac.dist.fun.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(unifrac.dist.fun.br4$UniFrac_distance)
unifrac.dist.fun.br4$Comparison <- "Branch_4"
mean(unifrac.dist.fun.br4$UniFrac_distance) #0.666533

unifrac.dist.fun.br5 <- subset(unifrac.dist.fun.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(unifrac.dist.fun.br5$UniFrac_distance)
unifrac.dist.fun.br5$Comparison <- "Branch_5"
mean(unifrac.dist.fun.br5$UniFrac_distance) #0.6373793

unifrac.dist.fun.br6 <- subset(unifrac.dist.fun.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(unifrac.dist.fun.br6$UniFrac_distance)
unifrac.dist.fun.br6$Comparison <- "Branch_6"
mean(unifrac.dist.fun.br6$UniFrac_distance) #0.7145618

unifrac.dist.fun.br7 <- subset(unifrac.dist.fun.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(unifrac.dist.fun.br7$UniFrac_distance)
unifrac.dist.fun.br7$Comparison <- "Branch_7"
mean(unifrac.dist.fun.br7$UniFrac_distance) # 0.6561975

unifrac.dist.fun.br8 <- subset(unifrac.dist.fun.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(unifrac.dist.fun.br8$UniFrac_distance)
unifrac.dist.fun.br8$Comparison <- "Branch_8"
mean(unifrac.dist.fun.br8$UniFrac_distance) #0.6203458

unifrac.dist.within.br.f <- rbind(unifrac.dist.fun.br2, unifrac.dist.fun.br3,unifrac.dist.fun.br4,
                                  unifrac.dist.fun.br5,unifrac.dist.fun.br6,unifrac.dist.fun.br7,
                                  unifrac.dist.fun.br8)
unifrac.dist.within.br.f$Group <- "Within_branch"

##Among branch
unifrac.dist.fun.br2.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
unifrac.dist.fun.br2.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
unifrac.dist.fun.br2.other <- rbind(unifrac.dist.fun.br2.other.1, unifrac.dist.fun.br2.other.2)
length(unifrac.dist.fun.br2.other$UniFrac_distance)
unifrac.dist.fun.br2.other$Comparison <- "Br2_Others"
mean(unifrac.dist.fun.br2.other$UniFrac_distance) #0.7111592

unifrac.dist.fun.br3.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
unifrac.dist.fun.br3.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
unifrac.dist.fun.br3.other <- rbind(unifrac.dist.fun.br3.other.1, unifrac.dist.fun.br3.other.2)
length(unifrac.dist.fun.br3.other$UniFrac_distance)
unifrac.dist.fun.br3.other$Comparison <- "Br3_Others"
mean(unifrac.dist.fun.br3.other$UniFrac_distance) #0.6809353

unifrac.dist.fun.br4.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
unifrac.dist.fun.br4.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
unifrac.dist.fun.br4.other <- rbind(unifrac.dist.fun.br4.other.1, unifrac.dist.fun.br4.other.2)
length(unifrac.dist.fun.br4.other$UniFrac_distance)
unifrac.dist.fun.br4.other$Comparison <- "Br4_Others"
mean(unifrac.dist.fun.br4.other$UniFrac_distance) #0.6908744

unifrac.dist.fun.br5.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
unifrac.dist.fun.br5.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
unifrac.dist.fun.br5.other <- rbind(unifrac.dist.fun.br5.other.1, unifrac.dist.fun.br5.other.2)
length(unifrac.dist.fun.br5.other$UniFrac_distance)
unifrac.dist.fun.br5.other$Comparison <- "Br5_Others"
mean(unifrac.dist.fun.br5.other$UniFrac_distance) # 0.6832181

unifrac.dist.fun.br6.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
unifrac.dist.fun.br6.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
unifrac.dist.fun.br6.other <- rbind(unifrac.dist.fun.br6.other.1, unifrac.dist.fun.br6.other.2)
length(unifrac.dist.fun.br6.other$UniFrac_distance)
unifrac.dist.fun.br6.other$Comparison <- "Br6_Others"
mean(unifrac.dist.fun.br6.other$UniFrac_distance) #0.7075816

unifrac.dist.fun.br7.other.1 <- subset(unifrac.dist.fun.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
unifrac.dist.fun.br7.other.2 <- subset(unifrac.dist.fun.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
unifrac.dist.fun.br7.other <- rbind(unifrac.dist.fun.br7.other.1, unifrac.dist.fun.br7.other.2)
length(unifrac.dist.fun.br7.other$UniFrac_distance)
unifrac.dist.fun.br7.other$Comparison <- "Br7_Others"
mean(unifrac.dist.fun.br7.other$UniFrac_distance) #0.6944304

unifrac.dist.among.br.f <- rbind(unifrac.dist.fun.br2.other, unifrac.dist.fun.br3.other,unifrac.dist.fun.br4.other,
                                 unifrac.dist.fun.br5.other,unifrac.dist.fun.br6.other,unifrac.dist.fun.br7.other)
unifrac.dist.among.br.f$Group <- "Among_branch"


unifrac.dist.fun.br <- rbind(unifrac.dist.within.br.f,unifrac.dist.among.br.f)

### Statistical analysis
shapiro.test(unifrac.dist.fun.br$UniFrac_distance) #p-value < 2.2e-16
var.test(unifrac.dist.fun.br$UniFrac_distance~ unifrac.dist.fun.br$Group)
#wilcox.test(unifrac.dist.within.br.f$UniFrac_distance, unifrac.dist.among.br.f$UniFrac_distance, alternative = "two.sided") #p-value = 2.05e-06
mean(unifrac.dist.within.br.f$UniFrac_distance) #0.6109445
mean(unifrac.dist.among.br.f$UniFrac_distance) #0.6426395

unifrac.dist.fun.br$Group <- factor(unifrac.dist.fun.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(unifrac.dist.fun.br, aes(x=Group, y= UniFrac_distance, fill = Group))+
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

### Weighted Unifrac
wunifrac.dist.bac<-UniFrac(bac.clean.log, weighted = T, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.bac <-as.matrix(wunifrac.dist.bac)

wunifrac.dist.bac <-as.matrix(wunifrac.dist.bac)
wunifrac.dist.bac.melt <- melt(as.matrix(wunifrac.dist.bac), na.rm = T)
head(wunifrac.dist.bac.melt)
names(wunifrac.dist.bac.melt)[3] <- "WUniFracdistance"
wunifrac.dist.bac.melt <- subset(wunifrac.dist.bac.melt, WUniFracdistance != 0)

#### Distance within each branch
Br2<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_2")]
Br3<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_3")]
Br4<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_4")]
Br5<-b.meta$SampleID[which(b.meta$Branch_number  == "Branch_5")]
Br6<-b.meta$SampleID[which( b.meta$Branch_number  == "Branch_6")]
Br7<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_7")]
Br8<-b.meta$SampleID[which( b.meta$Branch_number   == "Branch_8")]

## Within branch
wunifrac.dist.bac.br2 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(wunifrac.dist.bac.br2$WUniFracdistance)
wunifrac.dist.bac.br2$Comparison <- "Branch_2"
mean(wunifrac.dist.bac.br2$WUniFracdistance) #0.1669823

wunifrac.dist.bac.br3 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(wunifrac.dist.bac.br3$WUniFracdistance)
wunifrac.dist.bac.br3$Comparison <- "Branch_3"
mean(wunifrac.dist.bac.br3$WUniFracdistance) #0.1157155

wunifrac.dist.bac.br4 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(wunifrac.dist.bac.br4$WUniFracdistance)
wunifrac.dist.bac.br4$Comparison <- "Branch_4"
mean(wunifrac.dist.bac.br4$WUniFracdistance) # 0.123696

wunifrac.dist.bac.br5 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(wunifrac.dist.bac.br5$WUniFracdistance)
wunifrac.dist.bac.br5$Comparison <- "Branch_5"
mean(wunifrac.dist.bac.br5$WUniFracdistance) #0.1422527

wunifrac.dist.bac.br6 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(wunifrac.dist.bac.br6$WUniFracdistance)
wunifrac.dist.bac.br6$Comparison <- "Branch_6"
mean(wunifrac.dist.bac.br6$WUniFracdistance) #0.1240462

wunifrac.dist.bac.br7 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(wunifrac.dist.bac.br7$WUniFracdistance)
wunifrac.dist.bac.br7$Comparison <- "Branch_7"
mean(wunifrac.dist.bac.br7$WUniFracdistance) #0.1191414

wunifrac.dist.bac.br8 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(wunifrac.dist.bac.br8$WUniFracdistance)
wunifrac.dist.bac.br8$Comparison <- "Branch_8"
mean(wunifrac.dist.bac.br8$WUniFracdistance) #0.4321489

wunifrac.dist.within.br <- rbind(wunifrac.dist.bac.br2, wunifrac.dist.bac.br3,wunifrac.dist.bac.br4,
                             wunifrac.dist.bac.br5,wunifrac.dist.bac.br6,wunifrac.dist.bac.br7,
                             wunifrac.dist.bac.br8)
wunifrac.dist.within.br$Group <- "Within_branch"

##Among branch
wunifrac.dist.bac.br2.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.bac.br2.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
wunifrac.dist.bac.br2.other <- rbind(wunifrac.dist.bac.br2.other.1, wunifrac.dist.bac.br2.other.2)
length(wunifrac.dist.bac.br2.other$WUniFracdistance)
wunifrac.dist.bac.br2.other$Comparison <- "Br2_Others"
mean(wunifrac.dist.bac.br2.other$WUniFracdistance) #0.1794707

wunifrac.dist.bac.br3.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.bac.br3.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
wunifrac.dist.bac.br3.other <- rbind(wunifrac.dist.bac.br3.other.1, wunifrac.dist.bac.br3.other.2)
length(wunifrac.dist.bac.br3.other$WUniFracdistance)
wunifrac.dist.bac.br3.other$Comparison <- "Br3_Others"
mean(wunifrac.dist.bac.br3.other$WUniFracdistance) #0.1626714

wunifrac.dist.bac.br4.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
wunifrac.dist.bac.br4.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
wunifrac.dist.bac.br4.other <- rbind(wunifrac.dist.bac.br4.other.1, wunifrac.dist.bac.br4.other.2)
length(wunifrac.dist.bac.br4.other$WUniFracdistance)
wunifrac.dist.bac.br4.other$Comparison <- "Br4_Others"
mean(wunifrac.dist.bac.br4.other$WUniFracdistance) #0.1734819

wunifrac.dist.bac.br5.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
wunifrac.dist.bac.br5.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
wunifrac.dist.bac.br5.other <- rbind(wunifrac.dist.bac.br5.other.1, wunifrac.dist.bac.br5.other.2)
length(wunifrac.dist.bac.br5.other$WUniFracdistance)
wunifrac.dist.bac.br5.other$Comparison <- "Br5_Others"
mean(wunifrac.dist.bac.br5.other$WUniFracdistance) #0.2021167

wunifrac.dist.bac.br6.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
wunifrac.dist.bac.br6.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
wunifrac.dist.bac.br6.other <- rbind(wunifrac.dist.bac.br6.other.1, wunifrac.dist.bac.br6.other.2)
length(wunifrac.dist.bac.br6.other$WUniFracdistance)
wunifrac.dist.bac.br6.other$Comparison <- "Br6_Others"
mean(wunifrac.dist.bac.br6.other$WUniFracdistance) #0.2147671

wunifrac.dist.bac.br7.other.1 <- subset(wunifrac.dist.bac.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
wunifrac.dist.bac.br7.other.2 <- subset(wunifrac.dist.bac.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
wunifrac.dist.bac.br7.other <- rbind(wunifrac.dist.bac.br7.other.1, wunifrac.dist.bac.br7.other.2)
length(wunifrac.dist.bac.br7.other$WUniFracdistance)
wunifrac.dist.bac.br7.other$Comparison <- "Br7_Others"
mean(wunifrac.dist.bac.br7.other$WUniFracdistance) #0.2971208

wunifrac.dist.among.br <- rbind(wunifrac.dist.bac.br2.other, wunifrac.dist.bac.br3.other,wunifrac.dist.bac.br4.other,
                            wunifrac.dist.bac.br5.other,wunifrac.dist.bac.br6.other,wunifrac.dist.bac.br7.other)
wunifrac.dist.among.br$Group <- "Among_branch"


wunifrac.dist.bac.br <- rbind(wunifrac.dist.within.br,wunifrac.dist.among.br)

### Statistical analysis
shapiro.test(wunifrac.dist.bac.br$WUniFracdistance) #p-value < 2.2e-16
var.test(wunifrac.dist.bac.br$WUniFracdistance~ wunifrac.dist.bac.br$Group)

#wilcox.test(wunifrac.dist.within.br$WUniFracdistance, wunifrac.dist.among.br$WUniFracdistance, alternative = "two.sided") #p-value = 2.05e-06
mean(wunifrac.dist.within.br$WUniFracdistance) #0.1748547
mean(wunifrac.dist.among.br$WUniFracdistance) #0.1865292

wunifrac.dist.bac.br$Group <- factor(wunifrac.dist.bac.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(wunifrac.dist.bac.br, aes(x=Group, y= WUniFracdistance, fill = Group))+
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

#### UniFrac dissimilarity - Fungi
wunifrac.dist.fun<-UniFrac(fun.clean.log, weighted = T, normalized = T, parallel = FALSE, fast = TRUE)
wunifrac.dist.fun <-as.matrix(wunifrac.dist.fun)

wunifrac.dist.fun.melt <- melt(as.matrix(wunifrac.dist.fun), na.rm = T)
head(wunifrac.dist.fun.melt)
names(wunifrac.dist.fun.melt)[3] <- "WUniFracdistance"
wunifrac.dist.fun.melt <- subset(wunifrac.dist.fun.melt, WUniFracdistance != 0)

#### Distance within each branch
Br2<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_2")]
Br3<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_3")]
Br4<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_4")]
Br5<-f.meta$SampleID[which(f.meta$Branch_number  == "Branch_5")]
Br6<-f.meta$SampleID[which( f.meta$Branch_number  == "Branch_6")]
Br7<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_7")]
Br8<-f.meta$SampleID[which( f.meta$Branch_number   == "Branch_8")]

## Within branch
wunifrac.dist.fun.br2 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br2 & Var2 %in% Br2)
length(wunifrac.dist.fun.br2$WUniFracdistance)
wunifrac.dist.fun.br2$Comparison <- "Branch_2"
mean(wunifrac.dist.fun.br2$WUniFracdistance) #0.0227567

wunifrac.dist.fun.br3 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br3 & Var2 %in% Br3)
length(wunifrac.dist.fun.br3$WUniFracdistance)
wunifrac.dist.fun.br3$Comparison <- "Branch_3"
mean(wunifrac.dist.fun.br3$WUniFracdistance) #0.01348914

wunifrac.dist.fun.br4 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br4 & Var2 %in% Br4)
length(wunifrac.dist.fun.br4$WUniFracdistance)
wunifrac.dist.fun.br4$Comparison <- "Branch_4"
mean(wunifrac.dist.fun.br4$WUniFracdistance) #0.01787269

wunifrac.dist.fun.br5 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br5 & Var2 %in% Br5)
length(wunifrac.dist.fun.br5$WUniFracdistance)
wunifrac.dist.fun.br5$Comparison <- "Branch_5"
mean(wunifrac.dist.fun.br5$WUniFracdistance) # 0.0210473

wunifrac.dist.fun.br6 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br6 & Var2 %in% Br6)
length(wunifrac.dist.fun.br6$WUniFracdistance)
wunifrac.dist.fun.br6$Comparison <- "Branch_6"
mean(wunifrac.dist.fun.br6$WUniFracdistance) #0.02785231

wunifrac.dist.fun.br7 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br7 & Var2 %in% Br7)
length(wunifrac.dist.fun.br7$WUniFracdistance)
wunifrac.dist.fun.br7$Comparison <- "Branch_7"
mean(wunifrac.dist.fun.br7$WUniFracdistance) # 0.01559169

wunifrac.dist.fun.br8 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br8 & Var2 %in% Br8)
length(wunifrac.dist.fun.br8$WUniFracdistance)
wunifrac.dist.fun.br8$Comparison <- "Branch_8"
mean(wunifrac.dist.fun.br8$WUniFracdistance) #0.0233345

wunifrac.dist.within.br.f <- rbind(wunifrac.dist.fun.br2, wunifrac.dist.fun.br3,wunifrac.dist.fun.br4,
                               wunifrac.dist.fun.br5,wunifrac.dist.fun.br6,wunifrac.dist.fun.br7,
                               wunifrac.dist.fun.br8)
wunifrac.dist.within.br.f$Group <- "Within_branch"

##Among branch
wunifrac.dist.fun.br2.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br2 & Var2 %in% c(Br3,Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.fun.br2.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br3,Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br2)
wunifrac.dist.fun.br2.other <- rbind(wunifrac.dist.fun.br2.other.1, wunifrac.dist.fun.br2.other.2)
length(wunifrac.dist.fun.br2.other$WUniFracdistance)
wunifrac.dist.fun.br2.other$Comparison <- "Br2_Others"
mean(wunifrac.dist.fun.br2.other$WUniFracdistance) #0.0225992

wunifrac.dist.fun.br3.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br3 & Var2 %in% c(Br4,Br5,Br6,Br7,Br8))
wunifrac.dist.fun.br3.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br4,Br5,Br6,Br7,Br8) & Var2 %in% Br3)
wunifrac.dist.fun.br3.other <- rbind(wunifrac.dist.fun.br3.other.1, wunifrac.dist.fun.br3.other.2)
length(wunifrac.dist.fun.br3.other$WUniFracdistance)
wunifrac.dist.fun.br3.other$Comparison <- "Br3_Others"
mean(wunifrac.dist.fun.br3.other$WUniFracdistance) #0.01921879

wunifrac.dist.fun.br4.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br4 & Var2 %in% c(Br5,Br6,Br7,Br8))
wunifrac.dist.fun.br4.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br5,Br6,Br7,Br8) & Var2 %in% Br4)
wunifrac.dist.fun.br4.other <- rbind(wunifrac.dist.fun.br4.other.1, wunifrac.dist.fun.br4.other.2)
length(wunifrac.dist.fun.br4.other$WUniFracdistance)
wunifrac.dist.fun.br4.other$Comparison <- "Br4_Others"
mean(wunifrac.dist.fun.br4.other$WUniFracdistance) #0.02174544

wunifrac.dist.fun.br5.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br5 & Var2 %in% c(Br6,Br7,Br8))
wunifrac.dist.fun.br5.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br6,Br7,Br8) & Var2 %in% Br5)
wunifrac.dist.fun.br5.other <- rbind(wunifrac.dist.fun.br5.other.1, wunifrac.dist.fun.br5.other.2)
length(wunifrac.dist.fun.br5.other$WUniFracdistance)
wunifrac.dist.fun.br5.other$Comparison <- "Br5_Others"
mean(wunifrac.dist.fun.br5.other$WUniFracdistance) #0.02311173

wunifrac.dist.fun.br6.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br6 & Var2 %in% c(Br7,Br8))
wunifrac.dist.fun.br6.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br7,Br8) & Var2 %in% Br6)
wunifrac.dist.fun.br6.other <- rbind(wunifrac.dist.fun.br6.other.1, wunifrac.dist.fun.br6.other.2)
length(wunifrac.dist.fun.br6.other$WUniFracdistance)
wunifrac.dist.fun.br6.other$Comparison <- "Br6_Others"
mean(wunifrac.dist.fun.br6.other$WUniFracdistance) #0.02477315

wunifrac.dist.fun.br7.other.1 <- subset(wunifrac.dist.fun.melt, Var1 %in% Br7 & Var2 %in% c(Br8))
wunifrac.dist.fun.br7.other.2 <- subset(wunifrac.dist.fun.melt, Var1 %in% c(Br8) & Var2 %in% Br7)
wunifrac.dist.fun.br7.other <- rbind(wunifrac.dist.fun.br7.other.1, wunifrac.dist.fun.br7.other.2)
length(wunifrac.dist.fun.br7.other$WUniFracdistance)
wunifrac.dist.fun.br7.other$Comparison <- "Br7_Others"
mean(wunifrac.dist.fun.br7.other$WUniFracdistance) #0.02187585

wunifrac.dist.among.br.f <- rbind(wunifrac.dist.fun.br2.other, wunifrac.dist.fun.br3.other,wunifrac.dist.fun.br4.other,
                              wunifrac.dist.fun.br5.other,wunifrac.dist.fun.br6.other,wunifrac.dist.fun.br7.other)
wunifrac.dist.among.br.f$Group <- "Among_branch"


wunifrac.dist.fun.br <- rbind(wunifrac.dist.within.br.f,wunifrac.dist.among.br.f)

### Statistical analysis
shapiro.test(wunifrac.dist.fun.br$WUniFracdistance) # p-value < 2.2e-16
var.test(wunifrac.dist.fun.br$WUniFracdistance~ wunifrac.dist.fun.br$Group)

#wilcox.test(wunifrac.dist.within.br.f$WUniFracdistance, wunifrac.dist.among.br.f$WUniFracdistance, alternative = "two.sided") #p-value = 2.05e-06
mean(wunifrac.dist.within.br.f$WUniFracdistance) # 0.02027776
mean(wunifrac.dist.among.br.f$WUniFracdistance) #0.02187754

wunifrac.dist.fun.br$Group <- factor(wunifrac.dist.fun.br$Group, levels = c("Within_branch","Among_branch"))

ggplot(wunifrac.dist.fun.br, aes(x=Group, y= WUniFracdistance, fill = Group))+
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
