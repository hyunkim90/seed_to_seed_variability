#####################################################################
#Original scripts adapted from Lee et al 2017 & Tucker et al 2016####
#R scripts modified by M. Amine Hassani - hassani.medamine@gmail.com#
#####################################################################
install.packages("reldist")
install.packages("bipartite")
#required R packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "PMCMR", "reldist","vegan","bipartite","phangorn","metagenomeSeq")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

#required source scripts
source("C:/Users/Hyun Kim/Desktop/SNU/Endophytes/Analysis/Analysis with 2018 samples/MetacommunityDynamicsFctsOikos.R") #Tucker et al 2016
source("C:/Users/Hyun Kim/Desktop/SNU/Endophytes/Analysis/Analysis with 2018 samples/PANullDevFctsOikos.R") #Tucker et al 2016
set.seed(131)

#color_palette<-c("#cc2a36","#166590","#11902d","#999999")
color_palette<-color_compartment
theme_new <- theme (
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="none",
  axis.text.x = element_text(angle=90, vjust=1)
)
sample_data(bac.clean.log)
# uploading and preparing phyloseq objects
## Branch
Bnull<-list()
Bnull.out<-data.frame()
LIST<-list ("Branch_2", "Branch_3", "Branch_4", "Branch_5","Branch_6", "Branch_7","Branch_8")

physeq_B<-bac.clean.log
sd<-sample_data(bac.clean.log)
sd <- data.frame(sd)
tax<-tax_table(bac.clean.log)

###############################################
#This script took 06:21:17 to run in HPCluster#
###############################################

for( i in LIST ) 
{
  print(paste0("Computing null deviation for ", i, " samples " ))
  #	print( paste0("subset ", i, " samples" ))
  SUB=c( i )
  physeq_SUB=subset_samples( physeq_B, Branch_number %in% SUB ) #Compartment or Microhabitat
  comm=data.frame(otu_table(physeq_SUB))
  
  map=sd[colnames(comm),]
  rdp=tax[rownames(comm),]
  
  list.OTUs=row.names(comm)[rowSums(comm)>1]
  COMM=comm[rownames(comm) %in% list.OTUs, ]
  
  #	print("preparing the data and setting parameters ... ")
  comm.t=t(COMM)
  bbs.sp.site <- comm.t
  patches=nrow(bbs.sp.site)
  rand <- 1000
  
  null.alphas <- matrix(NA, ncol(comm.t), rand)
  null.alpha <- matrix(NA, ncol(comm.t), rand)
  expected_beta <- matrix(NA, 1, rand)
  null.gamma <- matrix(NA, 1, rand)
  null.alpha.comp <- numeric()
  bucket_bray_res <- matrix(NA, patches, rand)
  
  bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
  mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
  gamma <- ncol(bbs.sp.site) #gamma
  obs_beta <- 1-mean.alpha/gamma
  obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma
  
  #	print (" Generating null patches ...")
  for (randomize in 1:rand) 
  {  
    null.dist = comm.t
    for (species in 1:ncol(null.dist)) 
    {
      tot.abund = sum(null.dist[,species])
      null.dist[,species] = 0
      for (individual in 1:tot.abund) 
      {
        sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
        null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
      }
    }
    #		print ("Calculating null deviation for null patches ...")
    null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
    null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
    expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
    null.alpha <- mean(null.alphas[,randomize])
    null.alpha.comp <- c(null.alpha.comp, null.alpha)
    
    bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
    diag(bucket_bray) <- NA
    bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
  }
  ## Computing beta-diversity for observed communities
  beta_comm_abund <- vegdist(comm.t, "bray")
  res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
  diag(res_beta_comm_abund) <- NA
  # output beta diversity (Bray)
  beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
  # output abundance beta-null deviation
  bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)
  
  betanull.out=data.frame(I(beta_div_abund_stoch),I(bray_abund_null_dev),stringsAsFactors=FALSE)
  colnames(betanull.out)=c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation")
  
  Bnull[[i]]<-betanull.out
  Bnull.out<-rbind(Bnull.out, betanull.out)
  
  print( paste0("Computing null deviation for ", i, " samples done" ))
  print("##############################################################")
}
print("Computing BC null deviations completed")

write.csv(Bnull.out, "bnullout_bacteria_seed.csv") #Branch

# upload BC null deviations of bacteria
Bnull=read.csv("bnullout_bacteria_seed.csv", row.names=1, header=T)
DT.B=data.table(Bnull, keep.rownames=T, key="rn")
setnames(DT.B, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCF","NullB"))
sd <- data.frame(sd)
DT.sd=data.table(sd, keep.rownames=T, key="rn")
DT=DT.sd[DT.B]
DT<-DT[,-c("BCF")]
DT.m.b<-melt(DT, measure.vars =c("NullB"))
DT.m.b$variable<-as.factor(DT.m.b$variable)
DT.m.b$VAR<-paste0(DT.m.b$Concat,DT.m.b$variable)
DT.m.b$kingdom<-ifelse(DT.m.b$variable=="NullB","B","NA")

# preparing the plot
#Branch
pp2 <- ggplot(data=DT.m[DT.m$kingdom=="B",], aes(x=Branch_number, y=as.numeric(value)))
panel_b=pp2+geom_boxplot(aes(fill=Branch_number)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="Abundance Null Deviation")+theme_bw(base_size=10)+
  scale_fill_jco() + theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_b
#	pdf("Fig3b.pdf",paper="A4" ,useDingbats=FALSE)
gridExtra::grid.arrange(panel_b, nrow=2, ncol=2)
#	dev.off()


###### Fungi ######
## Compartment
Bnull<-list()
Bnull.out<-data.frame()
LIST<-list ("Branch_2", "Branch_3", "Branch_4", "Branch_5","Branch_6", "Branch_7","Branch_8")

physeq_F<-fun.clean.log
sd<-sample_data(fun.clean.log)
sd <- data.frame(sd)
tax<-tax_table(fun.clean.log)


###############################################
#This script took 06:21:17 to run in HPCluster#
###############################################

for( i in LIST ) 
{
  print(paste0("Computing null deviation for ", i, " samples " ))
  #	print( paste0("subset ", i, " samples" ))
  SUB=c( i )
  physeq_SUB=subset_samples( physeq_F, Branch_number %in% SUB ) #Compartment or Microhabitat
  comm=data.frame(otu_table(physeq_SUB))
  
  map=sd[colnames(comm),]
  rdp=tax[rownames(comm),]
  
  list.OTUs=row.names(comm)[rowSums(comm)>1]
  COMM=comm[rownames(comm) %in% list.OTUs, ]
  
  #	print("preparing the data and setting parameters ... ")
  comm.t=t(COMM)
  bbs.sp.site <- comm.t
  patches=nrow(bbs.sp.site)
  rand <- 1000
  
  null.alphas <- matrix(NA, ncol(comm.t), rand)
  null.alpha <- matrix(NA, ncol(comm.t), rand)
  expected_beta <- matrix(NA, 1, rand)
  null.gamma <- matrix(NA, 1, rand)
  null.alpha.comp <- numeric()
  bucket_bray_res <- matrix(NA, patches, rand)
  
  bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
  mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
  gamma <- ncol(bbs.sp.site) #gamma
  obs_beta <- 1-mean.alpha/gamma
  obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma
  
  #	print (" Generating null patches ...")
  for (randomize in 1:rand) 
  {  
    null.dist = comm.t
    for (species in 1:ncol(null.dist)) 
    {
      tot.abund = sum(null.dist[,species])
      null.dist[,species] = 0
      for (individual in 1:tot.abund) 
      {
        sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
        null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
      }
    }
    #		print ("Calculating null deviation for null patches ...")
    null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
    null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
    expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
    null.alpha <- mean(null.alphas[,randomize])
    null.alpha.comp <- c(null.alpha.comp, null.alpha)
    
    bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
    diag(bucket_bray) <- NA
    bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
  }
  ## Computing beta-diversity for observed communities
  beta_comm_abund <- vegdist(comm.t, "bray")
  res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
  diag(res_beta_comm_abund) <- NA
  # output beta diversity (Bray)
  beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
  # output abundance beta-null deviation
  bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)
  
  betanull.out=data.frame(I(beta_div_abund_stoch),I(bray_abund_null_dev),stringsAsFactors=FALSE)
  colnames(betanull.out)=c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation")
  
  Bnull[[i]]<-betanull.out
  Bnull.out<-rbind(Bnull.out, betanull.out)
  
  print( paste0("Computing null deviation for ", i, " samples done" ))
  print("##############################################################")
}
print("Computing BC null deviations completed")

write.csv(Bnull.out, "bnullout_fungi_seed.csv")

# upload BC null deviations of bacteria
Bnull=read.csv("bnullout_fungi_seed.csv", row.names=1, header=T)
DT.B=data.table(Bnull, keep.rownames=T, key="rn")
setnames(DT.B, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCF","NullB"))
sd <- data.frame(sd)
DT.sd=data.table(sd, keep.rownames=T, key="rn")
DT=DT.sd[DT.B]
DT<-DT[,-c("BCF")]
DT.m.f<-melt(DT, measure.vars =c("NullB"))
DT.m.f$variable<-as.factor(DT.m.f$variable)
DT.m.f$VAR<-paste0(DT.m.f$Concat,DT.m.f$variable)
DT.m.f$kingdom<-ifelse(DT.m.f$variable=="NullB","F","NA")

DT.m <- rbind(DT.m.b, DT.m.f)

# preparing the plot
#Compartment
pp2 <- ggplot(data=DT.m[DT.m$kingdom=="F",], aes(x=Branch_number, y=as.numeric(value)))
panel_b=pp2+geom_boxplot(aes(fill=Branch_number)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="Abundance Null Deviation")+theme_bw(base_size=10)+
  scale_fill_jco() +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_b
#	pdf("Fig3b.pdf",paper="A4" ,useDingbats=FALSE)
gridExtra::grid.arrange(panel_b, nrow=2, ncol=2)
dev.off()

write.csv(DT.m, "Beta null deviation_bacteria and fungi_log.csv")
shapiro.test(DT.m$value)
pp2 <- ggplot(data=DT.m, aes(x=Branch_number, y=as.numeric(value)))
panel_b=pp2+geom_boxplot(aes(fill=kingdom)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="Abundance Null Deviation")+theme_bw(base_size=10)+
  scale_fill_jco() +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  ggsignif::geom_signif(comparisons = list(c("B","F")), test = "t.test",na.rm = F)
panel_b
#	pdf("Fig3b.pdf",paper="A4" ,useDingbats=FALSE)
gridExtra::grid.arrange(panel_b, nrow=2, ncol=2)
dev.off()



### Comparison bacteria and fungi at each seed level
### Bacteria
Bnull<-list()
Bnull.out<-data.frame()
physeq_B<-bac.clean.log
sd<-sample_data(bac.clean.log)
sd <- data.frame(sd)
tax<-tax_table(bac.clean.log)

physeq_SUB <- bac.clean.log


### Fungi
Bnull<-list()
Bnull.out<-data.frame()
physeq_B<-fun.clean.log
sd<-sample_data(fun.clean.log)
sd <- data.frame(sd)
tax<-tax_table(fun.clean.log)

physeq_SUB <- fun.clean.log


## compute beta null deviation
comm=data.frame(otu_table(physeq_SUB))

map=sd[colnames(comm),]
rdp=tax[rownames(comm),]

list.OTUs=row.names(comm)[rowSums(comm)>1]
COMM=comm[rownames(comm) %in% list.OTUs, ]

#	print("preparing the data and setting parameters ... ")
comm.t=t(COMM)
bbs.sp.site <- comm.t
patches=nrow(bbs.sp.site)
rand <- 1000

null.alphas <- matrix(NA, ncol(comm.t), rand)
null.alpha <- matrix(NA, ncol(comm.t), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
gamma <- ncol(bbs.sp.site) #gamma
obs_beta <- 1-mean.alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

#	print (" Generating null patches ...")
for (randomize in 1:rand) 
{  
  null.dist = comm.t
  for (species in 1:ncol(null.dist)) 
  {
    tot.abund = sum(null.dist[,species])
    null.dist[,species] = 0
    for (individual in 1:tot.abund) 
    {
      sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }
  }
  #		print ("Calculating null deviation for null patches ...")
  null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
  null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
  expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
  null.alpha <- mean(null.alphas[,randomize])
  null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
  bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
  diag(bucket_bray) <- NA
  bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
}
## Computing beta-diversity for observed communities
beta_comm_abund <- vegdist(comm.t, "bray")
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
diag(res_beta_comm_abund) <- NA
# output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
# output abundance beta-null deviation
bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)

betanull.out=data.frame(I(beta_div_abund_stoch),I(bray_abund_null_dev),stringsAsFactors=FALSE)
colnames(betanull.out)=c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation")

Bnull[[i]]<-betanull.out
Bnull.out<-rbind(Bnull.out, betanull.out)

print("Done")

write.csv(Bnull.out, "Bull_bacteria_all seeds.csv")
write.csv(Bnull.out, "Bull_fungi_all seeds.csv")

Bnull=read.csv("Bull_bacteria_all seeds.csv", row.names=1, header=T)
DT.B=data.table(Bnull, keep.rownames=T, key="rn")
setnames(DT.B, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCF","NullB"))
sd <- data.frame(sd)
DT.sd=data.table(sd, keep.rownames=T, key="rn")
DT=DT.sd[DT.B]
DT<-DT[,-c("BCF")]
DT.m.1<-melt(DT, measure.vars =c("NullB"))
DT.m.1$variable<-as.factor(DT.m.1$variable)
DT.m.1$VAR<-paste0(DT.m.1$Concat,DT.m.1$variable)
DT.m.1$kingdom<-ifelse(DT.m.1$variable=="NullB","B","NA")


Bnull=read.csv("Bull_fungi_all seeds.csv", row.names=1, header=T)
DT.B=data.table(Bnull, keep.rownames=T, key="rn")
setnames(DT.B, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCF","NullB"))
sd <- data.frame(sd)
DT.sd=data.table(sd, keep.rownames=T, key="rn")
DT=DT.sd[DT.B]
DT<-DT[,-c("BCF")]
DT.m.2<-melt(DT, measure.vars =c("NullB"))
DT.m.2$variable<-as.factor(DT.m.2$variable)
DT.m.2$VAR<-paste0(DT.m.2$Concat,DT.m.2$variable)
DT.m.2$kingdom<-ifelse(DT.m.2$variable=="NullB","F","NA")


DT.m <- rbind(DT.m.1,DT.m.2)

write.csv(DT.m,"Raw table for plotting beta null deviation.csv")

shapiro.test(DT.m$value)

pp2 <- ggplot(data=DT.m, aes(x=kingdom, y=as.numeric(value)))
panel_b=pp2+geom_boxplot(aes(fill=kingdom)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="Abundance Null Deviation")+theme_bw(base_size=10)+
  scale_fill_jco() +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+theme(aspect.ratio = 2)+
  ggsignif::stat_signif(comparisons = list(c("B", "F")),test = "t.test")
panel_b
#	pdf("Fig3b.pdf",paper="A4" ,useDingbats=FALSE)
gridExtra::grid.arrange(panel_b, nrow=2, ncol=2)
dev.off()

