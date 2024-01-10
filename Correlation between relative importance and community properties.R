#####
bac.phylo.bin.DR <- subset(bac.phylo.bin.melt.mean, Index == "DR")
bac.phylo.bin.DL <- subset(bac.phylo.bin.melt.mean, Index == "DL")
bac.phylo.bin.HoS <- subset(bac.phylo.bin.melt.mean, Index == "HoS")
bac.phylo.bin.HeS <- subset(bac.phylo.bin.melt.mean, Index == "HeS")
bac.phylo.bin.HD <- subset(bac.phylo.bin.melt.mean, Index == "HD")


fun.phylo.bin.DR <- subset(fun.phylo.bin.melt.mean, Index == "DR")
fun.phylo.bin.DL <- subset(fun.phylo.bin.melt.mean, Index == "DL")
fun.phylo.bin.HoS <- subset(fun.phylo.bin.melt.mean, Index == "HoS")
fun.phylo.bin.HeS <- subset(fun.phylo.bin.melt.mean, Index == "HeS")
fun.phylo.bin.HD <- subset(fun.phylo.bin.melt.mean, Index == "HD")


### Abundance and prevalence table
### Bacteria
### Abundance
melt.bac.tab <- psmelt(bac.clean.ss.f)
head(melt.bac.tab)
melt.bac.tab$Bin <- 0
for (i in as.character(bacterial.bins$ID)){
  melt.bac.tab$Bin[which(melt.bac.tab$OTU == i)] <- bacterial.bins$Bin[which(bacterial.bins$ID == i)]
}
bin.abun <- melt.bac.tab %>% group_by(Bin) %>% summarise(sumAbund = sum(Abundance))
bin.abun$MeanAbund <- bin.abun$sumAbund/sum(bin.abun$sumAbund)

###Prevalence
bac.prev.tab
head(bac.prev.tab)
names(bac.prev.tab)[2] <- "ID"
bac.prev.bin.tab <- merge(bacterial.bins  , bac.prev.tab, by = c("ID" = "ID"))
bac.prev.bin.tab.2 <- bac.prev.bin.tab  %>% group_by(Bin) %>% summarise(meanPrev= mean(Prevalence))


##Merge abundance table and prevalence table
bac.abund.prev <- merge(bin.abun, bac.prev.bin.tab.2, by = "Bin")

###Fungi
###Abundance
melt.fun.tab <- psmelt(fun.clean.ss.f)
head(melt.fun.tab)
melt.fun.tab$Bin <- 0
for (i in as.character(fungal.bins$ID)){
  melt.fun.tab$Bin[which(melt.fun.tab$OTU == i)] <- fungal.bins$Bin[which(fungal.bins$ID == i)]
}
bin.abun.f <- melt.fun.tab %>% group_by(Bin) %>% summarise(sumAbund = sum(Abundance))
bin.abun.f$MeanAbund <- bin.abun.f$sumAbund/sum(bin.abun.f$sumAbund)


##Prevalence
fun.prev.tab
head(fun.prev.tab)
names(fun.prev.tab)[2] <- "ID"
fun.prev.bin.tab <- merge(fungal.bins  , fun.prev.tab, by = c("ID" = "ID"))
fun.prev.bin.tab.2 <- fun.prev.bin.tab  %>% group_by(Bin) %>% summarise(meanPrev= mean(Prevalence))

##Merge abundance table and prevalence table
fun.abund.prev <- merge(bin.abun.f, fun.prev.bin.tab.2, by = "Bin")


### Correlation with relative abundance
bac.phylo.bin.HoS.PCT <- merge(bac.phylo.bin.HoS, bac.abund.prev, by = "Bin")
bac.phylo.bin.HeS.PCT <- merge(bac.phylo.bin.HeS, bac.abund.prev, by = "Bin")
bac.phylo.bin.DL.PCT <- merge(bac.phylo.bin.DL, bac.abund.prev, by = "Bin")
bac.phylo.bin.HD.PCT <- merge(bac.phylo.bin.HD, bac.abund.prev, by = "Bin")
bac.phylo.bin.DR.PCT <- merge(bac.phylo.bin.DR, bac.abund.prev, by = "Bin")

b.HoS.RA<-ggplot(data=bac.phylo.bin.HoS.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HoS.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  0.579, P = 1.49e-05, estimate = 2.92696

b.HeS.RA<-ggplot(data=bac.phylo.bin.HeS.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Heterogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HeS.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  -0.00691, P = 0.36727, estimate = 0.21495

b.DL.RA<-ggplot(data=bac.phylo.bin.DL.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Dispersal limitation")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.DL.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  0.309, P = 0.00347, estimate = -3.01758


b.HD.RA<-ggplot(data=bac.phylo.bin.HD.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Homogenizing dispersal")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HD.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  0.0443, P = 0.170, estimate = 0.0250632


b.DR.RA<-ggplot(data=bac.phylo.bin.DR.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.DR.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  -0.0452, P = 0.829, estimate = -0.14938


### Correlation with prevalence
b.HoS.prev<-ggplot(data=bac.phylo.bin.HoS.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HoS.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.209, P = 0.0163, estimate = 3.53507

b.HeS.prev<-ggplot(data=bac.phylo.bin.HeS.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Heterogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HeS.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.00657, P = 0.297, estimate = 0.46804

b.DL.prev<-ggplot(data=bac.phylo.bin.DL.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Dispersal limitation")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.DL.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.384, P = 0.000957, estimate = -6.26676


b.HD.prev<-ggplot(data=bac.phylo.bin.HD.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Homogenizing dispersal")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HD.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.114, P = 0.0637, estimate = 0.0627366


b.DR.prev<-ggplot(data=bac.phylo.bin.DR.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.DR.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.0977, P = 0.0801, estimate = 2.20091



###Fungi
### Correlation with relative abundance
fun.phylo.bin.HoS.PCT <- merge(fun.phylo.bin.HoS, fun.abund.prev, by = "Bin")
fun.phylo.bin.HeS.PCT <- merge(fun.phylo.bin.HeS, fun.abund.prev, by = "Bin")
fun.phylo.bin.DL.PCT <- merge(fun.phylo.bin.DL, fun.abund.prev, by = "Bin")
fun.phylo.bin.HD.PCT <- merge(fun.phylo.bin.HD, fun.abund.prev, by = "Bin")
fun.phylo.bin.DR.PCT <- merge(fun.phylo.bin.DR, fun.abund.prev, by = "Bin")

f.HoS.RA<-ggplot(data=fun.phylo.bin.HoS.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HoS.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  -0.0511, P = 0.644, estimate = 0.23605

f.HeS.RA<-ggplot(data=fun.phylo.bin.HeS.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Heterogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HeS.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  -0.0154, P = 0.398, estimate = -0.07224

f.DL.RA<-ggplot(data=fun.phylo.bin.DL.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Dispersal limitation")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DL.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  0.453, P = 0.00183, estimate = -2.11924


f.HD.RA<-ggplot(data=fun.phylo.bin.HD.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Homogenizing dispersal")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HD.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  -0.0578, P = 0.7276, estimate = -0.006802


f.DR.RA<-ggplot(data=fun.phylo.bin.DR.PCT, aes(y=MeanRelImportance, x= MeanAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin abundance %")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DR.PCT, MeanRelImportance~MeanAbund)
summary(test2) # R-sq.(adj) =  0.427, P = 0.00264, estimate = 1.96223


### Correlation with prevalence
f.HoS.prev<-ggplot(data=fun.phylo.bin.HoS.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HoS.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  -0.0412, P =  0.554, estimate = 0.30210

f.HeS.prev<-ggplot(data=fun.phylo.bin.HeS.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Heterogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HeS.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  -0.0314, P = 0.485, estimate = -0.05995

f.DL.prev<-ggplot(data=fun.phylo.bin.DL.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Dispersal limitation")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DL.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.504, P = 0.000846, estimate = -2.22323


f.HD.prev<-ggplot(data=fun.phylo.bin.HD.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Homogenizing dispersal")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.HD.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.00439, P = 0.317, estimate = 0.019259


f.DR.prev<-ggplot(data=fun.phylo.bin.DR.PCT, aes(y=MeanRelImportance, x= meanPrev))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Bin prevalence %")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DR.PCT, MeanRelImportance~meanPrev)
summary(test2) # R-sq.(adj) =  0.426, P = 0.00269, estimate = 1.96182


### UniFrac distance
##Define functions
Wunifrac.bins<-function(list.bins, physeqs, keywords){
  bac.clean.nolog.bin <- prune_taxa(subset(list.bins, Bin == keywords)$ID,physeqs)
  dist.unifrac.bin<-UniFrac(bac.clean.nolog.bin, weighted=T, normalized=TRUE, parallel=FALSE, fast=TRUE)
  dist.unifrac.bin<-as.matrix(dist.unifrac.bin)
  
  get_lower_tri<-function(cormat){
    cormat[base::upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  dist.test.bin.lower<-get_lower_tri(dist.unifrac.bin)
  dist.test.bin.melt <- melt(as.matrix(dist.test.bin.lower), na.rm = T)
  head(dist.test.bin.melt)
  names(dist.test.bin.melt)[3] <- "wUnifrac"
  dist.test.bin.melt <- subset(dist.test.bin.melt, wUnifrac != 0)
  meanValue<-mean(dist.test.bin.melt$wUnifrac)#0.6605315 #0.1672967
  return(meanValue)
}

#### Unweighted unifrac
UWunifrac.bins<-function(list.bins, physeqs, keywords){
  bac.clean.nolog.bin <- prune_taxa(subset(list.bins, Bin == keywords)$ID,physeqs)
  dist.unifrac.bin<-UniFrac(bac.clean.nolog.bin, weighted=F, normalized=TRUE, parallel=FALSE, fast=TRUE)
  dist.unifrac.bin<-as.matrix(dist.unifrac.bin)
  
  get_lower_tri<-function(cormat){
    cormat[base::upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  dist.test.bin.lower<-get_lower_tri(dist.unifrac.bin)
  dist.test.bin.melt <- melt(as.matrix(dist.test.bin.lower), na.rm = T)
  head(dist.test.bin.melt)
  names(dist.test.bin.melt)[3] <- "wUnifrac"
  dist.test.bin.melt <- subset(dist.test.bin.melt, wUnifrac != 0)
  meanValue<-mean(dist.test.bin.melt$wUnifrac)#0.6605315 #0.1672967
  return(meanValue)
}

### Unifrac distance of bacteria
Bin.list.bac<-unique(bacterial.bins$Bin)

df.unifrac.bac <- data.frame(matrix(ncol =3, nrow = 23))
names(df.unifrac.bac)[1] <- "Bin"
names(df.unifrac.bac)[2] <- "WUniFrac"
names(df.unifrac.bac)[3] <- "UWUniFrac"

df.unifrac.bac$Bin <- Bin.list.bac
for(i in as.character(Bin.list.bac)){
  Wuni.Value<-Wunifrac.bins(bacterial.bins, bac.clean.nolog, keywords = i)
  df.unifrac.bac$WUniFrac[which(df.unifrac.bac$Bin == i)] <- as.numeric(Wuni.Value)
}

for(i in as.character(Bin.list.bac)){
  UWuni.Value<-UWunifrac.bins(bacterial.bins, bac.clean.nolog, keywords = i)
  df.unifrac.bac$UWUniFrac[which(df.unifrac.bac$Bin == i)] <- as.numeric(UWuni.Value)
}

### Unifrac distance of fungi
Bin.list.fun<-unique(fungal.bins$Bin)

df.unifrac.fun <- data.frame(matrix(ncol =3, nrow = 17))
names(df.unifrac.fun)[1] <- "Bin"
names(df.unifrac.fun)[2] <- "WUniFrac"
names(df.unifrac.fun)[3] <- "UWUniFrac"

df.unifrac.fun$Bin <- Bin.list.fun
for(i in as.character(Bin.list.fun)){
  Wuni.Value<-Wunifrac.bins(fungal.bins, fun.clean.nolog, keywords = i)
  df.unifrac.fun$WUniFrac[which(df.unifrac.fun$Bin == i)] <- as.numeric(Wuni.Value)
}

for(i in as.character(Bin.list.fun)){
  UWuni.Value<-UWunifrac.bins(fungal.bins, fun.clean.nolog, keywords = i)
  df.unifrac.fun$UWUniFrac[which(df.unifrac.fun$Bin == i)] <- as.numeric(UWuni.Value)
}

###Merge tables
bac.phylo.bin.HoS.PCT <- merge(bac.phylo.bin.HoS.PCT, df.unifrac.bac, by = "Bin")
bac.phylo.bin.HeS.PCT <- merge(bac.phylo.bin.HeS.PCT, df.unifrac.bac, by = "Bin")
bac.phylo.bin.DL.PCT <- merge(bac.phylo.bin.DL.PCT, df.unifrac.bac, by = "Bin")
bac.phylo.bin.HD.PCT <- merge(bac.phylo.bin.HD.PCT, df.unifrac.bac, by = "Bin")
bac.phylo.bin.DR.PCT <- merge(bac.phylo.bin.DR.PCT, df.unifrac.bac, by = "Bin")


fun.phylo.bin.HoS.PCT <- merge(fun.phylo.bin.HoS.PCT, df.unifrac.fun, by = "Bin")
fun.phylo.bin.HeS.PCT <- merge(fun.phylo.bin.HeS.PCT, df.unifrac.fun, by = "Bin")
fun.phylo.bin.DL.PCT <- merge(fun.phylo.bin.DL.PCT, df.unifrac.fun, by = "Bin")
fun.phylo.bin.HD.PCT <- merge(fun.phylo.bin.HD.PCT, df.unifrac.fun, by = "Bin")
fun.phylo.bin.DR.PCT <- merge(fun.phylo.bin.DR.PCT, df.unifrac.fun, by = "Bin")

b.HoS.wuni<-ggplot(data=bac.phylo.bin.HoS.PCT, aes(y=MeanRelImportance, x= WUniFrac))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Weighted UniFrac distance")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  #scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HoS.PCT, MeanRelImportance~WUniFrac)
summary(test2) #R-sq.(adj) =  0.555, P = 2.72e-05, estimate =  -0.83463


b.HeS.wuni<-ggplot(data=bac.phylo.bin.HeS.PCT, aes(y=MeanRelImportance, x= WUniFrac))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Weighted UniFrac distance")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  #scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=bac.phylo.bin.HeS.PCT, MeanRelImportance~WUniFrac)
summary(test2) #R-sq.(adj) =  0.555, P = 2.72e-05, estimate =  -0.83463

f.DR.wuni<-ggplot(data=fun.phylo.bin.DR.PCT , aes(y=MeanRelImportance, x= WUniFrac))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Weighted UniFrac distance")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  #scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DR.PCT, MeanRelImportance~WUniFrac)
summary(test2) #R-sq.(adj) =  0.0801, P = 0.143, estimate = -0.18496



### Composite figure-Main

pdf("220518_Community properties and assembly mechanism_main figure.pdf", width = 20, height =10)
ggpubr::ggarrange(b.DL.RA,f.DL.RA,b.DL.prev,f.DL.prev,b.HoS.RA, f.DR.RA,b.HoS.wuni, f.DR.wuni,
          labels = c("a","b","c","d","e","f","g","h"), ncol = 4, nrow=2)

dev.off()


### Composite figure-supple
pdf("220502_Community properties and assembly mechanism_abundance_supple figure.pdf", width = 20, height =10)
ggpubr::ggarrange(b.HeS.RA,b.HD.RA,b.DR.RA,f.HoS.RA,f.HeS.RA,f.HD.RA,
                  labels = c("a","b","c","d","e","f"), ncol = 3, nrow=2)

dev.off()


pdf("220502_Community properties and assembly mechanism_prevalence_supple figure.pdf", width = 20, height =10)
ggpubr::ggarrange(b.HeS.prev,b.HD.prev,b.DR.prev,f.HoS.prev,f.HeS.prev,f.HD.prev,
                  labels = c("a","b","c","d","e","f"), ncol = 3, nrow=2)

dev.off()



### Abundance deviation
### Abundance
melt.bac.tab <- psmelt(bac.clean.rel)
head(melt.bac.tab)
melt.bac.tab$Bin <- 0
for (i in as.character(bacterial.bins$ID)){
  melt.bac.tab$Bin[which(melt.bac.tab$OTU == i)] <- bacterial.bins$Bin[which(bacterial.bins$ID == i)]
}
bin.abun.sd <- melt.bac.tab %>% group_by(Bin, OTU) %>% summarise(SdAbund = sd(Abundance))
bin.abun.sd.2 <- bin.abun.sd %>% group_by(Bin) %>% summarise(meanSdAbund = mean(SdAbund))

melt.fun.tab <- psmelt(fun.clean.rel)
head(melt.fun.tab)
melt.fun.tab$Bin <- 0
for (i in as.character(fungal.bins$ID)){
  melt.fun.tab$Bin[which(melt.fun.tab$OTU == i)] <- fungal.bins$Bin[which(fungal.bins$ID == i)]
}
bin.abun.sd.f <- melt.fun.tab %>% group_by(Bin, OTU) %>% summarise(SdAbund = sd(Abundance))
bin.abun.sd.f.2 <- bin.abun.sd.f %>% group_by(Bin) %>% summarise(meanSdAbund = mean(SdAbund))


bac.phylo.bin.HoS.PCT <- merge(bac.phylo.bin.HoS.PCT, bin.abun.sd.2, by = "Bin")
bac.phylo.bin.HeS.PCT <- merge(bac.phylo.bin.HeS.PCT, bin.abun.sd.2, by = "Bin")
bac.phylo.bin.DL.PCT <- merge(bac.phylo.bin.DL.PCT, bin.abun.sd.2, by = "Bin")
bac.phylo.bin.HD.PCT <- merge(bac.phylo.bin.HD.PCT, bin.abun.sd.2, by = "Bin")
bac.phylo.bin.DR.PCT <- merge(bac.phylo.bin.DR.PCT, bin.abun.sd.2, by = "Bin")


fun.phylo.bin.HoS.PCT <- merge(fun.phylo.bin.HoS.PCT, bin.abun.sd.f.2 , by = "Bin")
fun.phylo.bin.HeS.PCT <- merge(fun.phylo.bin.HeS.PCT, bin.abun.sd.f.2 , by = "Bin")
fun.phylo.bin.DL.PCT <- merge(fun.phylo.bin.DL.PCT, bin.abun.sd.f.2 , by = "Bin")
fun.phylo.bin.HD.PCT <- merge(fun.phylo.bin.HD.PCT, bin.abun.sd.f.2 , by = "Bin")
fun.phylo.bin.DR.PCT <- merge(fun.phylo.bin.DR.PCT, bin.abun.sd.f.2 , by = "Bin")


f.DR.wuni<-ggplot(data=fun.phylo.bin.DR.PCT , aes(y=MeanRelImportance, x= meanSdAbund))+
  geom_smooth(method="gam", span=0.9)+ geom_point() + theme_bw()+
  xlab("Abundance deviation")+ylab("Ecological drift")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  #scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-mgcv::gam(data=fun.phylo.bin.DR.PCT, MeanRelImportance~meanSdAbund)
summary(test2) #R-sq.(adj) =  0.0801, P = 0.143, estimate = -0.18496


b.HoS.dev<-ggplot(data=bac.phylo.bin.HoS.PCT , aes(y=MeanRelImportance, x= meanSdAbund))+
  geom_smooth(method="lm", span=0.9)+ geom_point() + theme_bw()+
  xlab("Abundance deviation")+ylab("Homogeneous selection")+
  #scale_color_manual(values = c("deepskyblue","#ff00ff", "gold", "#00ff00"))+ scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_y_continuous(labels = scales::percent, limits = c(-0.25, 1))+
  #scale_x_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"))+ theme(legend.position = "none")+
  theme(aspect.ratio = 1)

test2<-lm(data=bac.phylo.bin.HoS.PCT, MeanRelImportance~meanSdAbund)
summary(test2) #R-sq.(adj) =  0.0801, P = 0.143, estimate = -0.18496



###Supplementary Tables
write.csv(bac.phylo.bin.HoS.PCT, "220523_bac.phylo.bin.HoS.PCT.csv")
write.csv(bac.phylo.bin.HeS.PCT, "220523_bac.phylo.bin.HeS.PCT.csv")
write.csv(bac.phylo.bin.DL.PCT, "220523_bac.phylo.bin.DL.PCT.csv")
write.csv(bac.phylo.bin.HD.PCT, "220523_bac.phylo.bin.HD.PCT.csv")
write.csv(bac.phylo.bin.DR.PCT, "220523_bac.phylo.bin.DR.PCT.csv")

write.csv(fun.phylo.bin.HoS.PCT, "220523_fun.phylo.bin.HoS.PCT.csv")
write.csv(fun.phylo.bin.HeS.PCT, "220523_fun.phylo.bin.HeS.PCT.csv")
write.csv(fun.phylo.bin.DL.PCT, "220523_fun.phylo.bin.DL.PCT.csv")
write.csv(fun.phylo.bin.HD.PCT, "220523_fun.phylo.bin.HD.PCT.csv")
write.csv(fun.phylo.bin.DR.PCT, "220523_fun.phylo.bin.DR.PCT.csv")
