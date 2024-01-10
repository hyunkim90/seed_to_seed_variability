##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
##Unconstrained PCoA
color_branch <- c("Branch_2" = "darkslategray4","Branch_3" = "dodgerblue4","Branch_4" = "steelblue3",
                  "Branch_5" = "salmon3", "Branch_6" = "tan1","Branch_7" = "lavenderblush4", "Branch_8" = "khaki3")

bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray')
bray2.bac <-  ordinate(bac.clean.log, 'CAP', 'bray',~ Branch_number)

# write.csv(bray1.bac$vectors, "Bacteria_PCoA_2017.csv")
# write.csv(bray2.bac$vectors, "Bacteria_PCoA_2018.csv")


plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='Branch_number', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  stat_ellipse(type = "t", linetype = 2, level = 0.95)

plot_ordination(bac.clean.log, bray2.bac, type = "samples", color='Branch_number', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+stat_ellipse(type = "t", linetype = 2, level = 0.95)

##### Fungal community#####
bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray')
bray2.fun <-  ordinate(fun.clean.log, 'CAP', 'bray',~Branch_number)

# write.csv(bray1.fun$vectors, "fungi_PCoA_2017.csv")
# write.csv(bray2.fun$vectors, "fungi_PCoA_2018.csv")

plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='Branch_number', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  stat_ellipse(type = "t", linetype = 2, level = 0.95)

plot_ordination(fun.clean.log, bray2.fun, type = "samples", color='Branch_number', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  stat_ellipse(type = "t", linetype = 2, level = 0.95)



##PERMANOVA
f.otu <- otu_table(fun.clean.log)
f.meta <- data.frame(sample_data(fun.clean.log))

f.permanova <- adonis(formula = t(f.otu) ~ (Branch_number), data = f.meta, permutations=9999, method = "bray")
f.permanova

b.otu <- otu_table(bac.clean.log)
b.meta <- data.frame(sample_data(bac.clean.log))

b.permanova <- adonis(formula = t(b.otu) ~ (Branch_number), data = b.meta, permutations=9999, method = "bray")
b.permanova






### Time series (DNA samples in which 3 seeds were pooled)
bac.time.clean.log
sample_data(bac.time.clean.log)
bray1.bac.time <-  ordinate(bac.time.clean.log, 'PCoA', 'bray')
plot_ordination(bac.time.clean.log, bray1.bac.time, type = "samples", color='Replication', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


bray1.fun.time <-  ordinate(fun.time.clean.log, 'PCoA', 'bray')
plot_ordination(fun.time.clean.log, bray1.fun.time, type = "samples", color='Replication', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


##PERMANOVA
f.otu <- otu_table(fun.time.clean.log)
f.meta <- data.frame(sample_data(fun.time.clean.log))
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication*Plot), data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication*Plot), data = f.meta, permutations=9999, method = "jaccard")
f.permanova

b.otu <- otu_table(bac.time.clean.log)
b.meta <- data.frame(sample_data(bac.time.clean.log))
## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication*Plot), data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication*Plot), data = b.meta, permutations=9999, method = "jaccard")
b.permanova

### 2018 samples after harvest (grouped by each plant)
bac.time.log.18.harvest <- subset_samples(bac.time.clean.log, Year == "year_2018")
bac.time.log.18.harvest <- subset_samples(bac.time.log.18.harvest, Age == "141days")

bray2.bac.time <-  ordinate(bac.time.log.18.harvest, 'PCoA', 'bray')
plot_ordination(bac.time.log.18.harvest, bray2.bac.time, type = "samples", color='Plot', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


fun.time.log.18.harvest <- subset_samples(fun.time.clean.log, Year == "year_2018")
fun.time.log.18.harvest <- subset_samples(fun.time.log.18.harvest, Age == "141days")

bray2.fun.time <-  ordinate(fun.time.log.18.harvest, 'PCoA', 'bray')
plot_ordination(fun.time.log.18.harvest, bray2.fun.time, type = "samples", color='Plot', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


### 2017 samples after harvest (grouped by each plant)
bac.time.log.17.harvest <- subset_samples(bac.time.clean.log, Year == "year_2017")
bac.time.log.17.harvest <- subset_samples(bac.time.log.17.harvest, Age == "140days")

bray2.bac.time <-  ordinate(bac.time.log.17.harvest, 'PCoA', 'bray')
plot_ordination(bac.time.log.17.harvest, bray2.bac.time, type = "samples", color='Plot', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


fun.time.log.17.harvest <- subset_samples(fun.time.clean.log, Year == "year_2017")
fun.time.log.17.harvest <- subset_samples(fun.time.log.17.harvest, Age == "140days")

bray2.fun.time <-  ordinate(fun.time.log.17.harvest, 'PCoA', 'bray')
plot_ordination(fun.time.log.17.harvest, bray2.fun.time, type = "samples", color='Plot', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ #scale_color_manual(values=color_branch)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## Perform pairwise comparisons
library(RVAideMemoire)
fun.time.meta <- sample_data(fun.time.log.18.harvest)
fun.time.dist.all <- vegdist(t(otu_table(fun.time.log.18.harvest)), method="bray")

bac.time.meta <- sample_data(bac.time.log.18.harvest)
bac.time.dist.all <- vegdist(t(otu_table(bac.time.log.18.harvest)), method="bray")


## Perform BETADISP test of multivariate dispersion
bdi_fun.time.dist.all <- betadisper(fun.time.dist.all,fun.time.meta$Plot,type = "centroid")
permutest(bdi_fun.time.dist.all,pairwise=T,permutations=how(nperm=9999))
boxplot(bdi_fun.time.dist.all)


bdi_bac.time.dist.all <- betadisper(bac.time.dist.all,bac.time.meta$Plot,type = "centroid")
permutest(bdi_bac.time.dist.all,pairwise=T,permutations=how(nperm=9999))
boxplot(bdi_bac.time.dist.all)


####
capscale()

f.otu <- otu_table(fun.clean.log)
f.meta <- data.frame(sample_data(fun.clean.log))

f.capscale <- capscale(formula = t(f.otu) ~ (Branch_number), data = f.meta,  dist = "bray")
summary(f.capscale)

b.otu <- otu_table(bac.clean.log)
b.meta <- data.frame(sample_data(bac.clean.log))

b.capscale <- capscale(formula = t(b.otu) ~ (Branch_number), data = b.meta,  dist = "bray")
b.capscale

