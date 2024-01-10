#### Core and Pan seed microbiota
bac.clean.ss.f
fun.clean.ss.f


### Abundance to prevalence table
otu.G.0 <- otu_table(bac.clean.ss.f)
df.otu.G.0 <- data.frame(otu.G.0)
head(df.otu.G.0)

df.otu.G.0[df.otu.G.0>0] <-1

df.otu.sum<-data.frame(rowSums(df.otu.G.0))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$OTU <- rownames(df.otu.sum)
bac.prev.tab <- merge(bac.list, df.otu.sum, by = c("OTU" = "OTU"))

bac.prev.tab<- bac.prev.tab %>% arrange(desc(Prevalence))


otu.G.0 <- otu_table(fun.clean.ss.f)
df.otu.G.0 <- data.frame(otu.G.0)
head(df.otu.G.0)

df.otu.G.0[df.otu.G.0>0] <-1

df.otu.sum<-data.frame(rowSums(df.otu.G.0))
names(df.otu.sum)[1] <- "sumPresence"
df.otu.sum$Prevalence <- df.otu.sum$sumPresence/70
df.otu.sum$OTU <- rownames(df.otu.sum)
fun.prev.tab <- merge(fun.list, df.otu.sum, by = c("OTU" = "OTU"))

fun.prev.tab<- fun.prev.tab %>% arrange(desc(Prevalence))


### Mean relative abundance
bac.clean.rel <- transform(bac.clean.ss.f, transform = "compositional")
melt.bac.rel <- psmelt(bac.clean.rel)
melt.bac.rel.mean <- melt.bac.rel %>% group_by(OTU) %>% summarise(MeanAbund = mean(Abundance), SdAbund = sd(Abundance))
melt.bac.rel.mean$logMeanAbund <- log10(melt.bac.rel.mean$MeanAbund)

bac.prev.tab <- merge(bac.prev.tab, melt.bac.rel.mean, by = c("OTU" = "OTU"))

bac.prev.tab<- bac.prev.tab %>% arrange(desc(Prevalence))


fun.clean.rel <- transform(fun.clean.ss.f, transform = "compositional")
melt.fun.rel <- psmelt(fun.clean.rel)
melt.fun.rel.mean <- melt.fun.rel %>% group_by(OTU) %>% summarise(MeanAbund = mean(Abundance), SdAbund = sd(Abundance))
melt.fun.rel.mean$logMeanAbund <- log10(melt.fun.rel.mean$MeanAbund)

fun.prev.tab <- merge(fun.prev.tab, melt.fun.rel.mean, by = c("OTU" = "OTU"))

fun.prev.tab<- fun.prev.tab %>% arrange(desc(Prevalence))

#write.csv(bac.prev.tab, "220727_prevalence and abundance table_bacteria.csv")
#write.csv(fun.prev.tab, "220727_prevalence and abundance table_fungi.csv")


####
bac.prev.tab <- read.csv("220727_prevalence and abundance table_bacteria_edit.csv")
fun.prev.tab <- read.csv("220727_prevalence and abundance table_fungi_edit.csv")

bac.prev.tab$Category <- factor(bac.prev.tab$Category, levels = c("Prevalent","Rare","Others"))
fun.prev.tab$Category <- factor(fun.prev.tab$Category, levels = c("Prevalent","Rare","Others"))
### Plot
ggplot(bac.prev.tab, aes(x=logMeanAbund, y=Prevalence, color = Category))+geom_point(size = 2.5, fill = "black")+
  geom_hline(yintercept = 0.2, linetype='dashed', color='black', size = 0.75)+
  geom_hline(yintercept = 0.8,linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept = log10(0.01),linetype='dashed', color='black', size = 0.75)+
  geom_errorbar(aes(x=logMeanAbund, y=Prevalence, xmin=logMeanAbund-log10(SdAbund), xmax = logMeanAbund+log10(SdAbund)))+
  mytheme_2d+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_x_continuous(breaks=seq(0,-12,-2))+
  theme(aspect.ratio = 1)
  
ggplot(fun.prev.tab, aes(x=logMeanAbund, y=Prevalence,color = Category))+geom_point(size = 2.5, fill = "black")+
  geom_hline(yintercept = 0.2, linetype='dashed', color='black', size = 0.75)+
  geom_hline(yintercept = 0.8,linetype='dashed', color='black', size = 0.75)+
  geom_vline(xintercept = log10(0.01),linetype='dashed', color='black', size = 0.75)+
  geom_errorbar(aes(x=logMeanAbund, y=Prevalence, xmin=logMeanAbund-log10(SdAbund), xmax = logMeanAbund+log10(SdAbund)))+
  mytheme_2d+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_x_continuous(breaks=seq(0,-12,-2))+theme(aspect.ratio = 1)




#### Prevalence-based approach
otu.G.0 <- otu_table(bac.clean.ss)
df.otu.G.0 <- data.frame(otu.G.0)
head(df.otu.G.0)

Group <- vector()
Core_number <- vector()
Core_list <- vector()
for (i in (0:ncol(df.otu.G.0))){
  core.0.bac <- core_members(bac.clean.ss, detection = 0, prevalence = i/70, include.lowest = F)
  l1<-length(core.0.bac)
  l2 <- paste(unique(bac.list$OTU_id[which(bac.list$OTU %in% core.0.bac)]), collapse = " ")
  Group <- append(Group, i)
  Core_number <- append(Core_number, l1)
  Core_list <- append(Core_list, l2)
}
core_preval_table.bac <- cbind(Group, Core_number, Core_list)
core_preval_table.bac <- cbind(Group, Core_number)
core_preval_table.bac <- data.frame(core_preval_table.bac)

write.csv(core_preval_table.bac, "Pan-core bac.csv")


otu.G.0 <- otu_table(bac.clean.ss)
df.otu.G.0 <- data.frame(otu.G.0)
head(df.otu.G.0)

Group <- vector()
Core_number <- vector()
Core_list <- vector()
for (i in (0:ncol(df.otu.G.0))){
  core.0.fun <- core_members(fun.clean.ss, detection = 0, prevalence = i/70, include.lowest = F)
  l1<-length(core.0.fun)
  l2 <- paste(unique(fun.list$OTU_id[which(fun.list$OTU %in% core.0.fun)]), collapse = " ")
  Group <- append(Group, i)
  Core_number <- append(Core_number, l1)
  Core_list <- append(Core_list, l2)
}
core_preval_table.fun <- cbind(Group, Core_number, Core_list)
core_preval_table.fun <- cbind(Group, Core_number)
core_preval_table.fun <- data.frame(core_preval_table.fun)

write.csv(core_preval_table.fun, "Pan-core fun.csv")


fun.core.pre<-read.csv("Pan-core fun.csv")
bac.core.pre<-read.csv("Pan-core bac.csv")

ggplot(bac.core.pre, aes(x=Group, y=Core_number))+geom_point(size = 1, color = "black")+
  geom_smooth(method="gam", span=0.9)+
ggplot(fun.core.pre, aes(x=Group, y=Core_number))+geom_point(size = 1, color = "black")+
  geom_smooth(method="gam", span=0.9)

dev.off()

#### 
bac.clean.ss.f


list.seed <- c("R1", "R2", "R3", "R4", "R5", "R6","R7","R8","R9",
               "R10", "R11", "R12","R13","R14","R15","R16", "R17", "R18","R19","R20","R21",
               "R22", "R23", "R24","R25","R26","R27","R28","R29","R30","R31","R32","R33",
               "R34","R35","R36","R37","R38","R39","R40","R41","R42","R43","R44","R45","R46",
               "R47","R48","R49","R50","R51","R52","R53","R54","R55","R56","R57","R58","R59","R60",
               "R61","R62","R63","R64","R65","R66","R67","R68","R69","R70")
sample_data(bac.clean.ss.f)

bac.clean.rel

for (i in list.seed){
  bac.clean.ss.f.wo1 <- subset_samples(bac.clean.rel, SampleID != i)
  core.bac.1 <- core_members(bac.clean.ss.f.wo1, detection = 0, prevalence = 1, include.lowest = F)
  print(length(core.fun.1))
  
}

for (i in list.seed){
  fun.clean.ss.f.wo1 <- subset_samples(fun.clean.rel, SampleID != i)
  core.fun.1 <- core_members(fun.clean.rel, detection = 0, prevalence = 0.9, include.lowest = F)
  print(length(core.fun.1))
  
}

### Core heatmaps
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(bac.clean.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()


# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

#Deletes "ASV" from taxa_names, e.g. ASV1 --> 1
#taxa_names(ps.m3.rel) = taxa_names(ps.m3.rel) %>% str_replace("ASV", "")
# Also define gray color palette

gray <- gray(seq(0,1,length=5))

p1 <- plot_core(bac.clean.rel,
                plot.type = "heatmap",
                colours = gray,
                prevalences = c(0, prevalences),
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

library(viridis)
print(p1 + scale_fill_viridis())



prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

#Deletes "ASV" from taxa_names, e.g. ASV1 --> 1
#taxa_names(ps.m3.rel) = taxa_names(ps.m3.rel) %>% str_replace("ASV", "")
# Also define gray color palette

gray <- gray(seq(0,1,length=5))

p1 <- plot_core(fun.clean.rel,
                plot.type = "heatmap",
                colours = gray,
                prevalences = prevalences,
                detections = c(0,detections), min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

library(viridis)
print(p1 + scale_fill_viridis())


### Core microbiota predicted using COREMIC (http://core-mic.com/) ### Not working...
## Reference: COREMIC: a web-tool to search for a niche associated CORE MICrobiome, 2018, PeerJ
### Preparing input files
# Data file
## Abundance and taxonomy table
# bac.asv.tab <- read.table("asv_table_final_bac.txt", sep = "\t", header = T)
# bac.asv.tab.2 <- subset(bac.asv.tab, OTU_id %in% rownames(otu_table(bac.clean.ss.f)))
# 
# for(i in as.character(bac.asv.tab.2$OTU_id)){
#     bac.asv.tab.2$OTU_id[which(bac.asv.tab.2$OTU_id == i)] <- bac.list$number[which(bac.list$OTU == i)]
#   }
# 
# write.table(bac.asv.tab.2, "asv_table_final_bac_edit.tsv", sep = "\t", quote = F)
# 
# fun.asv.tab <- read.table("asv_table_final_fun.txt", sep = "\t", header = T)
# fun.asv.tab.2 <- subset(fun.asv.tab, OTU_id %in% rownames(otu_table(fun.clean.ss.f)))
# for(i in as.character(fun.asv.tab.2$OTU_id)){
#   fun.asv.tab.2$OTU_id[which(fun.asv.tab.2$OTU_id == i)] <- fun.list$number[which(fun.list$OTU == i)]
# }
# 
# write.table(fun.asv.tab.2, "asv_table_final_fun_edit.tsv", sep = "\t", quote = F)
# 
# ### Covert table to Biom (in qiime 2)
# #biom convert -i asv_table_final_bac_edit.txt -o bac_COREMIC_datafile.biom --table-type="OTU table" --to-json
# #biom convert -i asv_table_final_fun_edit.txt -o fun_COREMIC_datafile.biom --table-type="OTU table" --to-json
# 
# ##Metadata (Group file)
# b.meta <- sample_data(bac.clean.ss.f)
# f.meta <- sample_data(fun.clean.ss.f)
# 
# write.table(b.meta, "sample data_bacteria.tsv", sep = "\t", quote = F)
# write.table(f.meta, "sample data_fungi.tsv", sep = "\t", quote = F)

core.bac <- core_members(bac.clean.ss, detection = 0, prevalence = 0.79, include.lowest = F)
core.fun <- core_members(fun.clean.ss, detection = 0, prevalence = 0.79, include.lowest = F)

fun.list$OTU_id[which(fun.list$OTU %in% core.fun)]




### Rarefaction curve: pan and core microbiota in seed

