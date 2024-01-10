
###### 2. make barplot with phyloseq
### Phyloseqs
bac.clean.ss
fun.clean.ss

sum(otu_table(fun.clean.ss))


order.sample <- c("R1", "R2", "R3", "R4", "R5", "R6","R7","R8","R9",
                  "R10", "R11", "R12","R13","R14","R15","R16", "R17", "R18","R19","R20","R21",
                  "R22", "R23", "R24","R25","R26","R27","R28","R29","R30","R31","R32","R33",
                  "R34","R35","R36","R37","R38","R39","R40","R41","R42","R43","R44","R45","R46",
                  "R47","R48","R49","R50","R51","R52","R53","R54","R55","R56","R57","R58","R59","R60",
                  "R61","R62","R63","R64","R65","R66","R67","R68","R69","R70")


## let's get phylum sir
df.phylum <- bac.clean.ss %>%
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
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

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

## Branch

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(Branch_number) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

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
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Branch_number, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  theme(aspect.ratio = 1.5)+
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Deltaproteobacteria" = "darkolivegreen4","Proteobacteria" = "#003333","Actinobacteriota"="indianred2","Acidobacteriota"="lightsalmon4",
                                "Firmicutes" ="tan1","unidentified" = "black", "Low abundance" = "light grey")) +
  
  
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



#### Fungi
## Phylum level
##Grouped by sample
df.phylum.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()

df.phylum.fun %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
levels(df.phylum.fun$Phylum) = c(levels(df.phylum.fun$Phylum), 'Low abundance')

# we need to group by samples
df.phylum.fun.rel <- df.phylum.fun %>%  
  group_by(SampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.fun.rel[df.phylum.fun.rel$RelAbundance < 0.1,]$Phylum <- 'Low abundance'
unique(df.phylum.fun$Phylum)

ord.f <- df.phylum.fun.rel %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Phylum
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.phylum.fun.rel$Phylum <- factor(df.phylum.fun.rel$Phylum, levels = vec.reorder.f) 
df.phylum.fun.rel$SampleID <-factor(df.phylum.fun.rel$SampleID, levels = order.sample)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.fun.rel.p1 <- ggplot(df.phylum.fun.rel, aes(x=SampleID, y = RelAbundance, fill = Phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Ascomycota" = "#11335F", "Basidiomycota"= "#BE4146",
                               "unidentified" ="#000000",
                               "Low abundance"="#BFBEBE")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.phylum.fun.rel.p1


dev.off()

###class
df.class.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

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

### Branch
df.class.fun.rel <- df.class.fun %>%  
  group_by(Branch_number) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 0.01,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)

df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=Branch_number, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Dothideomycetes" = "#5195D1","Ustilaginomycetes"="#ffcc33","Sordariomycetes"= "#1E63AF",
                               "Tremellomycetes"="#BE4146","Cystobasidiomycetes" = "#A871AE","Microbotryomycetes" ="#DC9A9E",
                               "Agaricomycetes" = "#CC6C71",
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
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 1.5)
df.class.fun.rel.p1

dev.off()


write.csv(df.class.fun.rel,"Relative abundance table_fungi_class_211214_latest.csv")
write.csv(df.phylum.rel,"Relative abundance table_bacteria_class and phyla_211028.csv")


### Time series
order.sample.time <- c("R1", "R2", "R3", "R4", "R5", "R6","R7","R8","R9",
                  "R10", "R11", "R12","R13","R14","R15","R16", "R17", "R18","R19","R20","R21",
                  "R22", "R23", "R24","R25","R26","R27","R28","R29","R30","R31","R32","R33",
                  "R34","R35","R36","R37","R38","R39","R40","R41","R42","R43","R44","R45","R46",
                  "R47","R48","R49","R50","R51","R52","R53","R54","R55","R56","R57","R58","R59","R60",
                  "R61","R62","R63","R64","R65","R66","R67","R68","R69","R70")
bac.time.meta <- read.table(file = './Time series data/bacteria/sample_metadata.tsv', sep = '\t', header = TRUE)
rownames(bac.time.meta) <- bac.time.meta$SampleID
sample_data(bac.time.clean.ss) <- sample_data(bac.time.meta)

## let's get phylum sir
df.phylum <- bac.time.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
#df.phylum$Phylum2[which(df.phylum$Class=="Betaproteobacteria")] <- "Betaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

head(df.phylum)
df.phylum$SampleID_2 <- factor(df.phylum$SampleID_2, levels = order.sample)

df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))
unique(df.phylum$Phylum2)

levels(df.phylum$Phylum2)
levels(df.phylum$Phylum2) = c(levels(df.phylum$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(SampleID_2) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

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
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=SampleID_2, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Deltaproteobacteria" = "darkolivegreen4", "Actinobacteriota"="indianred2","Proteobacteria" = "#003333",
                               "Bacteroidota"="steelblue1", "Firmicutes" ="tan1", 
                               "unidentified" = "black", "Low abundance" = "light grey")) +
  
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

##Fungi
###class
fun.time.meta <- read.table(file = './Time series data/fungi/dynamic DB_time series/sample_metadata.tsv', sep = '\t', header = TRUE)
rownames(fun.time.meta) <- fun.time.meta$SampleID
sample_data(fun.time.clean.ss) <- sample_data(fun.time.meta)


df.class.fun <- fun.time.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.time.rel <- df.class.fun %>%  
  group_by(SampleID_2) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.class.fun.time.rel[df.class.fun.time.rel$RelAbundance < 0.01,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.time.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
. <- append(vec.Low.f, vec.order.f)


df.class.fun.time.rel$Class <- factor(df.class.fun.time.rel$Class, levels = vec.reorder.f) 
df.class.fun.time.rel$SampleID_2 <- factor(df.class.fun.time.rel$SampleID_2, levels = order.sample.time) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.fun.time.rel.p1 <- ggplot(df.class.fun.time.rel, aes(x=SampleID_2, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Dothideomycetes" = "#5195D1","Ustilaginomycetes"="#ffcc33","Sordariomycetes"= "#1E63AF",
                               "Tremellomycetes"="#BE4146","Microbotryomycetes" ="#DC9A9E","Cystobasidiomycetes" = "#A871AE",
                               "Agaricomycetes" = "#CC6C71","Eurotiomycetes"= "#6DA9DC","Agaricostilbomycetes"="#003333",
                               "Malasseziomycetes" = "#99cc99","Saccharomycetes"="#CCCC99", "Spiculogloeomycetes"="#cc66ff",
                               "Leotiomycetes"= "#11335F",  "Lecanoromycetes"="#00cccc",
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
df.class.fun.time.rel.p1

dev.off()



df.class.fun.time.rel.2 <-df.class.fun.time.rel %>% group_by(Sample, Class) %>% summarise(SumRA = sum(RelAbundance))
df.class.bac.time.rel.2 <-df.phylum.rel %>% group_by(Sample, Phylum2) %>% summarise(SumRA = sum(RelAbundance))

write.csv(df.class.fun.time.rel.2,"df.class.fun.time.rel.2.csv")
write.csv(df.class.bac.time.rel.2,"df.class.bac.time.rel.2.csv")


###order
df.order.fun <- fun.time.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()

df.order.fun %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))
levels(df.order.fun$Order) = c(levels(df.order.fun$Order), 'Low abundance')

# we need to group by samples
df.order.fun.time.rel <- df.order.fun %>%  
  group_by(SampleID_2) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.order.fun.time.rel[df.order.fun.time.rel$RelAbundance < 0.01,]$Order <- 'Low abundance'
unique(df.order.fun$Order)

df.order.fun.time.rel$SampleID_2 <- factor(df.order.fun.time.rel$SampleID_2, levels = order.sample.time) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'
df.order.fun.time.rel.2 <-df.order.fun.time.rel %>% group_by(Sample, Order) %>% summarise(SumRA = sum(RelAbundance))

write.csv(df.order.fun.time.rel.2,"Distribution of fungal orders in pooled seeds.csv")



##### Relative abundance table
bac.clean.ss.rel
otu.bac.rel<-otu_table(bac.clean.ss.rel)
otu.bac.rel <- data.frame(otu.bac.rel)
otu.bac.rel$Total <- rowSums(otu.bac.rel)
otu.bac.rel$OTU <- rownames(otu.bac.rel)

bac.list.sub <- bac.list %>% select(OTU, OTU_id, number, Phylum,Class,Order,Family,Genus)

otu.bac.rel.tab<- merge(otu.bac.rel,bac.list.sub, by ="OTU")
write.csv(otu.bac.rel.tab,"Distribution of bacterial ASVs in indivudal seeds.csv")


otu.fun.rel<-otu_table(fun.clean.ss.rel)
otu.fun.rel <- data.frame(otu.fun.rel)
otu.fun.rel$Total <- rowSums(otu.fun.rel)
otu.fun.rel$OTU <- rownames(otu.fun.rel)

fun.list.sub <- fun.list %>% select(OTU, OTU_id, number, Phylum,Class,Order,Family,Genus)

otu.fun.rel.tab<- merge(otu.fun.rel,fun.list.sub, by ="OTU")
write.csv(otu.fun.rel.tab,"Distribution of fungal ASVs in indivudal seeds.csv")




df.phylum <- bac.clean.ss %>%
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
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.01,]$Phylum2 <- 'Low abundance'


df.phylum.rel.2 <-df.phylum.rel %>% group_by(Sample, Phylum2) %>% summarise(SumRA = sum(RelAbundance))
write.csv(df.phylum.rel.2, "Table S3_bacterial phyla and class.csv")



df.phylum <- bac.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

head(df.phylum)
df.phylum$SampleID <- factor(df.phylum$SampleID, levels = order.sample)

df.phylum %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
unique(df.phylum$Class)

levels(df.phylum$Class)
levels(df.phylum$Class) = c(levels(df.phylum$Class), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(SampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.01,]$Class <- 'Low abundance'


df.phylum.rel.2 <-df.phylum.rel %>% group_by(Sample, Class) %>% summarise(SumRA = sum(RelAbundance))
write.csv(df.phylum.rel.2, "Table S3_bacterial class.csv")