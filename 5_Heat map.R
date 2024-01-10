#### Construction of heat map to show seed to seed variability in time series data
## Extract relative abundance table
otu.bac.time<-otu_table(bac.time.clean.rel)
otu.bac.time<-data.frame(otu.bac.time)

bac.time.list$Phylum[is.na(bac.time.list$Phylum)] <- "Unidentified"
bac.time.list$Class[is.na(bac.time.list$Class)] <- "Unidentified"
bac.time.list$Phylum2 <- bac.time.list$Phylum
bac.time.list$Phylum2[which(bac.time.list$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
bac.time.list$Phylum2[which(bac.time.list$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
bac.time.list$Phylum2[which(bac.time.list$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
bac.time.list$Phylum2[which(bac.time.list$Phylum == "Proteobacteria" & bac.time.list$Class == "Unidentified")] <- "Unidentified Proteobacteria"

unique(bac.time.list$Phylum2)

bac.time.list.phyla2 <- bac.time.list[-c(2:10)]
otu.bac.time$OTU <- rownames(otu.bac.time)
otu.bac.time.phyla <- merge(otu.bac.time, bac.time.list.phyla2)

write.xlsx(otu.bac.time.phyla, "Heat map rawdata_20211216_bacteria.xlsx")


otu.fun.time<-otu_table(fun.time.clean.rel)
otu.fun.time<-data.frame(otu.fun.time)

fun.time.list$Phylum[is.na(fun.time.list$Phylum)] <- "Unidentified"
fun.time.list$Class[is.na(fun.time.list$Class)] <- "Unidentified"
fun.time.list$Class2 <- fun.time.list$Class
fun.time.list$Class2[which(fun.time.list$Class == "Unidentified"& fun.time.list$Phylum == "Ascomycota")] <- "Unidentified Ascomycota"
fun.time.list$Class2[which(fun.time.list$Class == "Unidentified"& fun.time.list$Phylum == "Basidiomycota")] <- "Unidentified Basidiomycota"

unique(fun.time.list$Class2)

fun.time.list.class2 <- fun.time.list[-c(3:10)]
otu.fun.time$OTU <- rownames(otu.fun.time)
otu.fun.time.class <- merge(otu.fun.time, fun.time.list.class2)

write.xlsx(otu.fun.time.class, "Heat map rawdata_20211216_fungi.xlsx")




fun.time.list$OTU_id[which(fun.time.list$OTU %in%c("6d4c2e120911d01408a45f767588460e", "194a950bddd6144c5d7bdcf39ae09110"))]



### Supplementary Table
bac.heat.tab <- read.xlsx("Heat map rawdata_20211216_bacteria.xlsx",1)
fun.heat.tab <- read.xlsx("Heat map rawdata_20211216_fungi.xlsx",1)

bac.heat.tab.2<-merge(bac.heat.tab, bac.time.list, by = "OTU")
write.xlsx(bac.heat.tab.2, "Supplementary Table_220520_bacteria.xlsx")

fun.heat.tab.2<-merge(fun.heat.tab, fun.time.list, by = "OTU")
write.xlsx(fun.heat.tab.2, "Supplementary Table_220520_fungi.xlsx")


##### Individual seed samples
otu.bac<-otu_table(bac.clean.rel)
otu.bac<-data.frame(otu.bac)

bac.list$Phylum[is.na(bac.list$Phylum)] <- "Unidentified"
bac.list$Class[is.na(bac.list$Class)] <- "Unidentified"
bac.list$Phylum2 <- bac.list$Phylum
bac.list$Phylum2[which(bac.list$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
bac.list$Phylum2[which(bac.list$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
bac.list$Phylum2[which(bac.list$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
bac.list$Phylum2[which(bac.list$Phylum == "Proteobacteria" & bac.list$Class == "Unidentified")] <- "Unidentified Proteobacteria"

unique(bac.list$Phylum2)

bac.list.phyla2 <- bac.list[-c(2:10)]
otu.bac$OTU <- rownames(otu.bac)
otu.bac.phyla <- merge(otu.bac, bac.list.phyla2)

write.csv(otu.bac.phyla, "Heat map rawdata_20230119_bacteria.csv")


otu.fun<-otu_table(fun.clean.rel)
otu.fun<-data.frame(otu.fun)

fun.list$Phylum[is.na(fun.list$Phylum)] <- "Unidentified"
fun.list$Class[is.na(fun.list$Class)] <- "Unidentified"
fun.list$Class2 <- fun.list$Class
fun.list$Class2[which(fun.list$Class == "Unidentified"& fun.list$Phylum == "Ascomycota")] <- "Unidentified Ascomycota"
fun.list$Class2[which(fun.list$Class == "Unidentified"& fun.list$Phylum == "Basidiomycota")] <- "Unidentified Basidiomycota"

unique(fun.list$Class2)

fun.list.class2 <- fun.list[-c(3:10)]
otu.fun$OTU <- rownames(otu.fun)
otu.fun.class <- merge(otu.fun, fun.list.class2)

write.csv(otu.fun.class, "Heat map rawdata_20230119_fungi.csv")


### Supplementary Table
bac.ind.heat.tab <- read.csv("Heat map rawdata_20230119_bacteria.csv")
fun.ind.heat.tab <- read.csv("Heat map rawdata_20230119_fungi.csv")

bac.ind.heat.2<-merge(bac.ind.heat.tab, bac.list, by = "OTU")
write.csv(bac.ind.heat.2, "Supplementary Table_230119_bacteria.csv")

fun.ind.heat.2<-merge(fun.ind.heat.tab, fun.list, by = "OTU")
write.csv(fun.ind.heat.2, "Supplementary Table_230119_fungi.csv")

