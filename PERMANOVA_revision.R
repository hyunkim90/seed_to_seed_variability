##PERMANOVA
f.otu <- otu_table(fun.time.clean.log)
f.meta <- data.frame(sample_data(fun.time.clean.log))
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ (Year+Replication+Plot), data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication*Plot), data = f.meta, permutations=9999, method = "jaccard")
f.permanova

b.otu <- otu_table(bac.time.clean.log)
b.meta <- data.frame(sample_data(bac.time.clean.log))

## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ (Year+Replication+Plot), data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication*Plot), data = b.meta, permutations=9999, method = "jaccard")
b.permanova


##PERMANOVA
f.otu <- otu_table(fun.time.clean.ss)
f.meta <- data.frame(sample_data(fun.time.clean.ss))
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ (Year+Replication+Plot), data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication*Plot), data = f.meta, permutations=9999, method = "jaccard")
f.permanova

b.otu <- otu_table(bac.time.clean.ss)
b.meta <- data.frame(sample_data(bac.time.clean.ss))
## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ (Year+Replication+Plot), data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Year*Replication), data = b.meta, permutations=9999, method = "jaccard")
b.permanova


##PERMANOVA
f.otu <- otu_table(fun.time.clean.ss)
f.meta <- data.frame(sample_data(fun.time.clean.ss))
f.meta$Replication2 <- paste0(f.meta$Replication,"_",f.meta$Plot)
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication2), data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication2), data = f.meta, permutations=9999, method = "jaccard")
f.permanova

anosim(x = t(f.otu), grouping = f.meta$Plot, permutations = 9999, distance = "bray") 

b.otu <- otu_table(bac.time.clean.ss)
b.meta <- data.frame(sample_data(bac.time.clean.ss))
b.meta$Replication2 <- paste0(b.meta$Replication,"_",b.meta$Plot)
## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication2), data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication2), data = b.meta, permutations=9999, method = "jaccard")
b.permanova


f.otu <- otu_table(fun.time.clean.log)
f.meta <- data.frame(sample_data(fun.time.clean.log))
f.meta$Replication2 <- paste0(f.meta$Replication,"_",f.meta$Plot)
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication2), data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ (Replication2), data = f.meta, permutations=9999, method = "jaccard")
f.permanova

b.otu <- otu_table(bac.time.clean.log)
b.meta <- data.frame(sample_data(bac.time.clean.log))
b.meta$Replication2 <- paste0(b.meta$Replication,"_",b.meta$Plot)
## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication2), data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication2), data = b.meta, permutations=9999, method = "jaccard")
b.permanova



bac.time.clean.log.18


b.otu <- otu_table(bac.time.clean.log.18)
b.meta <- data.frame(sample_data(bac.time.clean.log.18))

## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ Plot, data = b.meta, permutations=9999, method = "bray")
b.permanova
## Jaccard
b.permanova <- adonis2(formula = t(b.otu) ~ (Replication2), data = b.meta, permutations=9999, method = "jaccard")
b.permanova

fun.time.clean.log.18

f.otu <- otu_table(fun.time.clean.log.18)
f.meta <- data.frame(sample_data(fun.time.clean.log.18))
## Bray-Curtis
f.permanova <- adonis2(formula = t(f.otu) ~ Plot, data = f.meta, permutations=9999, method = "bray")
f.permanova
## Jaccard
f.permanova <- adonis2(formula = t(f.otu) ~ Plot, data = f.meta, permutations=9999, method = "jaccard")
f.permanova



bac.time.clean.log.18.2<-subset_samples(bac.time.clean.log, Year != "year_2017")
b.otu <- otu_table(bac.time.clean.log.18.2)
b.meta <- data.frame(sample_data(bac.time.clean.log.18.2))

## Bray-Curtis
b.permanova <- adonis2(formula = t(b.otu) ~ Developmental_stage, data = b.meta, permutations=9999, method = "bray")
b.permanova
