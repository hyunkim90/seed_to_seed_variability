alpha.all <-microbiome::dominance(bac.rarefied.18, index = "all")
head(alpha.all)
mean(alpha.all$gini) # 0.9841488
sd(alpha.all$gini) #0.008513602
alpha.all$Cond <- "Individual"


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

##27 samples
random.sample<-c("R34","R29","R17","R13","R43","R23","R12","R4","R44","R8","R32","R14","R22","R21","R56","R60","R67","R58","R28","R45","R27",
                 "R11","R41","R37","R35","R7","R10")

alpha.all.sub<-subset(alpha.all, rownames(alpha.all) %in% random.sample)
DominTab.sub <- alpha.all.sub[c(7,8)]

DominTab.2 <- rbind(Domintab.pooled,DominTab.sub)
write.csv(DominTab.2,"230119_Dominance_Pooled vs Ind_27samples.csv")


