### 

##Abundance table
### Absolute
table1<-data.frame(otu_table(bac.clean.ss.f))
write.csv(table1, "asv_abund_all_bacteria.csv")
taxtable1<-data.frame(tax_table(bac.clean.ss.f))
write.csv(taxtable1, "asv_tax_all_bacteria.csv")
metatable1<-data.frame(sample_data(bac.clean.ss.f))
write.csv(metatable1, "asv_meta_all_bacteria.csv")


table1<-data.frame(otu_table(fun.clean.ss.f))
write.csv(table1, "asv_abund_all_fungi.csv")
taxtable1<-data.frame(tax_table(fun.clean.ss.f))
write.csv(taxtable1, "asv_tax_all_fungi.csv")

## Normalized
table2<-data.frame(otu_table(bac.clean.nolog))
write.csv(table2, "asv_abund_norm_bacteria.csv")
taxtable2<-data.frame(tax_table(bac.clean.nolog))
write.csv(taxtable2, "asv_tax_norm_bacteria.csv")
metatable2<-data.frame(sample_data(bac.clean.nolog))
write.csv(metatable2, "asv_meta_norm_bacteria.csv")


table2<-data.frame(otu_table(fun.clean.nolog))
write.csv(table2, "asv_abund_norm_fungi.csv")
taxtable2<-data.frame(tax_table(fun.clean.nolog))
write.csv(taxtable2, "asv_tax_norm_fungi.csv")

fun.clean.rel <- transform(fun.clean.ss.f, transform = "compositional")
### Relative abundance
table3<-data.frame(otu_table(bac.clean.rel))
write.csv(table3, "asv_abund_rel_bacteria.csv")

table3<-data.frame(otu_table(fun.clean.rel))
write.csv(table3, "asv_abund_rel_fungi.csv")

##Classification of ASVs according to prevalence
## core, ASVs with 80% of prevalence; rare, ASVs with 20% of prevalence; others, remaining ASVs
bac.prev.tab$Type <- 0
fun.prev.tab$Type <- 0

bac.prev.tab$Type[which(bac.prev.tab$Prevalence <= 0.2)] <- "Rare"
fun.prev.tab$Type[which(fun.prev.tab$Prevalence <= 0.2)] <- "Rare"

bac.prev.tab$Type[which(bac.prev.tab$Prevalence >= 0.8)] <- "Core"
fun.prev.tab$Type[which(fun.prev.tab$Prevalence >= 0.8)] <- "Core"

bac.prev.tab$Type[which(bac.prev.tab$Type == 0)] <- "Others"
fun.prev.tab$Type[which(fun.prev.tab$Type == 0)] <- "Others"

write.csv(bac.prev.tab, "Type of bacterial ASVs.csv")
write.csv(fun.prev.tab, "Type of fungi ASVs.csv")
# version 2020.8.23
# version 2020.9.21, add classification information
# version 2021.1.7, add (step 15) icamp.cate to summarize for different categories of taxa, e.g. core versus rare taxa.
# version 2021.4.17, add step 9.5 and 9.6.

rm(list=ls())
t0=Sys.time() # to calculate time cost

# 1 # key parameter setting
prefix="SeedVariable_211229"  # prefix of the output file names. usually use a project ID.
rand.time=1000  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=5 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

# 2 # load R packages and data
library(iCAMP)
library(ape)
load.wd <- "C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation/Assembly mechanism/Input"
setwd(load.wd)

##Bacteria
comm=t(read.table("asv_abund_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=phy_tree(bac.clean.ss.f)
clas=read.table("asv_tax_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table("asv_meta_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

## Fungi
comm=t(read.table("asv_abund_all_fungi.txt", header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=phy_tree(fun.clean.ss.f)
clas=read.table("asv_tax_all_fungi.txt", header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table("asv_meta_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

### Normalized abundance
##Bacteria
comm=t(read.table("asv_abund_norm_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=phy_tree(bac.clean.nolog)
clas=read.table("asv_tax_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table("asv_meta_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

## Fungi
comm=t(read.table("asv_abund_norm_fungi.txt", header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=phy_tree(fun.clean.nolog)
clas=read.table("asv_tax_all_fungi.txt", header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table("asv_meta_all_bacteria.txt", header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

# 4 # match sample IDs in OTU table and treatment information table
#sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
#env=sampid.check$env # skip this if you do not have env.file

# 5 # match OTU IDs in OTU table and tree file
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed.
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 
##Abs
save.wd <- "C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation/Assembly mechanism/Bacteria_all"
save.wd <- "C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation/Assembly mechanism/Fungi_all"

##Norm
save.wd <- "C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation/Assembly mechanism/Bacteria_all/Normalized"
save.wd <- "C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation/Assembly mechanism/Fungi_all/Normalized"

setwd(save.wd)

if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  # output files:
  # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
  # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
  # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
  # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

####################
# you may skip step 7-8, if the "alternative way" based on stochasticity is applicable, as mentioned in the method part of iCAMP paper (Ning et al 2020 Nature Communications).
# 7 # assess niche preference difference between species
# env is required for this step.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large niche difference matrix. 
setwd(save.wd)
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)

# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 5 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)
# 8.2 # test within-bin phylogenetic signal.
sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)
# since this example small data is randomly generated, the correlation should be very weak.
# usually, you are looking for a binning setting lead to higher RAsig.abj (relative abundance of bins with significant phylogenetic signal) and relative high meanR (mean correlation coefficient across bins).
# see help document of the function "ps.bin" for the meaning of output.

####################
# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
# commonly use # set sig.index as Confidence instead of SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 12 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
# there are quite a few parameters in this function, please check the help document of "icamp.big".
# output files:
# Test.iCAMP.detail.rda: the object "icres" saved in R data format. it is a list object. The first element bNRIiRCa is the result of relative importance of each assembly process in each pairwise comparison (each turnover). The second element "detail" including binning information (named taxabin), phylogenetic and taxonomic metrics results in each bin (named like bNRIi, RCa, etc.), relative abundance of each bin (bin.weight), relative importance of each process in each turnover between communities (processes), input settings (setting), and input community data matrix (comm). See help document of the function icamp.big for more details.

##Abs
bac.all.abs.icres<-icres
fun.all.abs.icres<-icres

##Norm
bac.all.norm.icres<-icres
fun.all.norm.icres<-icres

##Norm_bin size 12 (confidence)
bac.all.norm.icres.2<-icres
fun.all.norm.icres.2<-icres

############################
# 9.2 to 9.4 are some optional special settings you may explore.
# 9.2 # explore different ways for null model significance test.
# 9.2.1 # set detail.null=TRUE, output all null values, to facilitate normality test and switch between different options
detail.null=TRUE
bin.size.limit = 5 
sig.index="SES.RC" # this is traditional way, with assumption that null values of phylogenetic metrics follow normal distribution. 
prefixb="TestB"

icres2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                        pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                        prefix = prefixb, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                        phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                        phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                        nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                        qp.save = FALSE, detail.null = detail.null, ignore.zero = TRUE, output.wd = save.wd, 
                        correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                        ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
# 9.2.2 # normality test
nntest=iCAMP::null.norm(icamp.output=icres2, p.norm.cut=0.05, detail.out=FALSE)
# output shows non-normal distribution ratio in each bin, i.e. the proportion of turnovers which have null values significantly deviated from normal distribution.
# if some ratio values are very high, may need to change to use "Confidence" as sig.index.

# 9.2.3 # change sig.index to "Confidence".
icres3=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "Confidence", detail.save = TRUE, detail.null = FALSE, conf.cut = 0.975)
head(icres3$CbMPDiCBraya)

# 9.2.4 # change sig.index to "RC" for both phylogenetic and taxonomic metrics.
icres4=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "RC", detail.save = TRUE, detail.null = FALSE, rc.cut = 0.95)
head(icres4$RCbMPDiRCbraya)

# 9.2.5 # the function can also change the significance threshold.
icres5=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "SES.RC", detail.save = TRUE, detail.null = FALSE, ses.cut = 1.64, rc.cut = 0.9)
head(icres5$bNRIiRCbraya)

# 9.3 # you may specify the relative abundance of each species in the regional pool, if it is not the same as the average relative abundance from the "comm" you input.
meta.ab=rep(1,ncol(comm)) # here i assume all the species actuall have the same relative abundance in the regional pool.
prefix2=paste0(prefix,".MetaCrct")
sig.index="Confidence" # see other options in help document of icamp.big.
icres.meta=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix2, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                            qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab=meta.ab)

# 9.4 # consider to omit small bins
# 9.4.1 # if you would like to omit small bins rather than merging them to nearest relatives, set omit.option as "test" to check what will be omitted.
omit.option = "test"
icres.omit=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                            qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)
# "test" will return a detailed table of omitted species.

# 9.4.2 # then set it as "omit" to omit the small bins.
omit.option = "omit"
icres.omit2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                             pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                             prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                             phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                             phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                             nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                             qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                             correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                             ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)
# In this simple example, since all bins are small, "omit" should return an error. In real data, this will go ahead to do iCAMP analysis with the strict bins which are big enough (>bin.size.limit).


# 9.6 # community data transformation and taxonomic dissimilarity index change
taxo.metric='euclidean'
transform.method='hellinger'
prefixtran=paste0(prefix,"Hel")
bin.size.limit = 5 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
icres7=iCAMP::icamp.big(comm=comm,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixtran,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric=taxo.metric, transform.method=transform.method,
                        logbase=2, dirichlet=FALSE)


###############################
# 10 # iCAMP bin level statistics
icres <- bac.all.abs.icres
icres <- fun.all.abs.icres

icres <- bac.all.norm.icres 
icres <- fun.all.norm.icres 

icres <- bac.all.norm.icres.2 
icres <- fun.all.norm.icres.2 

icbin=icamp.bins(icamp.detail = icres$detail,treat = treat,
                        clas=clas,silent=FALSE, boot = TRUE,
                        rand.time = rand.time,between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)

# output files:
# Test.iCAMP.Summary.rda: the object "icbin" saved in R data format. see help document of the function icamp.bins for description of each element in the object.
# Test.ProcessImportance_EachGroup.csv: Relative importance of each process in governing the turnovers in a group of samples.
# Test.ProcessImportance_EachBin_EachGroup.csv: Relative importance of each process in governing the turnovers of each bin among a group of samples.
# Test.ProcessImportance_EachTurnover.csv: Relative importance of each process in governing the turnovers between each pair of communities (samples).
# Test.BinContributeToProcess_EachGroup.csv: Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.
# Test.Taxon_Bin.csv: a matrix showing the bin ID and classification information for each taxon.
# Test.Bin_TopTaxon.csv: a matrix showing the bin relative abundance; the top taxon ID, percentage in bin, and classification; the most abundant name at each phylogeny level in the bin.

# 11 # Bootstrapping test
# please specify which column in the treatment information table.
i=1
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)

# output files:
# Test.iCAMP.Boot.Management.rda: the object "icboot" saved in R data format. see help document of the function icamp.boot for description of each element in the object.
# Test.BootSummary.Management.csv: a table to summarize bootstrapping results. see help document of the function icamp.boot for description of the output element "summary".
# Test.Compare.Management.csv: a table to summarize comparison index, effect size, and significance between each two groups. see help document of the function icamp.boot for description of the output element "compare".
# 12 # Other approach: QPEN (quantifying community assembly processes based on entire-community null model analysis)
# 12.1 # QPEN calculation
qpout=iCAMP::qpen(comm=comm,pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                  pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                  rand.time=rand.time, nworker=nworker,project=prefix,
                  wd=save.wd, save.bNTIRC=TRUE)

##Abs
classic.assem.imp.bac<-qpout
classic.assem.imp.fun<-qpout

##Normal
classic.assem.imp.bac.norm<-qpout
classic.assem.imp.fun.norm<-qpout

##Normal bin 12
classic.assem.imp.bac.norm.2<-qpout
classic.assem.imp.fun.norm.2<-qpout

write.csv(classic.assem.imp.bac$result, "Classical relative importance_rawdata_abs_bacteria.csv")
write.csv(classic.assem.imp.bac$ratio, "Classical relative importance_ratio_abs_bacteria.csv")

write.csv(classic.assem.imp.fun$result, "Classical relative importance_rawdata_abs_fungi.csv")
write.csv(classic.assem.imp.fun$ratio, "Classical relative importance_ratio_abs_fungi.csv")

write.csv(classic.assem.imp.bac.norm$result, "Classical relative importance_rawdata_norm_bacteria.csv")
write.csv(classic.assem.imp.bac.norm$ratio, "Classical relative importance_ratio_norm_bacteria.csv")

write.csv(classic.assem.imp.fun.norm$result, "Classical relative importance_rawdata_norm_fungi.csv")
write.csv(classic.assem.imp.fun.norm$ratio, "Classical relative importance_ratio_norm_fungi.csv")

write.csv(classic.assem.imp.bac.norm.2$result, "Classical relative importance_rawdata_norm_bacteria_bin12.csv")
write.csv(classic.assem.imp.bac.norm.2$ratio, "Classical relative importance_ratio_norm_bacteria_bin12.csv")

write.csv(classic.assem.imp.fun.norm.2$result, "Classical relative importance_rawdata_norm_fungi_bin12.csv")
write.csv(classic.assem.imp.fun.norm.2$ratio, "Classical relative importance_ratio_norm_fungi_bin12.csv")

# # 12.2 # significance test
# qptest=qpen.test(qpen.result = qpout,treat = treat,rand.time = rand.time,
#                  between.group = TRUE,out.detail=TRUE,silent=FALSE)
# write.csv(qptest$obs.summary,file = paste0(prefix,".QPEN.Index.Obs.Summary.csv"),row.names = FALSE)
# write.csv(qptest$boot.summary,file = paste0(prefix,".QPEN.Bootstrapping.Summary.csv"),row.names = FALSE)
# write.csv(qptest$compare,file = paste0(prefix,".QPEN.Comparison.Summary.csv"),row.names = FALSE)
# save(qptest,file = paste0(prefix,".QPEN.bootstrap.rda"))
# 
# # 13 # Other approach: Neutral taxa percentage
# snmout=iCAMP::snm.comm(comm = comm, treat = treat, 
#                        rand=rand.time, alpha=0.05)
# write.csv(snmout$stats,file = paste0(prefix,".NeutralModel.Stats.csv"))
# write.csv(snmout$ratio.summary,file = paste0(prefix,".NeutralModel.TypeRatio.csv"))
# 
# # 14 # Other approach: tNST and pNST (taxonomic and phylogenetic normalized stochasticity ratio)
# # need to install package NST if not yet
# if(!("NST" %in% installed.packages()[,"Package"])){install.packages("NST")}
# library(NST)
# i=1
# treat.use=treat[,i,drop=FALSE]
# 
# # 14.1a # tNST
# tnstout=NST::tNST(comm=comm, group=treat.use, dist.method="bray", 
#                   abundance.weighted=TRUE, rand=rand.time,  
#                   nworker=nworker, null.model="PF", output.rand = TRUE,
#                   SES = TRUE, RC = TRUE)
# write.csv(tnstout$index.grp,file = paste0(prefix,".tNST.summary.",colnames(treat)[i],".csv"))
# write.csv(tnstout$index.pair.grp,file = paste0(prefix,".tNST.pairwise.",colnames(treat)[i],".csv"))
# 
# # 14.1b # bootstrapping test for tNST
# tnst.bt=NST::nst.boot(nst.result=tnstout, group=treat.use,
#                       rand=rand.time, nworker=nworker)
# write.csv(tnst.bt$NST.summary,file = paste0(prefix,".tNST.bootstr.",colnames(treat)[i],".csv"))
# write.csv(tnst.bt$NST.compare,file = paste0(prefix,".tNST.compare.",colnames(treat)[i],".csv"))
# 
# # 14.2a # pNST
# pnstout=NST::pNST(comm=comm, pd.desc=pd.big$pd.file, pd.wd=pd.big$pd.wd, 
#                   pd.spname=pd.big$tip.label, group=treat.use, abundance.weighted=TRUE,
#                   rand=rand.time, phylo.shuffle=TRUE, nworker=nworker,
#                   output.rand = TRUE, SES=FALSE, RC=FALSE)
# write.csv(pnstout$index.grp,file = paste0(prefix,".pNST.summary.",colnames(treat)[i],".csv"))
# write.csv(pnstout$index.pair.grp,file = paste0(prefix,".pNST.pairwise.",colnames(treat)[i],".csv"))
# 
# pnst.bt=NST::nst.boot(nst.result=pnstout, group=treat.use,
#                       rand=rand.time, nworker=nworker)
# write.csv(pnst.bt$NST.summary,file = paste0(prefix,".pNST.bootstr.",colnames(treat)[i],".csv"))
# write.csv(pnst.bt$NST.compare,file = paste0(prefix,".pNST.compare.",colnames(treat)[i],".csv"))

# 15.1 # define the types of different taxa in category.txt
setwd(load.wd)

cate.file="Type of bacterial ASVs.txt"
cate.file="Type of fungi ASVs.txt"

cate=read.table(cate.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
cate=cate[which(rownames(cate) %in% colnames(comm)),,drop=FALSE] # remove unmatched taxa.
setwd(save.wd)

# 15.2
iccate=icamp.cate(icamp.bins.result = icbin,comm = comm,cate = cate,
                  treat = treat, silent = FALSE,between.group = TRUE)
write.csv(iccate$Ptuvx,file = paste0(prefix,".iCAMP.Process_EachTurnover_EachCategory.csv"))
write.csv(iccate$Ptx,file = paste0(prefix,".iCAMP.Process_EachGroup_EachCategory.csv"))

#### Summarize the ICAMP result and compare it to conventional QPEN 
bac.icamp<-bac.all.norm.icres.2$CbMPDiCBraya
fun.icamp<-fun.all.norm.icres.2$CbMPDiCBraya

bac.all.norm.icres.2$detail$taxabin

### Overall result
bac.icamp.num<-bac.icamp[-c(1,2)]
bac.icamp.summary<-colMeans(bac.icamp.num)

fun.icamp.num<-fun.icamp[-c(1,2)]
fun.icamp.summary<-colMeans(fun.icamp.num)

### Classical 
setwd("C:/Users/Hyun Kim/Desktop/SNU/SNU_Taeon/Seed variation")

bac.classic <- read.csv('./Assembly mechanism/Bacteria_all/Normalized/Classical relative importance_ratio_norm_bacteria_bin12.csv')
fun.classic <- read.csv('./Assembly mechanism/Fungi_all/Normalized/Classical relative importance_ratio_norm_fungi_bin12.csv')

### bind ICAMP and QPEN data frames
#Bacteria
df.bac.icamp.summary <- data.frame(bac.icamp.summary)
names(df.bac.icamp.summary)[1] <- "ICAMP"

df.bac.classic <- t(data.frame(bac.classic))
colnames(df.bac.classic) <- "QPEN"

df.bac.rel.imp.tab<-cbind(df.bac.icamp.summary,df.bac.classic)

#Fungi
df.fun.icamp.summary <- data.frame(fun.icamp.summary)
names(df.fun.icamp.summary)[1] <- "ICAMP"

df.fun.classic <- t(data.frame(fun.classic))
colnames(df.fun.classic) <- "QPEN"

df.fun.rel.imp.tab<-cbind(df.fun.icamp.summary,df.fun.classic)

### Melt data frames
df.bac.rel.imp.tab$Mechanism <- rownames(df.bac.rel.imp.tab)
df.fun.rel.imp.tab$Mechanism <- rownames(df.fun.rel.imp.tab)

df.bac.rel.imp.tab.melt <- melt(df.bac.rel.imp.tab)
names(df.bac.rel.imp.tab.melt)[2] <- "Tool"
names(df.bac.rel.imp.tab.melt)[3] <- "RelImportance"
df.bac.rel.imp.tab.melt$Kingdom <- "Bacteria"

df.fun.rel.imp.tab.melt <- melt(df.fun.rel.imp.tab)
names(df.fun.rel.imp.tab.melt)[2] <- "Tool"
names(df.fun.rel.imp.tab.melt)[3] <- "RelImportance"
df.fun.rel.imp.tab.melt$Kingdom <- "Fungi"

## Plotting
## Bacteria
df.bac.rel.imp.tab.melt$Mechanism <- factor(df.bac.rel.imp.tab.melt$Mechanism, levels = rev(c("Homogeneous.Selection","Heterogeneous.Selection",
                                                                                          "Dispersal.Limitation", "Homogenizing.Dispersal",
                                                                                          "Drift.and.Others")))
ggplot(df.bac.rel.imp.tab.melt, aes(x=Tool, y = RelImportance, fill = Mechanism)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Homogeneous.Selection"= "#3c93c2","Heterogeneous.Selection" = "#9ec9e2",
                               "Dispersal.Limitation" = "#fcde9c", "Homogenizing.Dispersal" = "#feb24c",
                               "Drift.and.Others"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


df.fun.rel.imp.tab.melt$Mechanism <- factor(df.fun.rel.imp.tab.melt$Mechanism, levels = rev(c("Homogeneous.Selection","Heterogeneous.Selection",
                                                                                              "Dispersal.Limitation", "Homogenizing.Dispersal",
                                                                                              "Drift.and.Others")))
ggplot(df.fun.rel.imp.tab.melt, aes(x=Tool, y = RelImportance, fill = Mechanism)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Homogeneous.Selection"= "#3c93c2","Heterogeneous.Selection" = "#9ec9e2",
                               "Dispersal.Limitation" = "#fcde9c", "Homogenizing.Dispersal" = "#feb24c",
                               "Drift.and.Others"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


dev.off()

### Bacteria and fungi_iCAMP
df.bac.rel.imp.tab.melt
df.fun.rel.imp.tab.melt

df.rel.imp.tab <- rbind(df.bac.rel.imp.tab.melt,df.fun.rel.imp.tab.melt)
df.rel.imp.tab$Mechanism <- factor(df.rel.imp.tab$Mechanism, levels = rev(c("Homogeneous.Selection","Heterogeneous.Selection",
                                                                                              "Dispersal.Limitation", "Homogenizing.Dispersal",
                                                                                              "Drift.and.Others")))
ggplot(subset(df.rel.imp.tab, Tool == "ICAMP"), aes(x=Kingdom, y = RelImportance, fill = Mechanism)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Homogeneous.Selection"= "#3c93c2","Heterogeneous.Selection" = "#9ec9e2",
                               "Dispersal.Limitation" = "#fcde9c", "Homogenizing.Dispersal" = "#feb24c",
                               "Drift.and.Others"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

ggplot(df.rel.imp.tab, aes(x=Tool, y = RelImportance, fill = Mechanism)) + 
  facet_wrap(~Kingdom)+
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Homogeneous.Selection"= "#3c93c2","Heterogeneous.Selection" = "#9ec9e2",
                               "Dispersal.Limitation" = "#fcde9c", "Homogenizing.Dispersal" = "#feb24c",
                               "Drift.and.Others"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


dev.off()

### Relative importance in each phylogenetic groups (bins)
##Bacteria
## Data 1: Tip information (taxonomic information)
bac.list.2 <- bac.list
bac.list.2$Phylum[is.na(bac.list.2$Phylum)] <- "Unidentified"
bac.list.2$Class[is.na(bac.list.2$Class)] <- "Unidentified"
bac.list.2$Phylum2 <- bac.list.2$Phylum
bac.list.2$Phylum2[which(bac.list.2$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
bac.list.2$Phylum2[which(bac.list.2$Phylum == "Proteobacteria" & bac.list.2$Class == "Unidentified")] <- "Unidentified Proteobacteria"

bac.list.2[-c(2:10)]
data1 <- bac.list.2[-c(2:10)]
names(data1)[1]<-"ID"

##top taxon in each bin
bac.phylo.bin.top.taxon<-read.csv("./Assembly mechanism/Bacteria_all/Normalized/SeedVariable.Bin_TopTaxon.csv")
data1$Shape <- 0
data1$Shape <- ifelse(data1$ID %in% bac.phylo.bin.top.taxon$TopTaxonID, "star","others")
unique(data1$Shape)

data1$Shape <- as.factor(data1$Shape)
## Data 2: Phylogenetic bins
bac.phylo.bin<-read.csv("./Assembly mechanism/Bacteria_all/Normalized/SeedVariable.Taxon_Bin.csv")
bac.phylo.bin <- bac.phylo.bin[c(1,2)]
rownames(bac.phylo.bin) <- bac.phylo.bin$ID
bac.phylo.bin <- bac.phylo.bin[-c(1)]

data2 <- bac.phylo.bin
data2$Bin <- factor(data2$Bin, levels = c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7",
                                          "Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14","Bin15",
                                          "Bin16","Bin17","Bin18","Bin19","Bin20","Bin21","Bin22","Bin23"))

tree <- phyloseq::phy_tree(bac.clean.ss.f)
# tree$Bin <- 0
# for (i in tree$tip.label){
#   tree$Bin[which(tree$tip.label == i)] <- as.character(data2$Bin[rownames(data2) == i])
# }

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=0,branch.length='none')
p <- p %<+% data1 + geom_star(
  mapping=aes(fill=Phylum2, starshape=factor(Shape)),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Deltaproteobacteria" = "darkolivegreen4", "Unidentified Proteobacteria" = "#003333","Actinobacteriota"="indianred2",
                               "Firmicutes" ="tan1","Bacteroidota"="steelblue1","Cyanobacteria"="#557882","Planctomycetota"="#dc00f4",
                               "Acidobacteriota"="#F5DC50","Patescibacteria"="#64508C","Bdellovibrionota"="#9BD2D2", "Deinococcota"="#ffaade",
                               "Unidentified" = "black"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)

p<-p + new_scale_fill()

p <- gheatmap(p, data2, offset=.8, width=.1,
                                      colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="Bin")

p

### Fungi


##Bacteria
## Data 1: Tip information (taxonomic information)
fun.list.2 <- fun.list
fun.list.2$Phylum[is.na(fun.list.2$Phylum)] <- "Unidentified"
fun.list.2$Class[is.na(fun.list.2$Class)] <- "Unidentified"
fun.list.2$Class2 <- fun.list.2$Class
fun.list.2$Class2[which(fun.list.2$Class == "Unidentified"& fun.list.2$Phylum == "Ascomycota")] <- "Unidentified Ascomycota"
fun.list.2$Class2[which(fun.list.2$Class == "Unidentified"& fun.list.2$Phylum == "Basidiomycota")] <- "Unidentified Basidiomycota"

fun.list.2[-c(2:10)]
data1 <- fun.list.2[-c(2:10)]
names(data1)[1]<-"ID"

##top taxon in each bin
fun.phylo.bin.top.taxon<-read.csv("./Assembly mechanism/Fungi_all/Normalized/SeedVariable.Bin_TopTaxon.csv")
data1$Shape <- 0
data1$Shape <- ifelse(data1$ID %in% fun.phylo.bin.top.taxon$TopTaxonID, "star","others")
unique(data1$Shape)

data1$Shape <- as.factor(data1$Shape)

## Data 2: Phylogenetic bins
fun.phylo.bin<-read.csv("./Assembly mechanism/Fungi_all/Normalized/SeedVariable.Taxon_Bin.csv")
fun.phylo.bin <- fun.phylo.bin[c(1,2)]
rownames(fun.phylo.bin) <- fun.phylo.bin$ID
fun.phylo.bin <- fun.phylo.bin[-c(1)]

data2 <- fun.phylo.bin
data2$Bin <- factor(data2$Bin, levels = c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7",
                                          "Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14",
                                          "Bin15","Bin16","Bin17"))

# The circular layout tree.
tree <- phyloseq::phy_tree(fun.clean.ss.f)
p <- ggtree(tree, layout="fan", size=0.15, open.angle=0,branch.length='none')
p <- p %<+% data1 + geom_star(
  mapping=aes(fill=Class2, starshape=Shape),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values = c("Dothideomycetes" = "#5195D1","Sordariomycetes"= "#1E63AF", "Eurotiomycetes"= "#6DA9DC",
                               "Leotiomycetes"= "#11335F","Saccharomycetes"="#758C32","Pezizomycetes"="#ECE9B5",
                               "Pezizomycotina_cls_Incertae_sedis"="#60848E","Unidentified Ascomycota" = "#13008A",
                               "Ustilaginomycetes"="#ffcc33","Tremellomycetes"="#BE4146","Cystobasidiomycetes" = "#A871AE",
                               "Microbotryomycetes" ="#DC9A9E","Agaricomycetes" = "#CC6C71","Malasseziomycetes" = "#99cc99",
                               "Agaricostilbomycetes"="#003333",
                               "Unidentified Basidiomycota"="#6A0014","Spizellomycetes" = "#D2C5BC",
                               "Unidentified" ="black"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)

p<-p + new_scale_fill()

p <- gheatmap(p, data2, offset=.8, width=.1,
              colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="Bin")

p


### Relative importance of ecological mechanism in each bin
##Bacteria
bac.phylo.bin<-read.csv("./Assembly mechanism/Bacteria_all/Normalized/SeedVariable.ProcessImportance_EachBin_EachGroup.2.csv")
bac.phylo.bin.melt <- melt(bac.phylo.bin)
head(bac.phylo.bin.melt)
bac.phylo.bin.melt<-subset(bac.phylo.bin.melt, value != "NaN")
bac.phylo.bin.melt.mean <- bac.phylo.bin.melt %>% group_by(Index, variable) %>% summarise(MeanRelImportance = mean(value))
head(bac.phylo.bin.melt.mean)
names(bac.phylo.bin.melt.mean)[2] <- "Bin"

write.csv(fun.phylo.bin.melt.mean,"fun.phylo.bin.melt.mean.csv")

data3 <- bac.phylo.bin.melt.mean
data3$Bin <- factor(data3$Bin, levels = c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7",
                                          "Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14","Bin15",
                                          "Bin16","Bin17","Bin18","Bin19","Bin20","Bin21","Bin22","Bin23"))

data3$Index <- factor(data3$Index, levels = rev(c("HoS","HeS","DL","HD","DR")))

ggplot(data3, aes(x=Bin, y = MeanRelImportance, fill = Index)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("HoS"= "#3c93c2","HeS" = "#9ec9e2",
                               "DL" = "#fcde9c", "HD" = "#feb24c",
                               "DR"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 1) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

##Fungi
fun.phylo.bin<-read.csv("./Assembly mechanism/Fungi_all/Normalized/SeedVariable.ProcessImportance_EachBin_EachGroup.2.csv")
fun.phylo.bin.melt <- melt(fun.phylo.bin)
head(fun.phylo.bin.melt)
fun.phylo.bin.melt<-subset(fun.phylo.bin.melt, value != "NaN")
fun.phylo.bin.melt.mean <- fun.phylo.bin.melt %>% group_by(Index, variable) %>% summarise(MeanRelImportance = mean(value))
head(fun.phylo.bin.melt.mean)
names(fun.phylo.bin.melt.mean)[2] <- "Bin"
data3 <- fun.phylo.bin.melt.mean
data3$Bin <- factor(data3$Bin, levels = c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7",
                                          "Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14",
                                          "Bin15","Bin16","Bin17"))

data3$Index <- factor(data3$Index, levels = rev(c("HoS","HeS","DL","HD","DR")))

ggplot(data3, aes(x=Bin, y = MeanRelImportance, fill = Index)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("HoS"= "#3c93c2","HeS" = "#9ec9e2",
                               "DL" = "#fcde9c", "HD" = "#feb24c",
                               "DR"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 1.2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


### Assembly of core, rare, and other ASVs
##Bacteria
bac.phylo.bin2<-read.csv("./Assembly mechanism/Bacteria_all/Normalized/SeedVariable.iCAMP.Process_EachTurnover_EachCategory.2.csv")
bac.phylo.bin2.melt <- melt(bac.phylo.bin2)
head(bac.phylo.bin2.melt)
bac.phylo.bin2.melt<-subset(bac.phylo.bin2.melt, value != "NaN")
bac.phylo.bin2.melt$variable <- as.character(bac.phylo.bin2.melt$variable)

library(stringr)
bac.phylo.bin2.melt<-bac.phylo.bin2.melt %>% 
  mutate(Type=bac.phylo.bin2.melt$variable %>%      
           str_split("\\.") %>%  
           sapply(`[`,1), 
         Mechanism=bac.phylo.bin2.melt$variable %>%      
           str_split("\\.") %>%   
           sapply(`[`,2),
         .before=variable)   


bac.phylo.bin2.melt.mean <- bac.phylo.bin2.melt %>% group_by(Type, Mechanism) %>% summarise(MeanRelImportance = mean(value))
head(bac.phylo.bin2.melt.mean)

bac.phylo.bin2.melt.mean <- subset(bac.phylo.bin2.melt.mean, Mechanism != "Stochasticity")

bac.phylo.bin2.melt.mean$Mechanism <- factor(bac.phylo.bin2.melt.mean$Mechanism, levels = rev(c("HoS","HeS","DL","HD","DR")))

ggplot(bac.phylo.bin2.melt.mean, aes(x=Type, y = MeanRelImportance, fill = Mechanism)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("HoS"= "#3c93c2","HeS" = "#9ec9e2",
                               "DL" = "#fcde9c", "HD" = "#feb24c",
                               "DR"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


##Fungi
fun.phylo.bin2<-read.csv("./Assembly mechanism/Fungi_all/Normalized/SeedVariable.iCAMP.Process_EachTurnover_EachCategory.2.csv")
fun.phylo.bin2.melt <- melt(fun.phylo.bin2)
head(fun.phylo.bin2.melt)
fun.phylo.bin2.melt<-subset(fun.phylo.bin2.melt, value != "NaN")
fun.phylo.bin2.melt$variable <- as.character(fun.phylo.bin2.melt$variable)

fun.phylo.bin2.melt<-fun.phylo.bin2.melt %>% 
  mutate(Type=fun.phylo.bin2.melt$variable %>%      
           str_split("\\.") %>%  
           sapply(`[`,1), 
         Mechanism=fun.phylo.bin2.melt$variable %>%      
           str_split("\\.") %>%   
           sapply(`[`,2),
         .before=variable)   


fun.phylo.bin2.melt.mean <- fun.phylo.bin2.melt %>% group_by(Type, Mechanism) %>% summarise(MeanRelImportance = mean(value))
head(fun.phylo.bin2.melt.mean)

fun.phylo.bin2.melt.mean <- subset(fun.phylo.bin2.melt.mean, Mechanism != "Stochasticity")

fun.phylo.bin2.melt.mean$Mechanism <- factor(fun.phylo.bin2.melt.mean$Mechanism, levels = rev(c("HoS","HeS","DL","HD","DR")))

ggplot(fun.phylo.bin2.melt.mean, aes(x=Type, y = MeanRelImportance, fill = Mechanism)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("HoS"= "#3c93c2","HeS" = "#9ec9e2",
                               "DL" = "#fcde9c", "HD" = "#feb24c",
                               "DR"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


dev.off()

### Within and among branches
bac.branch.raw.data <- read.csv("./Assembly mechanism/Bacteria_all/Normalized/SeedVariable.ProcessImportance_EachGroup.csv")
bac.branch.raw.data$Category <- 0
bac.branch.raw.data$Category <- ifelse(bac.branch.raw.data$Group %in% c("Branch_2", "Branch_3","Branch_4","Branch_5","Branch_6","Branch_7","Branch_8"), "Within","Among")

##Normality test
shapiro.test(bac.branch.raw.data$HeS) # p-value = 0.07693 #normal distribution
shapiro.test(bac.branch.raw.data$HoS) # p-value = 0.997 #normal distribution
shapiro.test(bac.branch.raw.data$DL) # p-value = 0.3549 #normal distribution
shapiro.test(bac.branch.raw.data$HD) # p-value = 5.665e-08 #non-normal distribution
shapiro.test(bac.branch.raw.data$DR) # p-value = 0.2904 #normal distribution

bac.branch.raw.data.melt <- melt(bac.branch.raw.data)

mean(bac.branch.raw.data$HoS[which(bac.branch.raw.data$Category == "Within")]) #0.6214761
mean(bac.branch.raw.data$HoS[which(bac.branch.raw.data$Category == "Among")]) #0.5856083


### Homogeneous selection
bac.HoS <- subset(bac.branch.raw.data.melt, variable == "HoS")
bac.HoS$Category <- factor(bac.HoS$Category, levels = c("Within","Among"))

ggplot(bac.HoS, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Homogeneous selection
bac.HeS <- subset(bac.branch.raw.data.melt, variable == "HeS")
bac.HeS$Category <- factor(bac.HeS$Category, levels = c("Within","Among"))

ggplot(bac.HeS, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Dispersal limitation
bac.DL <- subset(bac.branch.raw.data.melt, variable == "DL")
bac.DL$Category <- factor(bac.DL$Category, levels = c("Within","Among"))

ggplot(bac.DL, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Homogenizing dispersal
bac.HD <- subset(bac.branch.raw.data.melt, variable == "HD")
bac.HD$Category <- factor(bac.HD$Category, levels = c("Within","Among"))

ggplot(bac.HD, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Drift
bac.DR <- subset(bac.branch.raw.data.melt, variable == "DR")
bac.DR$Category <- factor(bac.DR$Category, levels = c("Within","Among"))

ggplot(bac.DR, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


##Fungi
fun.branch.raw.data <- read.csv("./Assembly mechanism/Fungi_all/Normalized/SeedVariable.ProcessImportance_EachGroup.csv")
fun.branch.raw.data$Category <- 0
fun.branch.raw.data$Category <- ifelse(fun.branch.raw.data$Group %in% c("Branch_2", "Branch_3","Branch_4","Branch_5","Branch_6","Branch_7","Branch_8"), "Within","Among")

##Normality test
shapiro.test(fun.branch.raw.data$HeS) #  p-value = 1.351e-06 #non-normal distribution
shapiro.test(fun.branch.raw.data$HoS) # p-value = 0.5857 #normal distribution
shapiro.test(fun.branch.raw.data$DL) # p-value = 0.1842 #normal distribution
shapiro.test(fun.branch.raw.data$HD) # p-value = 0.003323 #non-normal distribution
shapiro.test(fun.branch.raw.data$DR) #p-value = 0.9984 #normal distribution

fun.branch.raw.data.melt <- melt(fun.branch.raw.data)

### Homogeneous selection
fun.HoS <- subset(fun.branch.raw.data.melt, variable == "HoS")
fun.HoS$Category <- factor(fun.HoS$Category, levels = c("Within","Among"))

ggplot(fun.HoS, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Homogeneous selection
fun.HeS <- subset(fun.branch.raw.data.melt, variable == "HeS")
fun.HeS$Category <- factor(fun.HeS$Category, levels = c("Within","Among"))

ggplot(fun.HeS, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Dispersal limitation
fun.DL <- subset(fun.branch.raw.data.melt, variable == "DL")
fun.DL$Category <- factor(fun.DL$Category, levels = c("Within","Among"))

ggplot(fun.DL, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Homogenizing dispersal
fun.HD <- subset(fun.branch.raw.data.melt, variable == "HD")
fun.HD$Category <- factor(fun.HD$Category, levels = c("Within","Among"))

ggplot(fun.HD, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### Drift
fun.DR <- subset(fun.branch.raw.data.melt, variable == "DR")
fun.DR$Category <- factor(fun.DR$Category, levels = c("Within","Among"))

ggplot(fun.DR, aes(x=Category, y= value, fill = Category))+
  geom_boxplot()+
  ggsignif::stat_signif(comparisons = list(c("Within","Among")), test = "t.test")+
  theme(aspect.ratio = 2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



####Bins
bacterial.bins$ASV <- 0
for (i in as.character(bacterial.bins$ID)){
  bacterial.bins$ASV[which(bacterial.bins$ID == i)] <- bac.list$OTU_id[which(bac.list$OTU == i)]
}
bacterial.bins

fungal.bins$ASV <- 0
for (i in as.character(fungal.bins$ID)){
  fungal.bins$ASV[which(fungal.bins$ID == i)] <- fun.list$OTU_id[which(fun.list$OTU == i)]
}
fungal.bins

write.csv(bacterial.bins,"220523_bacterial.bins.csv")
write.csv(fungal.bins,"220523_fungal.bins.csv")



## Core, rare, and others
RelImp.core.b.f <- read.xlsx("Ecological process_core and rare.xlsx",1)
ggplot(RelImp.core.b.f, aes(x=Class, y = RelImportance, fill = Category)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("HoS"= "#3c93c2","HeS" = "#9ec9e2",
                               "DL" = "#fcde9c", "HD" = "#feb24c",
                               "DR"="#fc4e2a")) +
  
  xlab('')+ theme(aspect.ratio = 2) + 
  ylab("Relative importance \n") +
  facet_wrap(~Kingdom)
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

  dev.off()
