### Normalization

## Let's do normalization with CSS
## phyloseq to metagenomeSeq


#phy.clean? or phy.clean --> let's start from phy.clean

bac.clean.ss
fun.clean.ss

## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(bac.clean.ss) == 0)
taxa_sums(bac.clean.ss)

bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(bac.clean.ss) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.clean.ss))
sort(sample_sums(fun.clean.ss))
sort(sample_sums(bac.time.clean.ss))
sort(sample_sums(fun.time.clean.ss))
## only keep samples over 0 read
(filt.sample <- sample_sums(bac.clean.ss) > 0)
sum(sample_sums(bac.clean.ss) <= 0)  ## 1 sample discarded
bac.clean.ss.f <- prune_samples(filt.sample, bac.clean.ss)
bac.clean.ss.f  


sort(sample_sums(bac.clean.ss.f))

## CODE for CSS normalization using preloaded data
bac.clean.filt <- bac.clean.ss.f
bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog <- bac.clean.filt
otu_table(bac.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log <- bac.clean.filt
otu_table(bac.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

### fungal community ###
## only keep samples over 0 read
(filt.sample <- sample_sums(fun.clean.ss) > 0)
sum(sample_sums(fun.clean.ss) <= 10)  ## 1 sample discarded
fun.clean.ss.f <- prune_samples(filt.sample, fun.clean.ss)
fun.clean.ss.f  ## 979 samples <- 984 samples


sort(sample_sums(fun.clean.ss.f))

## CODE for CSS normalization using preloaded data
fun.clean.filt <- fun.clean.ss.f
fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog <- fun.clean.filt
otu_table(fun.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log <- fun.clean.filt
otu_table(fun.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

###Time series
#phy.clean? or phy.clean --> let's start from phy.clean

bac.time.clean.ss
fun.time.clean.ss

## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(bac.time.clean.ss) == 0)
taxa_sums(bac.time.clean.ss)

bac.time.clean.ss <- phyloseq::filter_taxa(bac.time.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(bac.time.clean.ss) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.time.clean.ss))

mean(sample_sums(fun.time.clean.ss))
sum(sample_sums(fun.time.clean.ss))
## only keep samples over 0 read
(filt.sample <- sample_sums(bac.time.clean.ss) > 0)
sum(sample_sums(bac.time.clean.ss) <= 0)  ## 1 sample discarded
bac.time.clean.ss.f <- prune_samples(filt.sample, bac.time.clean.ss)
bac.time.clean.ss.f  


sort(sample_sums(bac.time.clean.ss.f))

## CODE for CSS normalization using preloaded data
bac.time.clean.filt <- bac.time.clean.ss.f
bac.time.clean.filt <- phyloseq::filter_taxa(bac.time.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.time.clean.filt    ## use all samples

met.bac.time.clean <- phyloseq_to_metagenomeSeq(bac.time.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.time.clean)
p
met.bac.time.norm <- cumNorm(met.bac.time.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.time.norm)
sort(normFactors(met.bac.time.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.time.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.time.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.time.clean.nolog <- bac.time.clean.filt
otu_table(bac.time.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.time.clean.log <- bac.time.clean.filt
otu_table(bac.time.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

### fungal community ###
## only keep samples over 0 read
(filt.sample <- sample_sums(fun.time.clean.ss) > 0)
sum(sample_sums(fun.time.clean.ss) <= 10)  ## 1 sample discarded
fun.time.clean.ss.f <- prune_samples(filt.sample, fun.time.clean.ss)
fun.time.clean.ss.f  ## 979 samples <- 984 samples


sort(sample_sums(fun.time.clean.ss.f))

## CODE for CSS normalization using preloaded data
fun.time.clean.filt <- fun.time.clean.ss.f
fun.time.clean.filt <- phyloseq::filter_taxa(fun.time.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.time.clean.filt    ## use all samples

met.fun.time.clean <- phyloseq_to_metagenomeSeq(fun.time.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.time.clean)
p
met.fun.time.norm <- cumNorm(met.fun.time.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.time.norm)
sort(normFactors(met.fun.time.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.time.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.time.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.time.clean.nolog <- fun.time.clean.filt
otu_table(fun.time.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.time.clean.log <- fun.time.clean.filt
otu_table(fun.time.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



#### Rarefying
set.seed(8046)

bac.rarefied<-rarefy_even_depth(bac.clean.ss.f, sample.size = min(sample_sums(bac.clean.ss.f)),
                  rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

fun.rarefied<-rarefy_even_depth(fun.clean.ss.f, sample.size = min(sample_sums(fun.clean.ss.f)),
                                rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


bac.time.rarefied<-rarefy_even_depth(bac.time.clean.ss.f, sample.size = min(sample_sums(bac.time.clean.ss.f)),
                                rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

fun.time.rarefied<-rarefy_even_depth(fun.time.clean.ss.f, sample.size = min(sample_sums(fun.time.clean.ss.f)),
                                rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



