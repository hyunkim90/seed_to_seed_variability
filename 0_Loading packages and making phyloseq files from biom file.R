## import libraries : phyloseq and microbiome
#remotes::install_github("vmikk/metagMisc")

library(dplyr)
library(forcats) 
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(scales)
library(grid)
library(reshape2)
library(seqtime)
library(agricolae)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
library(multifunc)
library(FSA)
library(rcompanion)
library(seqinr)
library(metagMisc)
# Set plotting theme
theme_set(theme_bw())


####### For 70 seed samples
### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("./Bacteria/asv_table_final.biom")

### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
map <- read.table(file = './Bacteria/sample_metadata.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID
rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("./Bacteria/tree.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)
phy   ## 864 ASVs

## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy  ## 18811 ASVs

sort(colSums(otu_table(phy)))

## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("./Fungi/Latest DB/asv_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = './Fungi/sample_metadata.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("./Fungi/Latest DB/tree.nwk") #rooted tree
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 395 ASVs

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun
summarize_phyloseq(fun)
sort(colSums(otu_table(fun)))



####### For Time series seed samples
### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("./Time series data/bacteria/asv_table_final.biom")

### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
map <- read.table(file = './Time series data/bacteria/sample_metadata.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID
rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("./Time series data/bacteria/tree.nwk")
phy.time <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy.time)
phy.time   ## 368 ASVs

## changing rank names
colnames(tax_table(phy.time)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy.time  ## 368 ASVs

sort(colSums(otu_table(phy.time)))

## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("./Time series data/fungi/dynamic DB_time series/asv_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = './Time series data/fungi/dynamic DB_time series/sample_metadata.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("./Time series data/fungi/dynamic DB_time series/tree.nwk") #rooted tree
fun.time <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun.time)
fun.time   ## 372 ASVs

## changing rank names
colnames(tax_table(fun.time)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun
summarize_phyloseq(fun.time)
sort(colSums(otu_table(fun.time)))

