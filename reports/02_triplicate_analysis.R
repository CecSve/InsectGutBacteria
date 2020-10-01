# Triplicate analysis

# Load required libraries
library(tidyverse)
library(reshape2)
library(vegan)
library(lme4)
library(data.table)
library(ggplot2); theme_set(theme_bw(base_size = 14))
library(stringr)
library(DESeq2)
library(phyloseq)
library(funrar)
library(ggpubr)
library(ggrepel)
library(cowplot)
#library(MASS)

# set working directory
#setwd("H:/Documents/Speciale/New_DADA2_run_w_missing_samples")

# load data
data <-
  read.table(
    "cleaned-data/triplicate_otus.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# load metadata
metadata <-
  read.table(
    "cleaned-data/triplicate_metadata.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

### examine negatives #####
negs <- metadata %>% rownames_to_column(var = "pcrid") %>% filter(Sample == "Negative")
keep <- negs$pcrid
negs_asvs <- data[, (colnames(data) %in% keep)]
max(colSums(negs_asvs))
max(colSums(data))

mean(colSums(negs_asvs))
mean(colSums(data))

# subsetting data further by removing negatives
names <- metadata %>% filter(Sample == "Negative")
names <- droplevels(names)
names <- names$SpecimenID # save the samples names of the negaties for subsetting metadata further 

metadata <-
  metadata[!(metadata$SpecimenID %in% names),] # remove negatives from metadata
metadata <- droplevels(metadata)
metadata <- metadata %>% rownames_to_column(var = "pcrid") %>% filter(Stratum == "Gut" | Stratum == "Fluid") %>% column_to_rownames(var = "pcrid")
attach(metadata)
table(metadata$SpecimenID)

names <- rownames(metadata) # save sample names to keep
data <- data[, (colnames(data) %in% names)] # keep only samples without negatives
data <- droplevels(data)
rowSums(data)
mean(colSums(data))
# remove otus that only contain zeros
data <-
  data[apply(data[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)
attach(data)

## ------- Analysis --------
#Load data, taxnomy file and metadata into phyloseq formats
data_one <- data + 1
OTU <- otu_table(data_one, taxa_are_rows = TRUE) # changed taxa_are_rows to true for phyloseq to work
#TAX <- tax_table(taxonomy)
sampledata <- sample_data(metadata)

#Create create phyloseq object named physeq
physeq <- phyloseq(OTU, sampledata)
physeq <- merge_phyloseq(physeq, metadata)

#Convert your phyloseq object to deseq2 object. ~ is the experimental design formula. Can be more complicated. 
diagdds = phyloseq_to_deseq2(physeq, ~ SpecimenID) # have set specimenID since triplicates are based on those
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds, fitType = "parametric")
diagvst = getVarianceStabilizedData(diagdds)

physeq0 = physeq #This keeps your original physeq object before normalization
otu_table(physeq) <- otu_table(diagvst, taxa_are_rows = TRUE) #Replaces your OTU data in physeq with normalised data.

#Creates a data.frame with normalised data so you don't need to use phyloseq.
norm.data = t(as(otu_table(physeq), "matrix")) # notice it is tranversed
dim(norm.data)
min(norm.data)
max(norm.data)
median(norm.data)
norm.data[norm.data<0.875] <- 0 # something is strange with the zeros having values... setting median as min

rowSums(norm.data)
min(rowSums(norm.data))

#This command rarefies data to your lowest sampling depth, you can save this in phyloseq instead of the DeSeq2 data
rarified <- rarefy_even_depth(OTU, sample.size = min(sample_sums(OTU)), replace = FALSE, trimOTUs = FALSE, rngseed = 711) # set to 5000 but could be changed. Following samples were removed: S006 S026 S032 S046 and should be removed from metadata as well for richness analysis - now sat to min sample depth
rarefied.data = t(as(otu_table(rarified), "matrix")) 
rarefied.data <- t(rarefied.data)
rarefied.data <- as.data.frame(rarefied.data)

colSums(rarefied.data)
rowSums(rarefied.data) # both look right

############Convert data to relative abundances################# 
norm.data <- as.matrix(norm.data)
rel.data <- (make_relative(norm.data)*100)
dim(rel.data)

### nMDS: community composition ###########################
bray.data<-vegdist(norm.data, Type = "bray", binary = FALSE) 
meta_mds <- metaMDS(bray.data)
plot(meta_mds)

class(meta_mds)
methods(class="metaMDS")

Y <- scores(meta_mds, display="sites")
plot(Y, type="n")
text(Y[,1], Y[,2], rownames(Y), col="red")

all(rownames(norm.data) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data
# data for samples with read count over 50 (33 samples)
Y.sel <- Y[rowSums(norm.data) > 50, ] # select ones with abundance over 50
plot(Y.sel[,1], Y.sel[,2], type="n", xlim=c(-.8, .8), ylim=c(-.4, .4))
text(Y.sel[,1], Y.sel[,2], rownames(Y.sel), col="red")

Y <- data.frame(Y)
Y$abundance <- rowSums(norm.data)
Y$labels <- metadata$Stratum
Y$colour <- metadata$SpecimenID

plotcols <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour friendly

# plot all samples
ggplot(Y, aes(x=NMDS1, y=NMDS2, size=abundance, colour = colour, group = colour)) +
  geom_point() + geom_text_repel(aes(label=labels), show.legend = FALSE) + 
  theme_minimal() + scale_size_continuous(range = c(5,10)) + labs(colour = "Sample ID", size = "Relative ASV \nabundance") + geom_line(colour = "lightgrey", size = 1, show.legend = F) + scale_colour_manual(values = plotcols)

trip_plot <- ggplot(Y, aes(x=NMDS1, y=NMDS2, size=abundance, colour = colour, group = colour)) +
  geom_point(alpha=0.7, shape = 16) + geom_text_repel(aes(label=labels), show.legend = FALSE) + 
  theme_minimal() + scale_size_continuous(range = c(5,10)) + labs(colour = "Sample ID", size = "Relative ASV \nabundance") + geom_line(size = 1, alpha = 0.5, show.legend = F) + scale_colour_manual(values = plotcols)

save_plot("triplicate-plots/triplicate_nmds_18082020.png", trip_plot, base_height = 8, base_width = 12)

### PERMANOVA #################

# # Homogeneity of dispersion test - betadiversity
beta <- betadisper(bray.data, metadata$SpecimenID, bias.adjust = T) # ## Calculate multivariate dispersions
plot(beta)
boxplot(beta) ## Draw a boxplot of the distances to centroid for each group
beta.HSD <- TukeyHSD(beta) ## Tukey's Honest Significant Differences
plot(beta.HSD)
anova(beta) # so not equal variances but that's what we assume anyways - we assume the community composition variances
permutest(beta, permutations = 999) ## Permutation test for F
# our betadisper results are significant, meaning we cannot accept the null hypothesis that our groups have the same dispersions (homogeneity of dispersion among groups (specimenID in our case) is a condition (assumption) for adonis). This means we can be less confident that our adonis result is a real result, because heterogeneity in group dispersions is not clear

beta <- betadisper(bray.data, metadata$Stratum, bias.adjust = T) # ## Calculate multivariate dispersions
plot(beta)
boxplot(beta) ## Draw a boxplot of the distances to centroid for each group
beta.HSD <- TukeyHSD(beta) ## Tukey's Honest Significant Differences
plot(beta.HSD)
anova(beta) # so not equal variances but that's what we assume anyways - we assume the community composition variances
permutest(beta, permutations = 999) ## Permutation test for F

# Adonis test
adonis(bray.data ~ SpecimenID + Stratum, data = metadata, strata = metadata$Species) # if their dispersions are quite different, adonis will give you a significant p-value, thus, the result is heavily influenced not by the difference in composition between groups but by differences in composition within groups (heterogeneous dispersion, and thus a measure of betadiversity).


# 1) remove low sample reads --> rrarefy to something  --> (you may not need to remove the low ones because they will be removed automatically) --> bray --> mds 2) chris sends a new script (the right one) and you should run this first

# Extract one of the triplicate samples each based on highest read count
rich <- as.data.frame(t(data))
rich$readcount <- rowSums(rich)
test <- rich %>% select(readcount)
rich <- rownames_to_column(test, var = "sample")
met <- rownames_to_column(metadata, var = "sample")
met <- met %>% select(sample, SpecimenID, Stratum)
rich <- left_join(met, rich, by = "sample")
top_trip <- rich %>% group_by(SpecimenID,sample) %>% summarise(readcount = max(readcount)) %>% group_by(SpecimenID) %>% dplyr::slice(which.max(readcount))  # find out which triplicate sample has the highest readcount and save them

samples <- top_trip$sample
data <- data[, (colnames(data) %in% samples)]
data <- droplevels(data)

# save output
write.table(top_trip, file = "cleaned-data/triplicate_otus_highest_readcount.txt")
write.table(data, file = "cleaned-data/triplicate_otus_filt.highest_readcount_otutable.txt")
