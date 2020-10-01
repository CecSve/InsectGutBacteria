# Preparing data for analysis

setwd("H:/Documents/Speciale/New_DADA2_run_w_missing_samples")
# Load required libraries
library(tidyverse)
library(reshape2)
library(vegan)
library(lme4)
library(data.table)
library(ggplot2); theme_set(theme_bw(base_size = 14))
library(stringr)

# Load data
metadata <-
  read.table(
    "data/metadata_css_with_alltriplicates.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  ) #Loads your data
attach(metadata) #Makes it available to call columns

metadata_triplicates <- metadata[grep("^S",rownames(metadata)),]# make a dataframe for triplicate analysis
metadata_triplicates <- droplevels(metadata_triplicates)

metadata <- metadata[grep("^CSS",rownames(metadata)),] # subset the metadata to remove triplicates
metadata <- droplevels(metadata)

data <-
  read.table("data/dada2_seq_table_132.txt", sep = "\t", header = TRUE) #Loads your data with taxonomy attached
attach(data) #Makes it available to call columns

# ------ DATA WRANGLING -----------

# place the taxonomy in a seperate object
taxonomy <- data %>% select(X.OTUID, Taxonomy)
data <- data %>% select(-Taxonomy) # then remove the taxonomy column from the dataset

### ----- prepare otu data for triplicate analysis ------ 
data_triplicates <-
  data[, grep("^S", colnames(data))] # make a dataframe for triplicate analysis
data_triplicates$X.OTUID <- data$X.OTUID
data_triplicates <-
  column_to_rownames(data_triplicates, var = "X.OTUID")
setdiff(rownames(metadata_triplicates), colnames(data_triplicates)) # an extraction control did not pass the dada2 filtration and should be removed from metadata_triplicates
rowname_remove <- "S003"
metadata_triplicates <-
  metadata_triplicates[!(row.names(metadata_triplicates) %in% rowname_remove),]
metadata_triplicates <- droplevels(metadata_triplicates)

# examine wether there are reads in all samples and otus
colSums(data_triplicates) # all samples have reads
rowSums(data_triplicates) # not all otus have reads, since we removed the samples that had those otus associated

data_triplicates <- data_triplicates[rowSums(data_triplicates != 0) > 0, ] # removes samples/columns with no OTUs
data_triplicates <- droplevels(data_triplicates)
#names <- rownames(data_triplicates) # the seqvars that can be used for further analysis - only relevant if taxonomy is needed

# save output
write.table(data_triplicates, file = "cleaned-data/triplicate_otus.txt")
write.table(metadata_triplicates, file = "cleaned-data/triplicate_metadata.txt")

### -------- prepare otu data for gut/gut content analysis --------- 
# select data without triplicates
data <- column_to_rownames(data, var = "X.OTUID")
data <- data[,grep("^CSS",colnames(data))]

# these samples did not pass the dada2 filtration
setdiff(rownames(metadata), colnames(data)) 
#metadata %>% count(Stratum) # some of the extra samples are negatives that should be removed from metadata and some samples did not pass the DADA2 filter
metadata <- rownames_to_column(metadata, var = "otuid") # in order to not lose rownames, we transfer them into a new column
#metadata <- metadata %>% filter(Sample != 'Negative') # removing all extraction blanks and PCR negatives

# make a dataset for the other stratum that will not be examined further for now
metadata_otherstratum <-
  metadata %>% filter(
    Stratum == 'Eggs' |
      Stratum == 'Mites' |
      Stratum == 'Reproductive_organs' |
      Stratum == 'Reproductive_organs_eggs'
  )

# remove the other stratum
metadata <-
  metadata %>% filter(Stratum != 'Eggs') %>% filter(Stratum != 'Mites') %>% filter(Stratum != 'Reproductive_organs') %>% filter(Stratum != 'Reproductive_organs_eggs') # Removing egg, reproductive organes and mite samples

metadata <-
  column_to_rownames(metadata, var = "otuid") # remake the rownmaes
samples_remove <-
  setdiff(rownames(metadata), colnames(data)) # examine which samples are missing from OTU data that are present in metadata. These samples did not pass through the DADA2 pipeline

# the samples that didn't pass DADA2 needs to be removed from metadata
metadata <- rownames_to_column(metadata, var = "sampleid")
metadata <- metadata %>% filter(!sampleid %in% samples_remove)
metadata <- column_to_rownames(metadata, var = "sampleid")

samples_remove <- setdiff(colnames(data), rownames(metadata)) # samples present in OTU data but missing from metdadata. We just removed a bunch of negatives, so these samples represent those negatives (double checked with original metadata)

# remove the negatives based on their name from the otu table
data <- data[ , -which(names(data) %in% samples_remove)] # data and metadata now have equal sample sizes

# Wrangle the taxonomy
taxonomy <- column_to_rownames(taxonomy, var = "X.OTUID") # set the seqvar names as rownames
taxa <- str_split_fixed(taxonomy$Taxonomy, ";", 7) # split the taxonomy string into all taxonomic levels
taxa <- as.data.frame(taxa) # coerce into dataframe
rownames(taxa) <- rownames(taxonomy) # reset the rownames/seqvar 
colnames(taxa) <-
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# examine wether there are reads in all samples and otus
colSums(data) # all samples have reads
rowSums(data) # not all otus have reads, since we removed the samples that had those otus associated

data <- data[rowSums(data != 0) > 0, ] # removes samples/columns with no OTUs
names <- rownames(data) # the seqvars that can be used for further analysis

# match seqvars between data and taxonomy
taxonomy$seqvars <- row.names(taxonomy) # create a variable in taxonomy for seqvars
taxonomy <- taxonomy[match(names, rownames(taxonomy)), ] # make sure the SeqVars in the data and the taxonomy match
taxonomy <- taxonomy %>% select(-seqvars) # remove seqvars column after filtration

taxa$seqvars <- row.names(taxa)
taxa <- taxa[match(names, rownames(taxa)), ] # make sure the SeqVars in the data and the taxonomy match
taxa <- taxa %>% select(-seqvars) # remove seqvars column after filtration

data <- droplevels(data)
metadata <- droplevels(metadata)
taxonomy <- droplevels(taxonomy)
taxa <- droplevels(taxa)

# save outputs
write.table(data, file = "cleaned-data/otus.txt")
write.table(metadata, file = "cleaned-data/metadata.txt")
write.table(taxonomy, file = "cleaned-data/taxonomy_string.txt") # string format is required for phyloseq
write.table(taxa, file = "cleaned-data/taxonomy_split.txt")
