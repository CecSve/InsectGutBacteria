# Microbiome analysis

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
library(ggrepel)
library(cowplot)
library(hrbrthemes)
library(grid)
library(gridExtra)
#library(mvabund)
#library(MASS)

# set working directory
#setwd("H:/Documents/Speciale/New_DADA2_run_w_missing_samples")

# load data
data <-
  read.table(
    "cleaned-data/otus.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# load triplicate data with highest reads to be included in the analysis
trip_data <-
  read.table(
    "cleaned-data/triplicate_otus_filt.highest_readcount_otutable.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# load metadata
metadata <-
  read.table(
    "cleaned-data/metadata.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# load triplicate metadata
trip_metadata <-
  read.table(
    "cleaned-data/triplicate_metadata.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# select only the samples that are retained in the triplicate otu table
trip_metadata <-
  trip_metadata[(rownames(trip_metadata) %in% names(trip_data)),]
trip_metadata <-
  trip_metadata %>% rownames_to_column(var = "sample") %>% filter(!Stratum == "Eggs") %>% column_to_rownames(var = "sample") # remove eggs from the analysis

trip_data <-
  trip_data[, (names(trip_data) %in% rownames(trip_metadata))] # remove egg samples from otu data

# load split taxonmy
taxa <- 
  read.table(
    "cleaned-data/taxonomy_split.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# load string taxonomy
taxonomy <- 
  read.table(
    "cleaned-data/taxonomy_string.txt",
    sep = "",
    header = TRUE,
    row.names = 1
  ) 

# merge triplicate otutable with the other otutable
data <- rownames_to_column(data, var = "seqvars")
trip_data <- rownames_to_column(trip_data, var = "seqvars") 

comb_data <- left_join(data, trip_data, by = "seqvars") # be careful - if there's not the same seqvars in both tables, then unmatching otus will get lost when left_joining 
comb_data <- column_to_rownames(comb_data, var = "seqvars")

# merge metadata with triplicate metadata
metadata <- rownames_to_column(metadata, var = "sample")
trip_metadata <- rownames_to_column(trip_metadata, var = "sample")

comb_metadata <- rbind(metadata, trip_metadata) 
comb_metadata <- column_to_rownames(comb_metadata, var = "sample")

data <- droplevels(comb_data)
metadata <- droplevels(comb_metadata)
taxonomy <- droplevels(taxonomy)

# one outlier occurs in the mds, so that sample will be excluded
data <- data[, !colnames(data) %in% "CSS141"] # this sample is a gut content/fluid sample so it wil automatically be removed because only guts will be used for this analysis
metadata <-  metadata[!rownames(metadata) %in% "CSS141",]

data <- as.data.frame(data)
data[is.na(data)] <- 0 # NAs were introduced with triplicate samples and needs to be changed to zeros

# remove otus that only contain zeros
data <-
  data[apply(data[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

egg_otus <- setdiff(rownames(comb_data), rownames(data)) # save the otus that were only present in eggs if you want to investigate further, e.g. their taxonomy
egg_taxonomy <- taxa %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% egg_otus) %>% column_to_rownames(var = "seqvars")

occurences <-
  rownames(data) # save occurence IDs so they can be used to subset taxonomy

taxonomy <- taxonomy %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars")
taxa <- taxa %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars")

# only keep gut samples for analysis
metadata <- metadata %>% rownames_to_column(var = "sample") %>% filter(!Stratum == "Fluid") %>% column_to_rownames(var = "sample")
data <- data[, (names(data) %in% rownames(metadata))] # remove fluid samples from otu data
rowSums(data) 

# remove otus that only contain zeros
data <-data[apply(data[,-1], 1, function(x)! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

occurences <-rownames(data) # save occurence IDs so they can be used to subset taxonomy

taxonomy <- taxonomy %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars")
taxa <- taxa %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars")

### Analysis ###############
plyr::count(metadata$Species) 

qqnorm(log(colSums(data)))

# only keep gut samples from +9 sampled specimens
#metadata_res <-metadata %>% rownames_to_column(var = "sample") %>% filter(!Species %in% c("calcarata", "equestris", "limbatus")) %>%column_to_rownames(var = "sample")

#data_res <-data[, (names(data) %in% rownames(metadata_res))] # remove fluid samples from otu data
#rowSums(data_res) 

# remove otus that only contain zeros
#data_res <-data_res[apply(data_res[,-1], 1, function(x) ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

occurences <-rownames(data) # save occurence IDs so they can be used to subset taxonomy

taxa_res <-taxa %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% filter(!Kingdom %in% c("Eukaryota", "Archaea")) %>% column_to_rownames(var = "seqvars")

occurences <-rownames(taxa_res) # save occurence IDs so they can be used to subset taxonomy

taxonomy_res <- taxonomy %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars")

data_res <- data %>% rownames_to_column(var = "seqvars") %>% filter(seqvars %in% occurences) %>% column_to_rownames(var = "seqvars") # remove OTUs assigned to Eukaryota, Archaea from the otu table

data_res <- droplevels(data_res)
metadata_res <- droplevels(metadata)
taxa_res <- droplevels(taxa_res)
taxonomy_res <- droplevels(taxonomy_res)
  
### DeSeq normalization ######################
#Load data, taxnomy file and metadata into phyloseq formats
data_one <- data_res + 1
OTU <- otu_table(data_one, taxa_are_rows = T) # changed taxa_are_rows to true for phyloseq to work
sampledata <- sample_data(metadata_res)
taxa <- as.matrix(taxonomy_res) #Makes a dataframe
TAX <- tax_table(taxa)

#Create create phyloseq object named physeq
physeq <- phyloseq(OTU, sampledata, TAX)
physeq <- merge_phyloseq(physeq, metadata)

#Convert your phyloseq object to deseq2 object. ~ is the experimental design formula. Can be more complicated. 
diagdds = phyloseq_to_deseq2(physeq, ~ Species) # set what you wish to test
diagdds = DESeq(diagdds, test = "Wald", fitType = "local")
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds, fitType = "local")
diagvst = getVarianceStabilizedData(diagdds)

physeq0 = physeq #This keeps your original physeq object before normalization
otu_table(physeq) <- otu_table(diagvst, taxa_are_rows = TRUE) #Replaces your OTU data in physeq with normalised data.

#Creates a data.frame with normalised data so you don't need to use phyloseq.
norm.data = t(as(otu_table(physeq), "matrix")) # notice it is tranversed
dim(norm.data)
min(norm.data)
norm.data[norm.data<0] <- 0 # set negatives to zeros

rowSums(norm.data)
min(rowSums(norm.data))

#This command rarefies data to your lowest sampling depth, you can save this in phyloseq instead of the DeSeq2 data
rarified <- rarefy_even_depth(OTU, sample.size = min(sample_sums(OTU)), replace = FALSE, trimOTUs = FALSE, rngseed = 711) # set to 5000 but could be changed. Following samples were removed: S006 S026 S032 S046 and should be removed from metadata as well for richness analysis - now set to min sample depth
rarefied.data = t(as(otu_table(rarified), "matrix")) 
rarefied.data <- t(rarefied.data)
rarefied.data <- as.data.frame(rarefied.data)

colSums(rarefied.data)
rowSums(rarefied.data) # both look right

############Convert data to relative abundances#################  
norm.data <- as.matrix(norm.data)
rel.data <- (make_relative(norm.data)*100)
dim(rel.data)

# ----- chris' method - using rel.data ----------
bray.data<-vegdist(rel.data, Type = "bray", binary = F) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
mds<-monoMDS(bray.data) # Performs multidimensional scaling for the data
plot(mds)

Nmds1<-mds$points[,1] #Extracts dimension 1 as a vector
Nmds2<-mds$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

Nmds=cbind(Nmds1,Nmds2) #Makes a new sheet with dimension1 and 2
Nmds<-as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

sample.MDS <-ggplot(data=Nmds, aes(y=Nmds2, x=Nmds1, Type="p", colour = metadata_res$Species)) + 
  geom_point(size = 4) + labs(x="Dimension 1", y = "Dimension 2") + 
  theme(plot.title=element_text(hjust=0)) + stat_ellipse(aes(colour = metadata_res$Species), show.legend = F)
sample.MDS #Shows you our graph

# ----my method - using norm. data -----
bray.data<-vegdist(rel.data, Type = "bray", binary = F) 
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

# histogram of distance differences
# merge distance data and metadata
mdsdata <- Y %>% rownames_to_column(var = "PCRID")
meta <- metadata_res %>% rownames_to_column(var = "PCRID")

data <- merge(meta, mdsdata, by = "PCRID")

# add differences in distance column
data$distdif <- data$NMDS1 - data$NMDS2

#lmer on distdiff
library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

### model of community compostition: species level ############

distmodel <- lmer(distdif~ Species + (1|SpecimenID), data = data)
summary(distmodel)
drop1(distmodel)
r.squaredGLMM(distmodel)

ems <- emmeans(distmodel, "Species",  type = "response")
pairs(ems)

ems_contrasts = pairs(ems) %>%
  as.data.frame()

write.table(ems_contrasts,file="cleaned-data/ems_contrast_species_comcomp.txt",sep="\t",col.names=NA, row.names = T)

summary(ems)
plot(ems, comparisons = TRUE)

ems_means <- summary(ems) %>%
  as.data.frame()

write.table(ems_means,file="cleaned-data/ems_means_species_comcomp.txt",sep="\t",col.names=NA, row.names = T)

plotcols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour friendly

### species plot: emmeans #################
p <-
  ems_means %>% mutate(
    Species = fct_relevel(
      Species,
      "griseoaptera",
      "cyathigerum",
      "prasina",
      "cinerea",
      "equestris",
      "calcarata",
      "limbatus"
    )
  ) %>% ggplot(aes(Species, emmean))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

speciesplot <-
  d + geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE, colour = Species),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "Species", y = "", subtitle = "A") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                        labels = c(
                          "griseoaptera" = "Ph. griseoaptera",
                          "cyathigerum" = "E. cyathigerum",
                          "prasina" = "Pa. prasina",
                          "cinerea" = "Ne. cinerea",
                          "equestris" = "L. equestris",
                          "calcarata" = "S. calcarata",
                          "limbatus" = "Na. limbatus"
                        )
                      ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols)

# plot black and white
#this defines the elements to go in the plot, both the x and y and upper and lower CIs
a<- ems_means %>% mutate(
  Species = fct_relevel(
    Species,
    "griseoaptera",
    "cyathigerum",
    "prasina",
    "limbatus",
    "calcarata",
    "equestris",
    "cinerea"
  )
) %>% ggplot(aes(x=Species,y=emmean,ymax=upper.CL,ymin=lower.CL,size=2))
#this defines the plot type
b<-a+geom_pointrange()
#this flips the co-ordinates so your x axis becomes your y and vice versa
c<-b+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)
#all of this gets rid of the grey grid and legends
e<-d+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
#this sets x and y axis titles
f<-e+ labs(x = "Species\n", y = "", subtitle = "A")
#this sets axis label size
g<-f+theme(axis.text.x = element_text(size = 12, colour = 'black')) +theme(axis.text.y = element_text(size = 12, colour = 'black', face = "italic"))
#this sets axis title size and there is your finished summary plot!
g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
    "griseoaptera" = "Ph. griseoaptera",
    "cyathigerum" = "E. cyathigerum",
    "prasina" = "Pa. prasina",
    "cinerea" = "N. cinerea",
    "equestris" = "L. equestris",
    "calcarata" = "S. calcarata",
    "limbatus" = "Na. limbatus"
  )) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

p1 <- g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "griseoaptera" = "Ph. griseoaptera",
  "cyathigerum" = "E. cyathigerum",
  "prasina" = "Pa. prasina",
  "cinerea" = "N. cinerea",
  "equestris" = "L. equestris",
  "calcarata" = "S. calcarata",
  "limbatus" = "Na. limbatus"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

### correlation plot ##############
# there appears to be some correlation between the diet and habitat variables, and the correlation variables should be removed from the diet and habitat analysis
library(corrplot)
test <- data %>% select(PCRID, Habitat, Diet, distdif)
str(test)
test1 <- spread(test, Habitat, distdif, fill = 0)
test2 <- spread(test, Diet, distdif, fill = 0)
test3 <- merge(test1, test2, by = "PCRID")
test4 <- test3[, c(1, 4:7, 9:11)]

rownames(test4) <- NULL
test <- test4 %>% column_to_rownames(var = "PCRID")
p <- cor(test)

# add significance
res1 <- cor.mtest(test, conf.level = .95, na.action("na.omit"))
res2 <- cor.mtest(test, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = plotcols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE)

# so grasshoppers and scorpion bugs needs to be removed due to high correlation (the only species with that feeding mode and habitat combined - could be solved by adding species that fit in both categories)

### model of community compostition: diet level ############
environ_data <- data %>% filter(Species %in% c("calcarata", "cyathigerum", "limbatus", "prasina", "equestris"))

distmodel <- lmer(distdif~ Diet + (1|Species), data = environ_data) # species as random effect
summary(distmodel)
drop1(distmodel)
r.squaredGLMM(distmodel)

ems <- emmeans(distmodel, "Diet",  type = "response")
pairs(ems)

ems_contrasts = pairs(ems) %>%
  #confint() %>%
  as.data.frame()

write.table(ems_contrasts,file="cleaned-data/ems_contrast_diet_comcomp.txt",sep="\t",col.names=NA, row.names = T)

plot(ems, comparisons = TRUE)

ems_means <- summary(ems) %>%
  as.data.frame()

write.table(ems_means,file="cleaned-data/ems_means_diet_comcomp.txt",sep="\t",col.names=NA, row.names = T)

plotcols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour friendly

### diet plot: emmeans ########## 
p <-
  ems_means %>% mutate(
    Diet = fct_relevel(
      Diet,
      "herbivore",
      "predator"
    )
  ) %>% ggplot(aes(Diet, emmean))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

dietplot <-
  d + geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE, colour = Diet),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "Diet", y = "", subtitle = "B") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                        labels = c(
                          "herbivore" = "Herbivore",
                          "predator" = "Predator"
                        )
                      ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols[c(3,6)])

# plot black and white
#this defines the elements to go in the plot, both the x and y and upper and lower CIs
a<- ems_means %>% ggplot(aes(x=Diet,y=emmean,ymax=upper.CL,ymin=lower.CL,size=2))
#this defines the plot type
b<-a+geom_pointrange()
#this flips the co-ordinates so your x axis becomes your y and vice versa
c<-b+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)
#all of this gets rid of the grey grid and legends
e<-d+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
#this sets x and y axis titles
f<-e+ labs(x = "Diet\n", y = "", subtitle = "B")
#this sets axis label size
g<-f+theme(axis.text.x = element_text(size = 12, colour = 'black')) +theme(axis.text.y = element_text(size = 12, colour = 'black'))
#this sets axis title size and there is your finished summary plot!
g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "predator" = "Predator",
  "omnivore" = "Omnivore",
  "herbivore" = "Herbivore"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

p2 <- g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "predator" = "Predator",
  "omnivore" = "Omnivore",
  "herbivore" = "Herbivore"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

### model of community composition: habitat level ############
table(data$Habitat)
environ_data <- data %>% filter(Species %in% c("calcarata", "cyathigerum", "limbatus", "prasina", "equestris"))

distmodel <- lmer(distdif~ Habitat + (1|Diet), data = environ_data) # ecosystem as random effect
summary(distmodel)
drop1(distmodel)
r.squaredGLMM(distmodel)

ems <- emmeans(distmodel, "Habitat",  type = "response")
pairs(ems)

ems_contrasts = pairs(ems) %>%
  #confint() %>%
  as.data.frame()

write.table(ems_contrasts,file="cleaned-data/ems_contrast_habitat_comcomp.txt",sep="\t",col.names=NA, row.names = T)

plot(ems, comparisons = TRUE)

ems_means <- summary(ems) %>%
  as.data.frame()

write.table(ems_means,file="cleaned-data/ems_means_habitat_comcomp.txt",sep="\t",col.names=NA, row.names = T)

#eff_size(ems, sigma = sigma(distmodel), edf = Inf)
#eff_size(pairs(ems), sigma = sigma(distmodel), edf = Inf, method = "identity")
#pwpp(ems, type = "response")

plot(ems, comparisons = TRUE)
plot(ems)

plotcols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour friendly

### habitat plot: emmeans ########## 
table(environ_data$Habitat)

p <-
  ems_means %>% mutate(
    Habitat = fct_relevel(
      Habitat,
      "Vincetoxicum_hirundinaria",
      "Grassland_Freshwater",
      "Grassland",
      "Aegopodium_podagrarium"
    )
  ) %>% ggplot(aes(Habitat, emmean))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

habitatplot <-
  d + geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE, colour = Habitat),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "Habitat", y = "", subtitle = "C") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                        labels = c(
                          "Aegopodium_podagrarium" = "Aegopodium podagrarium",
                          "Grassland" = "Grassland",
                          "Grassland_Freshwater" = "Freshwater (nymph)\nGrassland (adult)",
                          "Vincetoxicum_hirundinaria" = "Vincetoxicum hirundinaria"
                        )
                      ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols[c(5,2,3,1)])

# plot black and white
#this defines the elements to go in the plot, both the x and y and upper and lower CIs
a<- ems_means %>% ggplot(aes(x=Habitat,y=emmean,ymax=upper.CL,ymin=lower.CL,size=2))
#this defines the plot type
b<-a+geom_pointrange()
#this flips the co-ordinates so your x axis becomes your y and vice versa
c<-b+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)
#all of this gets rid of the grey grid and legends
e<-d+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
#this sets x and y axis titles
f<-e+ labs(x = "Habitat\n", y = "\nEstimated marginal means and 95% CIs", subtitle = "C")
#this sets axis label size
g<-f+theme(axis.text.x = element_text(size = 12, colour = 'black')) +theme(axis.text.y = element_text(size = 12, colour = 'black'))
#this sets axis title size and there is your finished summary plot!
g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "Aegopodium_podagrarium" = "Aegopodium podagrarium",
  "Freshwater" = "Freshwater",
  "Grassland" = "Grassland",
  "Grassland_Freshwater" = "Freshwater (nymph)\nGrassland (adult)",
  "Vincetoxicum_hirundinaria" = "Vincetoxicum hirundinaria"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

p3 <- g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "Aegopodium_podagrarium" = "Aegopodium podagrarium",
  "Freshwater" = "Freshwater",
  "Grassland" = "Grassland",
  "Grassland_Freshwater" = "Freshwater (nymph)\nGrassland (adult)",
  "Vincetoxicum_hirundinaria" = "Vincetoxicum hirundinaria"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

### combine plots ############
# combine plots
#p4 <- plot_grid(p1, p2, p3, nrow = 3, align = "hv")
#save_plot("gut-plots/emmeans_community_composition_compiled.png", p4, base_height = 12, base_width = 8)

p4 <- plot_grid(speciesplot, dietplot, habitatplot, nrow = 3, align = "hv")

y.grob <- textGrob("Estimated marginal means of community composition and SE\n", gp = gpar(col="black", fontsize=14))

communityplots <- grid.arrange(arrangeGrob(p4, bottom = y.grob))
save_plot("gut-plots/emmeans_community_compiled_colours.png", communityplots, base_height = 12, base_width = 8)

### NOT USED IN MS: nMDS plots ############

Y$abundance <- rowSums(norm.data)
Y$labels <- metadata_res$Diet
Y$colour <- metadata_res$Species
Y$shape <- metadata_res$Diet

# plot all samples
ggplot(Y, aes(x=NMDS1, y=NMDS2, colour = colour, group = labels)) +
  geom_point() + theme_minimal() + labs(title = "NMDS plot of normalised samples gut bacteria by species and diet") + labs(colour = "Species") + stat_ellipse(aes(colour = labels), show.legend = F)

# with fill --> diet
Y$labels <- metadata_res$Diet
diet.plot <- ggplot(Y, aes(x=NMDS1, y=NMDS2, group = labels))
diet.plot <- diet.plot + geom_point(aes(colour = labels, size = abundance), shape = 16, alpha = 0.7, show.legend = T) 
#diet.plot <- diet.plot + stat_ellipse(aes(fill = labels), geom = "polygon", alpha = 0.2) + scale_fill_manual(values = plotcols[c(1,3,4)])
diet.plot <- diet.plot + labs(size ="Relative SV \nabundance", colour ="Diet") + scale_colour_manual(values = plotcols)
diet.plot

# with fill --> species
Y$labels <- metadata_res$Species
species.plot <- ggplot(Y, aes(x=NMDS1, y=NMDS2, group = labels))
species.plot <- species.plot + geom_point(aes(colour = labels, size = abundance, shape = shape), alpha = 0.7, show.legend = T) 
#species.plot <- species.plot + stat_ellipse(aes(fill = shape), geom = "polygon", alpha = 0.3) + scale_fill_manual(values = plotcols) 
species.plot <- species.plot + labs(size ="Relative SV \nabundance", colour ="Insect \nspecies", shape = "Diet") + scale_colour_manual(values = plotcols)
species.plot

#plot_grid(diet.plot, species.plot)

#+ geom_line(colour = "lightgrey", size = 1, show.legend = F)
#+ geom_text_repel(aes(label=labels), show.legend = FALSE) 

mean(colSums(norm.data))
median(colSums(norm.data))
qqnorm(log(colSums(norm.data)))

# plot samples with read abundance over e.g. mean or median 
ggplot(Y[Y$abundance > 4.65, ], aes(x=NMDS1, y=NMDS2, size=abundance, colour = colour)) +
  geom_point() + geom_text_repel(aes(label=labels), show.legend = FALSE) + 
  theme_minimal() + scale_size_continuous(range = c(3,8))

# with fill --> Species
Y$labels <- metadata_res$Species
species.plot <- ggplot(Y[Y$abundance > 4.65, ], aes(x=NMDS1, y=NMDS2, group = labels))
species.plot <- species.plot + geom_point(aes(colour = labels, size = abundance), shape = 16, alpha = 0.7, show.legend = F) 
species.plot <- species.plot + stat_ellipse(aes(fill = labels), geom = "polygon", alpha = 0.3) 
species.plot <- species.plot + labs(fill = "Species")
species.plot

# support for nmds - not optmial for stats
envfit(meta_mds~metadata_res$SpecimenID)
envfit(meta_mds~metadata_res$Diet)
envfit(meta_mds~metadata_res$Ecosystem)
envfit(meta_mds~metadata_res$Habitat)
envfit(meta_mds~metadata_res$Order)
envfit(meta_mds~metadata_res$Family)
envfit(meta_mds~metadata_res$Genus)
envfit(meta_mds~metadata_res$Species)

# # Homogeneity of dispersion test - betadiversity
beta <- betadisper(bray.data, metadata_res$Species) # ## Calculate multivariate dispersions
plot(beta)
boxplot(beta) ## Draw a boxplot of the distances to centroid for each group
beta.HSD <- TukeyHSD(beta) ## Tukey's Honest Significant Differences. A Tukey's test can be done to see if and which groups differ in relation to their variances 
plot(beta.HSD)
anova(beta) # use group dispersions to perform an ANOVA test
permutest(beta, permutations = 999) ## Permutation test for F

# So our groups present homogeneity among group dispersions (compositions vary similarly) IF betadisper results are un-significant ("Null hypothesis of no difference in dispersion between groups"; https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/betadisper) while having significantly different compositions IF significant adonis results
# Adonis test

adonis(bray.data ~ Species, data = metadata_res)
adonis(bray.data ~ Diet, data = metadata_res)
#adonis(bray.data ~ Ecosystem, data = metadata_res)
#adonis(bray.data ~ Habitat, data = metadata_res)

### NOT INCLUDED IN MS: rel abundance ##################
# mvabund package paper: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00190.x 
mvadata <- mvabund(norm.data) #abundances table. -changed from rel.data to norm.data for consistency
metadata <- metadata_res # Creates X, a table a bunch of columns with independent variables in. i.e. metadata.
#foo <- mvformula(mvadata ~ metadata$Type + metadata$Region + metadata$Ground + metadata$Species) #Makes a formula called 'foo', that is spiddat against metadata.
foo <- mvformula(mvadata ~ metadata$Ecosystem + metadata$Habitat + 
                   metadata$Diet + metadata$Species) #Makes a formula called 'foo', that is spiddat against metadata.
lm.data <- manylm(foo)
glm.data <- manyglm(foo, family= "negative.binomial")
## Obtain a summary of test statistics using residual resampling:
#summary.glm <- summary(glm.data, nBoot=500)
##Calculate a ANOVA Table:
anova.lm <- anova(lm.data, nBoot = 500)
anova.glm <- anova(glm.data, nBoot = 500) 
plot(glm.data)
coefplot(glm.data)
coefficients(glm.data)

############Richness analyses###################
richness <- rowSums(t(rarefied.data) > 0) # OTU richness for each sample
metadata <- data.frame(metadata_res, richness)
qqnorm(log(richness))

### correlation plot ###############
metadata_rich <- metadata %>% rownames_to_column(var = "PCRID")
test <- metadata_rich %>% select(PCRID, Species, Habitat, Diet, richness)
str(test)
test1 <- spread(test, Habitat, richness, fill = 0)
test2 <- spread(test, Diet, richness, fill = 0)
test3 <- merge(test1, test2, by = "PCRID")
test4 <- test3[, c(1, 4:8, 11:13)]

rownames(test4) <- NULL
test <- test4 %>% column_to_rownames(var = "PCRID")
p <- cor(test)

# add significance
res1 <- cor.mtest(test, conf.level = .95, na.action("na.omit"))
res2 <- cor.mtest(test, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = plotcols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE)

# herbivore diet and aegopodium_podagrarium habitat are highly correlated, so Palomena prasina should be excluded from the analysis - not a problem in gls models!
#metadata_rich <- metadata %>% filter(!Species == "prasina")
#metadata_rich <- metadata_rich %>% filter(!Species == "griseoaptera")
#metadata_rich <- droplevels(metadata_rich)

### chris' richness analysis ############
# From Chris' script - added log transformation of richness but I am still uncertain of the correct model - species and ecosystem can not fit in the same model due to correlation
model1 <- glm(log(richness) ~ Diet+Species, data = metadata)

model3 <- aov(log(richness) ~ metadata$Species + metadata$Diet) # examine if another linear model gives the same approx results, it does

qqnorm(resid(model1))
summary(model1)
drop1(model1, test = "Chisq")
cbind(estimate=coef(model1),confint(model1))
t <- emmeans(model1, "Species", by = "Diet", type = "response")
contrast(t)

qqnorm(resid(model3))
drop1(model3, test = "Chisq")
t <- emmeans(model3, "Species", by = "Diet", type = "response")
contrast(t)

plot(t) + theme_bw() + 
  labs(x = "Estimated marginal mean (sequence variant richness)", y = "Species")

### boxplot of observed ASV richness #####
table(metadata$Species)
box1 <-
  metadata %>% mutate(
    Species = fct_relevel(
      Species,
      "griseoaptera",
      "cyathigerum",
      "prasina",
      "cinerea",
      "equestris",
      "calcarata",
      "limbatus"
    )
  ) %>% ggplot(aes(x = Diet, y = richness, fill = Species)) + geom_boxplot(show.legend = F) + 
  theme(plot.title=element_text(hjust=0)) + ggtitle('A') + 
  labs(x = NULL, y = "Observed ASV richness") + scale_y_log10() + scale_fill_manual(values = plotcols) + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
    labels = c(
      "calcarata" = "S. calcarata",
      "cinerea" = "Ne. cinerea",
      "equestris" = "L. equestris",
      "griseoaptera" = "Ph. griseoaptera",
      "limbatus" = "Na. limbatus",
      "prasina" = "Pa. prasina"
    )
  ) 

# ---- other richness visualisation ----

# boxplot of log transformed richness for diet and coloured by habitat
ggplot(data = metadata, aes(x = Diet, y = log(richness), fill = Ecosystem)) + geom_boxplot() + 
  theme(plot.title=element_text(hjust=0)) + ggtitle('c') + 
  labs(x = NULL, y = "Richness")

# boxplot of log transformed richness by insect species and coloured by habitat 
ggplot(data = metadata, aes(x = Species, y = log(richness), fill = Ecosystem)) + geom_boxplot() + 
  theme(plot.title=element_text(hjust=0)) + ggtitle('c') + 
  labs(x = NULL, y = "Richness") + theme(axis.text.x = element_text(face="bold", 
                                                                 size=10, angle=90))
### richness: stacked bar plots ###########
# get stats for taxonomic assignment
taxa_res %>% group_by(Class) %>% tally()
apply(taxa_res, 2, function(x){ length(which(!is.na(x))) })

taxonomy <- select(taxa_res, "Order") 

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(rarefied.data, var = "otuid")
taxonomy <- rownames_to_column(taxonomy, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

data <- left_join(otus, taxonomy, by = "otuid")
data <- column_to_rownames(data, var = "otuid")

# working with reshape2 to wrangle data from wide to long format
longdata <- reshape2::melt(data) 
longdata <- longdata[longdata$value!=0,]
longdata <- aggregate(. ~ Order + variable, data = longdata, FUN = sum) # if values in the insect order column and the sample (variable) column are identical, then sum up how many unique otus there were in the sample from that insect order (richness)

# First, include metadata and rename columns to be more descriptive of their content
meta <- metadata_res %>% rownames_to_column(var = "variable") %>% select(variable, SpecimenID, Family, Species, Diet)

str(longdata)
plot_data <- merge(meta, longdata, by.x = "variable", by.y = "variable")
plot_data <- as_tibble(plot_data)
plot_data <- plot_data %>% dplyr::rename(sampleid = variable, richness = value)

richness_test <- plot_data %>% dplyr::group_by(sampleid, Species, Diet) %>% dplyr::summarise(richness = sum(richness)) # since data is normalised it appear in relative abundance and does not make sense to use

mean(plot_data$richness)
min(plot_data$richness)
quantile(plot_data$richness)

# Load RColorBrewer
library(RColorBrewer)
# Classic palette BuPu, with 4 colors
coul <- brewer.pal(12, "Paired") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(39)

# Plot it
pie(rep(1, length(coul)), col = coul , main="") 

# Column/bar chart visualization
plot_data %>%  filter(richness > 9) %>% ggplot(aes(x = Species, y = richness)) + geom_bar(stat = "identity", aes(fill = Order), na.rm = T, position = "fill") + labs(x = "Insect species", y = "", title = "ASV richness by bacteria order") + theme_bw() + facet_grid(. ~ Diet, scales = "free") + labs(fill = "Bacteria order") + theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5), axis.text.y = element_text(face = "bold", size = 9), plot.title = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))+ guides(fill=guide_legend(ncol=2)) + scale_fill_manual(values = coul)

# Stacked bar chart oof OTU richness
g <- ggplot(plot_data, aes(Species, richness))
g + geom_bar(stat = "identity", aes(fill = Order), show.legend = F) + labs(x = "Insect species", y = "ASV richness", fill = "Bacteria order") + theme_classic2() + guides(fill=guide_legend(ncol=2)) #+ scale_fill_manual(values = custom_11) # position = "fill"in geom_bar gives relative

# that really didn't show us anaything usefull. Will try with ampvis
library(ampvis2)
taxonomy <- rownames_to_column(taxa_res, var = "otuid")
#taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package
rarefied.otus <- rownames_to_column(rarefied.data, var = "otuid") 
otutable <- left_join(otus, taxonomy, by = "otuid")# need to use unrarefied data and re-rarefy, otherwise it will be rarefied twice

otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames
# metadata must have sample data in the first column and column classes matter 
str(metadata_res)
metadata <- rownames_to_column(metadata_res, var = "otuid")

metadata <-
  metadata %>% dplyr::rename(InsectOrder = Order, InsectFamily = Family, InsectGenus = Genus, InsectSpecies = Species, InsectSuborder = Suborder, InsectInfraorder = Infraorder)

amp <- amp_load(otutable = otutable, metadata = metadata)

# heatmap vizulaization ----
# convert read abundances into percentages
#psn <- transform_sample_counts(ps, function(x) x/sum(x) * 100)

# heatmap (normalise gives the proportion of the sample)
amp_heatmap(
  data = amp,
  group_by = "InsectSpecies",
  normalise = T,
  tax_aggregate = "Genus",
  tax_add = "Family",
  plot_values = F,
  plot_colorscale = "log10",
  color_vector = c("#E1EAEC", "#989EA3", "#5C6063"),
  tax_empty = "remove",
  tax_show = 50,
)+ guides(fill = guide_legend(title = "Relative ASV abundance", reverse = T)) 

# venn diagram: https://madsalbertsen.github.io/ampvis2/reference/amp_venn.html 
amp_venn(amp,
         group_by = "InsectOrder",
         cut_a = 0.1,
         cut_f = 10)
# boxplot/point plot
amp_boxplot(
  amp,
  group_by = "InsectSpecies",
  sort_by = "mean",
  plot_type = "boxplot",
  tax_aggregate = "Family",
  tax_show = 25,
  tax_empty = "remove",
  plot_log = TRUE
) + ylab('Relative abundance') + guides(colour = guide_legend(title = "?", reverse = T)) + scale_color_manual(values = plotcols)

# ordination: https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html 
ord <- amp_ordinate(
  amp,
  type = "NMDS",
  distmeasure = "bray",
  transform = "none",
  sample_color_by = "InsectSpecies",
  sample_colorframe = T
)

ord + scale_color_manual(values = plotcols) + scale_fill_manual(values = plotcols) + guides(fill=guide_legend("Insect species"), color = FALSE)

### amp richness indices ############
# richness indices
amp_alphadiv <-
  amp_alphadiv(
    amp,
    measure = c("observed", "shannon", "simpson", "invsimpson"),
    richness = TRUE,
    rarefy = 1830
  )

# we can combine this with the other observed ASV richness in metadata
diversity_indices <- amp_alphadiv %>% select(otuid, Reads, ObservedOTUs, Shannon, Simpson, invSimpson, ACE, Chao1)

metadata <- metadata %>% rownames_to_column(var = "otuid")

richness_data <- merge(metadata, diversity_indices, by = "otuid")
metadata <- metadata %>% column_to_rownames(var = "otuid")

### boxplot of richness indices ############
# chao1 (richness estimate)
box2 <-
  richness_data %>% mutate(
    InsectSpecies = fct_relevel(
      InsectSpecies,
      "griseoaptera",
      "cyathigerum",
      "prasina",
      "cinerea",
      "equestris",
      "calcarata",
      "limbatus"
    )
  ) %>% ggplot(aes(x = Diet, y = Chao1, fill = InsectSpecies)) + geom_boxplot(show.legend = F) + 
  theme(plot.title=element_text(hjust=0)) + ggtitle('B') + 
  labs(x = NULL, y = "Chao1 richness estimate") + scale_y_log10() + scale_fill_manual(values = plotcols) + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
    labels = c(
      "calcarata" = "S. calcarata",
      "cinerea" = "Ne. cinerea",
      "equestris" = "L. equestris",
      "griseoaptera" = "Ph. griseoaptera",
      "limbatus" = "Na. limbatus",
      "prasina" = "Pa. prasina"
    )
  ) 

# Shannon (diversity index)
box3 <-
  richness_data %>% mutate(
    InsectSpecies = fct_relevel(
      InsectSpecies,
      "griseoaptera",
      "cyathigerum",
      "prasina",
      "cinerea",
      "equestris",
      "calcarata",
      "limbatus"
    )
  ) %>% ggplot(aes(x = Diet, y = Shannon, fill = InsectSpecies)) + geom_boxplot(show.legend = T) + 
  theme(plot.title=element_text(hjust=0)) + ggtitle('C') + 
  labs(x = NULL, y = "Shannon diversity index", fill = "Insect species") + scale_y_log10() + scale_fill_manual(values = plotcols) + theme(plot.subtitle = element_text(size = 20, face = "bold"), legend.position = "none", legend.title = element_text()) + scale_x_discrete(
    labels = c(
      "calcarata" = "S. calcarata",
      "cinerea" = "Ne. cinerea",
      "equestris" = "L. equestris",
      "griseoaptera" = "Ph. griseoaptera",
      "limbatus" = "Na. limbatus",
      "prasina" = "Pa. prasina"
    )
  )

# extract a legend that is laid out horizontally
legend_b <- get_legend(
  box3 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# arrange the three plots in a single row
boxplots <- plot_grid(
  box1 + theme(legend.position="none"),
  box2 + theme(legend.position="none"),
  box3 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 3
)
boxplots

#boxplots <- plot_grid(box1, box2, box3, nrow = 3, align = "hv")
# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
bplot <- plot_grid(boxplots, legend_b, ncol = 1, rel_heights = c(1, .1))

save_plot("gut-plots/amp_richness_boxplots_compiled.png", bplot, base_height = 12, base_width = 8)

# rarefaction curve on read abundance
# join otus and taxonomy
otus <- rownames_to_column(data_res, var = "otuid")
#taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy <- taxonomy %>% select(otuid, Kingdom, Phylum, Class, Order, Family, Genus, Species)
otutable <- left_join(otus, taxonomy, by = "otuid")
# revert otuids back to rownames if needed
otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(metadata)
metadata <- rownames_to_column(metadata, var = "otuid")

# load amp data
#otutable <- otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = metadata)
# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp,
    facet_by = "Diet",
    stepsize = 100,
    color_by = "InsectSpecies",
    facet_scales = "free"
  )

min(colSums(otus[, 2:64]))
max(colSums(otus[, 2:64]))

rarecurve + geom_line(size = 1.5, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 120000),
  breaks = seq(from = 0, to = 120000, by = 5000)
) + theme(legend.position = "bottom") + guides(colour=guide_legend(ncol=8)) + ylab("Number of observed SVs") + scale_colour_manual(values = plotcols)

# plot richness indices on log10 transformed y-axis - based on rarefied richness
amp_alphadiv %>% tidyr::gather("id", "Richness", 15:20) %>% ggplot(., aes(InsectSpecies, Richness, fill = InsectSpecies)) + geom_boxplot(show.legend = T) + facet_wrap(~ id, ncol = 3) + scale_y_log10() + theme_pubclean() + labs(fill = "Insect species") + ylab("Richness (log-transformed)") + theme(axis.text.x=element_blank(), legend.position = "bottom") + scale_fill_manual(values = plotcols)

# test if diversity means are different
anova <- aov(Chao1 ~ Diet + InsectSpecies, data = amp_alphadiv)
summary(anova) 
plot(anova, 1) # visualise residuals
plot(anova, 2) # investigate normality

# post hoc analysis
library(emmeans)
emmeans(anova, specs = pairwise ~InsectSpecies)

### gls on richness observed, estimated and diversity index #########
# On insect species
model1 <- glm(log(ObservedOTUs) ~ InsectSpecies, data = amp_alphadiv)
model2 <- glm(log(Chao1) ~ InsectSpecies, data = amp_alphadiv)
model3 <- glm(log(Shannon) ~ InsectSpecies, data = amp_alphadiv)

# observed richness
qqnorm(resid(model1))
summary(model1)
drop1(model1, test = "Chisq")
cbind(estimate=coef(model1),confint(model1))
t <- emmeans(model1, "InsectSpecies", type = "response")
contrast(t)
r.squaredGLMM(model1)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means

# estimated richness thats takes the importance of rare ASVs into account (Chao1)
qqnorm(resid(model2))
summary(model2)
drop1(model2, test = "Chisq")
cbind(estimate=coef(model2),confint(model2))
t <- emmeans(model2, "InsectSpecies", type = "response")
contrast(t)
r.squaredGLMM(model2)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means
# richness index
qqnorm(resid(model3))
summary(model3)
drop1(model3, test = "Chisq")
cbind(estimate=coef(model3),confint(model3))
t <- emmeans(model3, "InsectSpecies", type = "response")
contrast(t)
r.squaredGLMM(model3)
plot(t)

ems_means <- summary(t) %>%
  as.data.frame()

p <-
  ems_means %>% mutate(
    InsectSpecies = fct_relevel(
      InsectSpecies,
      "griseoaptera",
      "cyathigerum",
      "prasina",
      "cinerea",
      "equestris",
      "calcarata",
      "limbatus"
    )
  ) %>% ggplot(aes(InsectSpecies, response))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)

richplot_1 <- c + geom_pointrange(aes(ymin = response-SE, ymax = response+SE, colour = InsectSpecies),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "Species", y = "", subtitle = "A") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                        labels = c(
                          "griseoaptera" = "Ph. griseoaptera",
                          "cyathigerum" = "E. cyathigerum",
                          "prasina" = "Pa. prasina",
                          "cinerea" = "Ne. cinerea",
                          "equestris" = "L. equestris",
                          "calcarata" = "S. calcarata",
                          "limbatus" = "Na. limbatus"
                        )
                      ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols)

### Richness insect diet #################
model1 <- glm(log(ObservedOTUs) ~ Diet, data = amp_alphadiv)
model2 <- glm(log(Chao1) ~ Diet, data = amp_alphadiv)
model3 <- glm(log(Shannon) ~ Diet, data = amp_alphadiv)

# observed richness
qqnorm(resid(model1))
summary(model1)
drop1(model1, test = "Chisq")
cbind(estimate=coef(model1),confint(model1))
t <- emmeans(model1, "Diet", type = "response")
contrast(t)
r.squaredGLMM(model1)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means

# estimated richness thats takes the importance of rare ASVs into account (Chao1)
qqnorm(resid(model2))
summary(model2)
drop1(model2, test = "Chisq")
cbind(estimate=coef(model2),confint(model2))
t <- emmeans(model2, "Diet", type = "response")
contrast(t)
r.squaredGLMM(model2)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means

# richness index
qqnorm(resid(model3))
summary(model3)
drop1(model3, test = "Chisq")
cbind(estimate=coef(model3),confint(model3))
t <- emmeans(model3, "Diet", type = "response")
contrast(t)
r.squaredGLMM(model3)
plot(t, comparisons = T)

plot(t)

ems_means <- summary(t) %>%
  as.data.frame()

p <-
  ems_means %>% mutate(
    Diet = fct_relevel(
      Diet,
      "herbivore",
      "omnivore",
      "predator"
    )
  ) %>% ggplot(aes(Diet, response))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)

richnplot_2 <- c + geom_pointrange(aes(ymin = response-SE, ymax = response+SE, colour = Diet),
                    size = 1.5) + theme_minimal_grid() + theme(
                      legend.title = element_blank(),
                      legend.key = element_rect(size = 0.1),
                      legend.key.size = unit(1, 'cm')
                    ) + labs(x = "Diet", y = "", subtitle = "B") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                      labels = c(
                        "herbivore" = "Herbivore",
                        "omnivore" = "Omnivore",
                        "predator"="Predator"
                      )
                    ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols)

### habitat richness ####################
model1 <- glm(log(ObservedOTUs) ~ Habitat, data = amp_alphadiv)
model2 <- glm(log(Chao1) ~ Habitat, data = amp_alphadiv)
model3 <- glm(log(Shannon) ~ Habitat, data = amp_alphadiv)

# observed richness
qqnorm(resid(model1))
summary(model1)
drop1(model1, test = "Chisq")
cbind(estimate=coef(model1),confint(model1))
t <- emmeans(model1, "Habitat", type = "response")
contrast(t)
r.squaredGLMM(model1)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means

# estimated richness thats takes the importance of rare ASVs into account (Chao1)
qqnorm(resid(model2))
summary(model2)
drop1(model2, test = "Chisq")
cbind(estimate=coef(model2),confint(model2))
t <- emmeans(model2, "Habitat", type = "response")
contrast(t)
r.squaredGLMM(model2)
ems_means <- summary(t) %>%
  as.data.frame()
ems_means

# richness index
qqnorm(resid(model3))
summary(model3)
drop1(model3, test = "Chisq")
cbind(estimate=coef(model3),confint(model3))
t <- emmeans(model3, "Habitat", type = "response")
contrast(t)
r.squaredGLMM(model3)
plot(t, comparisons = T)

plot(t)

ems_means <- summary(t) %>%
  as.data.frame()
ems_means

p <-
  ems_means %>% mutate(
    Habitat = fct_relevel(
      Habitat,
      "Vincetoxicum_hirundinaria",
      "Grassland_Freshwater",
      "Grassland",
      "Freshwater",
      "Aegopodium_podagrarium"
    )
  ) %>% ggplot(aes(Habitat, response))

c<-p+coord_flip()+scale_size_area(max_size = 1.5)

richplot_3 <- c + geom_pointrange(aes(ymin = response-SE, ymax = response+SE, colour = Habitat),
                                   size = 1.5) + theme_minimal_grid() + theme(
                                     legend.title = element_blank(),
                                     legend.key = element_rect(size = 0.1),
                                     legend.key.size = unit(1, 'cm')
                                   ) + labs(x = "Habitat", y = "", subtitle = "C") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_x_discrete(
                                     labels = c(
                                       "Aegopodium_podagrarium" = "Aegopodium podagrarium",
                                       "Grassland" = "Grassland",
                                       "Grassland_Freshwater" = "Freshwater (nymph)\nGrassland (adult)",
                                       "Freshwater" = "Freshwater",
                                       "Vincetoxicum_hirundinaria" = "Vincetoxicum hirundinaria"
                                     )
                                   ) + theme(axis.text.x = element_text(size = 12, angle = 90), legend.position = "none") + scale_colour_manual(values = plotcols)

### combine plots ############
# combine plots
#p4 <- plot_grid(p1, p2, p3, nrow = 3, align = "hv")
#save_plot("gut-plots/emmeans_community_composition_compiled.png", p4, base_height = 12, base_width = 8)

p4 <- plot_grid(richplot_1, richnplot_2, richplot_3, nrow = 3, align = "hv")

y.grob <- textGrob("Estimated marginal means of Shannon Diversity Index and SE\n", gp = gpar(col="black", fontsize=14))

richnessplots <- grid.arrange(arrangeGrob(p4, bottom = y.grob))
save_plot("gut-plots/emmeans_richness_compiled_colurs_side.png", richnessplots, base_height = 12, base_width = 8)
