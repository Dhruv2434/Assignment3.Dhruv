# Load necessary libraries
library(rentrez)         #for downloading data from NCBI
library(tidyverse)       #for cleaning up the data 
library(Biostrings)      # For DNA sequence manipulation with DNAStringSet
library(ape)             # For phylogenetic analysis and tree visualization
library(phangorn)        # Additional phylogenetic analysis tools
library(dplyr)           # For data manipulation and filtering
library(muscle)          # For sequence alignment with MUSCLE
library(ape)            #for bootstrapping analysis of trees 
library(phytools)       #for converting trees to dendograms 
library(dendextend)     #form comparing dendograms 


#### - PART 1 - Original data files ----

# Defining file paths for COI and HLA genes
#setwd( "../ST assignment 2")
coi_path <- ("../ST assignment 2/COI")
hla_path <- ("../ST assignment 2/HLA")

# List of species for COI and HLA files

species_names <- c("Cat", "Chimpanzee", "GiantPanda", "Gibbon", "Gorilla", 
                   "GrayWolf", "GraySquirrel", "Human", "Mouse", "Rat", 
                   "RedFox", "Rhesus", "Tiger")

# Function to read all FASTA files in a given directory and create a data frame

read_fasta_files <- function(path, species_list, suffix) {
  sequences <- list()
  for (species in species_list) {
    file_path <- file.path(path, paste0(species, suffix, ".fasta"))
    seq_data <- readDNAStringSet(file_path)
    sequences[[species]] <- data.frame(
      processid = species, 
      nucleotides2 = as.character(seq_data),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, sequences)
}

# Read and combine COI and HLA data into data frames

dfCOI <- read_fasta_files(coi_path, species_names, "COI")

dfHLA <- read_fasta_files(hla_path, species_names, "HLA")


# Filter and clean data to retain only the sequences of interest

# Here we make sure each dataframe has unique sequences and processids for consistency

dfCOI <- dfCOI %>% filter(processid %in% species_names)

dfHLA <- dfHLA %>% filter(processid %in% species_names)

# I would add something in here to visualize the sequence lengths to make sure that they are all similar and that you don't have any major issues with the data! 

dfCOI <- dfCOI %>%
  mutate(seqlength = nchar(dfCOI$nucleotides2))

dfHLA <- dfHLA %>%
  mutate(seqlength = nchar(dfHLA$nucleotides2))

# visualizing the sequence lengths with a histogram

hist(x = dfCOI$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #most are between 0 and 5000 but this is one which is between 15000 and 20000, this may pose an issue- see part 2 for downloading data directly from NCBI 
hist(x = dfHLA$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #this one looks good! most of the sequences are between 0 and 6000 base pairs 

#visualizing the sequence length for each species with a bar plot

ggplot(dfCOI, aes(x= processid, y= seqlength))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Sequence length")+
  xlab("Species")

#We can see here that it is a Rhesus monkey which has the much longer sequence length, it may be useful to remove this species or see part 2 for downloading data from NCBI 

ggplot(dfHLA, aes(x= processid, y= seqlength))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Sequence length")+
  xlab("Species") 


# Converting sequences to DNAStringSet format and name them by processid
dfCOI$nucleotides2 <- DNAStringSet(dfCOI$nucleotides2)
names(dfCOI$nucleotides2) <- dfCOI$processid

dfHLA$nucleotides2 <- DNAStringSet(dfHLA$nucleotides2)
names(dfHLA$nucleotides2) <- dfHLA$processid

# Aligning COI sequences using MUSCLE
cat("Aligning COI sequences...\n")
dfCOI.alignment <- DNAStringSet(muscle::muscle(dfCOI$nucleotides2))

# Aligning HLA sequences using MUSCLE
cat("Aligning HLA sequences...\n")
dfHLA.alignment <- DNAStringSet(muscle::muscle(dfHLA$nucleotides2))


# Inspect the alignments
print(dfCOI.alignment)  # Shows aligned COI sequences
print(dfHLA.alignment)  # Shows aligned HLA sequences


# Building Phylogenetic Trees for COI and HLA
# Convert aligned sequences to phyDat format for tree-building

coi_phyDat <- phyDat(as.matrix(dfCOI.alignment), type = "DNA")

hla_phyDat <- phyDat(as.matrix(dfHLA.alignment), type = "DNA")


# Neighbor-Joining Tree for COI

cat("Building Neighbor-Joining Tree for COI...\n")
coi_dist <- dist.ml(coi_phyDat)           # Calculate distance matrix
coi_nj_tree <- NJ(coi_dist)                # Neighbor-Joining tree
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# I would preform a model test to determine the best model for your data this way you can determine which model works best! 

modelTest_coi <- modelTest(coi_phyDat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
View(modelTest_coi)

# Choosing the best model based on AIC since the sample size is small this is the better model to follow 
best_model_coi <- modelTest_coi[which.min(modelTest_coi$AIC), "Model"]
print(paste("Best modelAIC:", best_model_coi)) #GTR is the best model so I changed the next section for this model  


# Maximum Likelihood Tree for COI

cat("Building Maximum Likelihood Tree for COI...\n")
coi_ml_tree <- pml(coi_nj_tree, coi_phyDat)
coi_ml_tree <- optim.pml(coi_ml_tree, model = "HKY")
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Neighbor-Joining Tree for HLA

cat("Building Neighbor-Joining Tree for HLA...\n")
hla_dist <- dist.ml(hla_phyDat)           # Calculate distance matrix
hla_nj_tree <- NJ(hla_dist)                # Neighbor-Joining tree
plot(hla_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

#again do a model test to find the best model 

modelTest_hla <- modelTest(hla_phyDat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
View(modelTest_hla)

# Choosing the best model based on AIC since the sample size is small this is the better model to follow 
best_model_hla <- modelTest_hla[which.min(modelTest_coi$AIC), "Model"]
print(paste("Best modelAIC:", best_model_hla)) #GTR is the best model so I changed the next section for this model  

# Maximum Likelihood Tree for HLA

cat("Building Maximum Likelihood Tree for HLA...\n")
hla_ml_tree <- pml(hla_nj_tree, hla_phyDat)
hla_ml_tree <- optim.pml(hla_ml_tree, model = "HKY")
plot(hla_ml_tree$tree, main = "Maximum Likelihood Tree for HLA Sequences")

#Comparing Tree Topologies for COI and HLA
# Here, we calculate a distance metric between the COI and HLA trees

cat("Comparing COI and HLA tree topologies...\n")
tree_distance <- RF.dist(coi_nj_tree, hla_nj_tree)
cat("Robinson-Foulds distance between COI and HLA trees: ", tree_distance, "\n")

# I would add in boot strapping here so that you can determine which tree best represents the data beyond determining that the trees are different 

#calcutaing bootstrapping for COI

bs_coi <- bootstrap.pml(coi_ml_tree, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

#calculating bootstrapping for HLA 

bs_hla <- bootstrap.pml(hla_ml_tree, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

#graphing trees with bootstrapping values for both 

plotBS(midpoint(coi_ml_tree$tree), bs_coi, type = "phylogram", cex = 0.6)

plotBS(midpoint(hla_ml_tree$tree), bs_hla, type = "phylogram", cex = 0.6)

#bootstrapping values are much higher for the COI phylogeny which indicates that it is a much better representation of the phylogeny and that HLA likely is too polymorphic to represent genetic relationships 

# Visualizing Trees Side-by-Side for Comparison
par(mfrow = c(2, 2))  # Set up 2x2 plotting area

# Plot 1: Neighbor-Joining Tree for COI
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# Plot 2: Maximum Likelihood Tree for COI
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Plot 3: Neighbor-Joining Tree for HLA
plot(hla_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

# Plot 4: Maximum Likelihood Tree for HLA
plot(hla_ml_tree$tree, main = "Maximum Likelihood Tree for HLA Sequences")

# Save alignments and tree plots for future reference and analysis
writeXStringSet(dfCOI.alignment, file = "COI_aligned_sequences.fasta")
writeXStringSet(dfHLA.alignment, file = "HLA_aligned_sequences.fasta")
saveRDS(coi_nj_tree, file = "COI_NJ_Tree.rds")
saveRDS(coi_ml_tree, file = "COI_ML_Tree.rds")
saveRDS(hla_nj_tree, file = "HLA_NJ_Tree.rds")
saveRDS(hla_ml_tree, file = "HLA_ML_Tree.rds")

cat("Analysis complete. Results saved.\n")

#Instead of just plotting the trees next to each other it would be good to make a visual comparison of the trees 
#using dendextend to make a visual comparison of the two different genes for the two maximum likelihood trees

par(mfrow = c(1, 1))  # Set up 1x1 plotting area

# Converting coi tree from phylo to dendrogram
rooted_coi_tree <- midpoint.root(coi_ml_tree$tree) #rooting the tree 
coi_ultrametric <- chronos(rooted_coi_tree) #time-calibrate the tree
coi_dend <- as.dendrogram(coi_ultrametric) #turn the tree into a dendrogram 
plot(coi_dend, main = "COI Dendrogram")

#repeating for hla 
rooted_hla_tree <- midpoint.root(hla_ml_tree$tree) 
hla_ultrametric <- chronos(rooted_hla_tree, model = "relaxed") #allows for different rates of evolution along the tree which is needed for a gene as polymorphic as HLA
hla_dend <- as.dendrogram(hla_ultrametric)
plot(hla_dend, main = "HLA Dendrogram")

dl <- dendlist(coi_dend, hla_dend)
dl <- dendlist(highlight_branches_col(coi_dend, 
                                      values = c("#000000", "#0072B2","#CC79A7","#009E73","#E69F00")), 
               highlight_branches_col(hla_dend, values = c("#000000", "#0072B2","#CC79A7","#009E73","#E69F00"))) 
#here we colour the brnaches of the tree in the same format for both dendograms so that the branching of each tree can be compared 

tanglegram(dl, lab.cex = 1.3, margin_inner = 8, margin_outer = 2, common_subtrees_color_lines = TRUE,
           common_subtrees_color_lines_default_single_leaf_color = "grey", main_left = "COI", main_right = "HLA")

#here we compare the two trees with lines matching the species on each tree. In most cases these lines would be coloured when species are in the same location on the tree, this isn't shown here because none of the species are in the same location (or in the same group/clade)


#### - PART 2 - downloading data from NCBI ----

#I would download this data directly from NCBI and then filter for the species you want therefore you get a larger breadth of the sequences for the species and you can calculate the centroid sequence for each species before analysis

#loading source function to download larger data sets from NCBI

source("Entrez_Functions.R")

#create a search object for COI gene 

Mammalscoi_search <- entrez_search(db = "nuccore", term = "(Mammals[ORGN] AND COI[GENE] AND 200:1500[SLEN])", use_history = T)

Mammalscoi_search # 46879 sequences were found
Mammalscoi_search$web_history # A web history object was successfully made

#fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Mammals[ORGN] AND COI AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Mammals_COI")
#I uploaded the data files as well so if you end up running this I wouldnt run this step because it takes a while and its a lot of data 

#combing the many 100 Mammal_COI sequence files into one data frame and moving it into R

Mammals_COI <- MergeFastaFiles(filePattern = "Mammals_COI*")

#creating a search for the RAG1 gene I changed it from the HLA gene because from my research this gene is more specific to humans and there wasn't a sequence for each species you had from NCBI 

Mammalshla_search <- entrez_search(db = "nuccore", term = "(Mammals[ORGN] AND HLA AND 200:1500[SLEN])", use_history = T)

Mammalshla_search # 46879 sequences were found
Mammalshla_search$web_history # A web history object was successfully made

# fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Mammals[ORGN] AND HLA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "MammalsHLA")

# combing the many 100 RAG1 sequence files into one data frame and moving it into R

Mammals_HLA <- MergeFastaFiles(filePattern = "MammalsHLA*")

#create a list showing which Latin name corresponds to which common name ** I realized that gibbon is a family whist most others you have here are species so I chose a random species of Gibbon but you can obviously change which species you want to use 

latin_to_common <- c(
  "Felis catus" = "Domestic Cat",
  "Pan troglodytes" = "Chimpanzee",
  "Ailuropoda melanoleuca" = "Giant Panda",
  "Hylobates lar" = "White-handed Gibbon",
  "Gorilla gorilla" = "Gorilla",
  "Canis lupus" = "Gray Wolf",
  "Sciurus carolinensis" = "Gray Squirrel",
  "Homo sapiens" = "Human",
  "Mus musculus" = "Mouse",
  "Rattus norvegicus" = "Rat",
  "Vulpes vulpes" = "Red Fox",
  "Macaca mulatta" = "Rhesus",
  "Panthera tigris" = "Tiger")

#creating a function to clean up the NCBI file to only have the species of choice and add in the common name

missing.data <- 0.01

cleaning_sequencefile <- function(data) {
  result <- data %>%
    mutate(Title = gsub("\\bPREDICTED:\\s*", "", Title))  %>%  #remove predicted from before species so that these species aren't filtered out later on
    separate(Title,
             into = c("Code", "Genus", "Species", "rest"), 
             sep = " ",
             extra = "merge",
             fill = "right",
             convert = FALSE) %>% #separate the title column into separate groups 
    mutate(Species = paste(Genus, Species)) %>% #recombine genus and species to get full species name 
    mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>% #make sure all species names are formatted the same by finding the number of spaces
    select(Code, Species, Sequence, spaces) %>% #select for columns needed in later analysis 
    filter(spaces == 1, !is.na(spaces), !is.na(Species)) %>% #filter out missing data and for species with only one space
    mutate(Sequence2 = str_remove_all(Sequence, "^N+|N+$|-")) %>% #remove Ns at start and end of sequence and gaps
    filter(str_count(Sequence2, "N") <= (missing.data * str_count(Sequence))) %>% #remove sequences with a high proportion of Ns (0.01)
    filter(Species %in% c("Felis catus", "Pan troglodytes", "Ailuropoda melanoleuca", 
                          "Hylobates lar", "Gorilla gorilla", "Canis lupus", 
                          "Sciurus carolinensis", "Homo sapiens", "Mus musculus", 
                          "Rattus norvegicus", "Vulpes vulpes", "Macaca mulatta", "Panthera tigris"))%>% #filter for the species used in later analysis 
    mutate(Common_Name = recode(Species, !!!latin_to_common)) #add in a column for common names 
  
  return(result) 
}

#applying the function to both markers 

dfCOI_NCBI <- cleaning_sequencefile(Mammals_COI)

dfHLA_NCBI <- cleaning_sequencefile(Mammals_HLA)

#checking to make sure that all of the species you need are present 

length(unique(dfCOI_NCBI$Species)) #13 unique species 
list(unique(dfCOI_NCBI$Common_Name)) 

length(unique(dfHLA_NCBI$Species)) #13 unique species 
list(unique(dfHLA_NCBI$Common_Name)) 

#Creating a function to look at the number of sequences for each species and the average length of those sequences 

sequence_counts <- function(data) {
  result <- data %>%
    mutate(seqlength = nchar(data$Sequence)) %>% 
    group_by(Species) %>%
    summarize(sequence_count = n(),  avg_seqlength = mean(seqlength, na.rm = TRUE))
  
  return(result)
}

#applying function 

sequence_counts_COI_NCBI <- sequence_counts(dfCOI_NCBI)
View(sequence_counts_COI_NCBI)

sequence_counts_HLA_NCBI <- sequence_counts(dfHLA_NCBI)
View(sequence_counts_HLA_NCBI) #there are over 45000 sequences for humans and only around 100 for each other species 

#creating another column in the orginal file for sequence length 

dfCOI_NCBI <- dfCOI_NCBI %>%
  mutate(seqlength = nchar(dfCOI_NCBI$Sequence2))

dfHLA_NCBI <- dfHLA_NCBI %>%
  mutate(seqlength = nchar(dfHLA_NCBI$Sequence2))

# visualizing the sequence lengths with a histogram

hist(x = dfCOI_NCBI$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #most are between 0 and 5000 but this is one which is between 15000 and 20000, this may pose an issue- see part 2 for downloading data directly from NCBI 
hist(x = dfHLA_NCBI$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #this one looks good! most of the sequences are between 0 and 6000 base pairs 

#visualizing the sequence length for each species with a bar plot

ggplot(sequence_counts_COI_NCBI, aes(x= Species, y= avg_seqlength))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Sequence length")+
  xlab("Species")

#We can see here that it is a Rhesus monkey which has the much longer sequence length, it may be useful to remove this species or see part 2 for downloading data from NCBI 

ggplot(sequence_counts_HLA_NCBI, aes(x= Species, y= avg_seqlength))+
  geom_bar(stat= 'identity', fill=c("#0072B2"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title = element_text(margin = margin(t = 0, r = 10, b = 10, l = 0), size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))+
  ylab("Sequence length")+
  xlab("Species") 

#reduce the number of sequences for humans by

#randomly sampling for 200 sequences 

dfHLA_sub <- dfHLA_NCBI %>%
  filter(Species == "Homo sapiens") %>%
  slice_sample(n = 200) 

#recombing this with the other data

dfHLA_rest <- dfHLA_NCBI %>%
  filter(Species != "Homo sapiens") 

dfHLA_NCBI <- bind_rows(dfHLA_sub, dfHLA_rest)

#looking at the new number of sequences per species 

sequence_counts_HLA_sub <- sequence_counts(dfHLA_NCBI)
View(sequence_counts_HLA_sub)


