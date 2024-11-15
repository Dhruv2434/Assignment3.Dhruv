# Load necessary libraries
library(rentrez)         #for downloading data from NCBI
library(tidyverse)       #for cleaning up the data 
library(Biostrings)      # For DNA sequence manipulation with DNAStringSet
library(ape)             # For phylogenetic analysis and tree visualization
library(phangorn)        # Additional phylogenetic analysis tools
library(dplyr)           # For data manipulation and filtering
library(muscle)          # For sequence alignment with MUSCLE
library(ape)            #for bootstrapping analysis of trees 



#### - PART 1 - downloading data from NCBI ----

#I would maybe download this data directly from NCBI or BOLD and then filter for the species you want therefore you get a larger breadth of the sequences for the species and you can calculate the centroid sequence for each species before analysis

#loading source function to download larger data sets from NCBI

source("Entrez_Functions.R")

#create a search object for COI gene 

Mammalscoi_search <- entrez_search(db = "nuccore", term = "(Mammals[ORGN] AND COI[GENE] AND 200:1500[SLEN])", use_history = T)

Mammalscoi_search # 46879 sequences were found
Mammalscoi_search$web_history # A web history object was successfully made

#fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Mammals[ORGN] AND COI AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Mammals_COI")

#combing the many 100 Mammal_COI sequence files into one data frame and moving it into R

Mammals_COI <- MergeFastaFiles(filePattern = "Mammals_COI*")

#creating a search for the RAG1 gene I changed it from the HLA gene because from my research this gene is more specific to humans and there wasn't a sequence for each species you had from NCBI 

Mammalshla_search <- entrez_search(db = "nuccore", term = "(Mammals[ORGN] AND HLA AND 200:1500[SLEN])", use_history = T)

Mammalshla_search # 46879 sequences were found
Mammalshla_search$web_history # A web history object was successfully made

# fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Mammals[ORGN] AND HLA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "MammalsHLA")

# combing the many 100 RAG1 sequence files into one data frame and moving it into R

Mammals_HLAC <- MergeFastaFiles(filePattern = "MammalsHLA*")

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

missing.data <- 0.001

cleaning_sequencefile <- function(data) {
  result <- data %>%
    mutate(Title = gsub("\\bPREDICTED:\\s*", "", Title))  %>%
  separate(Title,
             into = c("Code", "Genus", "Species", "rest"),
             sep = " ",
             extra = "merge",
             fill = "right",
             convert = FALSE) %>%
    mutate(Species = paste(Genus, Species)) %>%
    mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>%
    select(Code, Species, Sequence, spaces) %>%
    filter(spaces == 1, !is.na(spaces), !is.na(Species)) %>%
    mutate(Sequence2 = str_remove_all(Sequence, "^N+|N+$|-")) %>%
    filter(str_count(Sequence2, "N") <= (missing.data * str_count(Sequence))) %>%
    filter(Species %in% c("Felis catus", "Pan troglodytes", "Ailuropoda melanoleuca", 
                          "Hylobates lar", "Gorilla gorilla", "Canis lupus", 
                          "Sciurus carolinensis", "Homo sapiens", "Mus musculus", 
                          "Rattus norvegicus", "Vulpes vulpes", "Macaca mulatta", "Panthera tigris"))%>%
    mutate(Common_Name = recode(Species, !!!latin_to_common))
  
  return(result) 
}


#applying the function to both markers 

dfCOI <- cleaning_sequencefile(Mammals_COI)

dfHLA <- cleaning_sequencefile(Mammals_HLA)

#checking to make sure that all of the species you need are present 

length(unique(dfCOI$Species)) # 13 unique species 
list(unique(dfCOI$Common_Name)) 

length(unique(dfHLA$Species)) #13 species are in the data set 
list(unique(dfHLA$Common_Name)) 

#Creating a function to look at the number of sequences for each species

sequence_counts <- function(data) {
  result <- data %>%
  mutate(seqlength = nchar(data$Sequence)) %>% 
  group_by(Species) %>%
  summarize(sequence_count = n(),  avg_seqlength = mean(seqlength, na.rm = TRUE))
  
return(result)
}

#applyting function 

sequence_counts_HLA <- sequence_counts(dfHLA)
View(sequence_counts_HLA) #there are over 45000 sequences for humans and only around 100 for each other species 

#reduce the number of sequences for humans by

#randomly sampling for 200 sequences 

dfHLA_sub <- dfHLA %>%
  filter(Species == "Homo sapiens") %>%
  slice_sample(n = 200) 

#recombing this with the other data

dfHLA_rest <- dfHLA %>%
  filter(Species != "Homo sapiens") 

dfHLA <- bind_rows(dfHLA_sub, dfHLA_rest)
  
#looking at the new snumber of sequences per species 

sequence_counts_HLA_sub <- sequence_counts(dfHLA)
View(sequence_counts_HLA_sub)
  
  
# I would maybe add something in here to visualize the sequence lengths to make sure that they are all similar and that you don't have any major issues with the data! 

# visualizing the sequence lengths with a histogram

hist(x = dfCOI$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #most are between 0 and 1000 
hist(x = dfHLA$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences") #this one looks good! most of the sequences are between 200 and 1600 base pairs 

 
# Converting sequences to DNAStringSet format and name them by processid
class(dfCOI$Sequence)
dfCOI$nucleotides2 <- DNAStringSet(dfCOI$Sequence)
class(dfCOI$nucleotides2)

class(dfHLA$Sequence2)
dfHLA$nucleotides2 <- DNAStringSet(dfHLA$Sequence2)
class(dfHLA$nucleotides2)

#### - PART 2 - Finding the centroid seqeunce ----

# Creating the function to find the centroid

calculate_centroid <- function(seqs) {
  # performing multiple sequence alignment
  alignment <- DNAStringSet(muscle::muscle(seqs, gapOpening = -3000))
  
  # converting the alignment to a DNAbin object to calculate distances
  alignment_dnabin <- as.DNAbin(alignment)
  
  # calculating pairwise distance matrix using TN93 model
  dist_matrix <- dist.dna(as.DNAbin(alignment), model = "TN93")
  
  # calculating centroid (sequence with the lowest sum of pairwise distances)
  centroid_index <- which.min(rowSums(as.matrix(dist_matrix)))
  
  # extracting the sequence of the centroid
  centroid_sequence <- as.character(alignment[centroid_index])
  
  return(centroid_sequence)
}

# Grouping sequences by species for COI 

grouped_sequences_COI <- split(dfCOI$nucleotides2, dfCOI$Common_Name)

# Grouping species for HLA 

grouped_sequences_HLA <- split(dfHLA$nucleotides2, dfHLA$Common_Name)

# Apply the centroid calculation for each specimens and save it as a data frame

centroids_COI <- lapply(grouped_sequences_COI, calculate_centroid)
class(centroids_COI)

# Apply the centroid calculation for each specimens and save it as a data frame

centroids_HLA <- lapply(grouped_sequences_HLA, calculate_centroid)
class(centroids_HLA)

# Convert centroids to a DNAStringSet
centroids_COIset <- DNAStringSet(sapply(centroids_COI, DNAString))

centroids_HLAset <- DNAStringSet(sapply(centroids_HLA, DNAString))

# Aligning COI sequences using MUSCLE
cat("Aligning COI sequences...\n")
dfCOI.alignment <- DNAStringSet(muscle::muscle(centroids_COIset))

# Aligning HLA sequences using MUSCLE
cat("Aligning HLA sequences...\n")
dfHLA.alignment <- DNAStringSet(muscle::muscle(centroids_HLAset))

# Inspect the alignments
print(dfCOI.alignment)  # Shows aligned COI sequences
print(dfRAG1.alignment)  # Shows aligned HLA sequences

#### - PART 3 - phlygenetic analysis ----

# Building Phylogenetic Trees for COI and HLA
# Convert aligned sequences to phyDat format for tree-building

coi_phyDat <- phyDat(as.matrix(dfCOI.alignment), type = "DNA")

rag1_phyDat <- phyDat(as.matrix(dfRAG1.alignment), type = "DNA")

# Neighbor-Joining Tree for COI

cat("Building Neighbor-Joining Tree for COI...\n")
coi_dist <- dist.ml(coi_phyDat, model = "F81")           # Calculate distance matrix
coi_nj_tree <- NJ(coi_dist)                # Neighbor-Joining tree
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# I would preform a model test to determine the best model for your data this way you can determine which model works best! 

modelTest_coi <- modelTest(coi_phyDat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
View(modelTest_coi)

# Choosing the best model based on AIC since the sample size is small this is the better model to follow 
best_model_coi <- modelTest_coi[which.min(modelTest_coi$AIC), "Model"]
print(paste("Best modelAIC:", best_model_coi)) #GTR is the best model 

# Maximum Likelihood Tree for COI

cat("Building Maximum Likelihood Tree for COI...\n")
coi_ml_tree <- pml(coi_nj_tree, coi_phyDat)
coi_ml_tree <- optim.pml(coi_ml_tree, model = "GTR")
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Neighbor-Joining Tree for RAG1

cat("Building Neighbor-Joining Tree for RAG1...\n")
rag1_dist <- dist.ml(rag1_phyDat, model = "F81")           # Calculate distance matrix
rag1_nj_tree <- NJ(rag1_dist)                # Neighbor-Joining tree
plot(rag1_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

#Again here is the model test to make sure that you have the best model for your data   

modelTest_rag1 <- modelTest(rag1_phyDat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
View(modelTest_rag1)

# Choosing the best model based on AIC since the sample size is small this is the better model to follow 
best_model_rag1AIC <- modelTest_rag1[which.min(modelTest_rag1$AIC), "Model"]
print(paste("Best model_rag1AIC:", best_model_rag1AIC)) #GTR is the best model 

# Maximum Likelihood Tree for RAG1

cat("Building Maximum Likelihood Tree for RAG1...\n")
rag1_ml_tree <- pml(rag1_nj_tree, rag1_phyDat)
rag1_ml_tree <- optim.pml(rag1_ml_tree, model = "GTR")
plot(rag1_ml_tree$tree, main = "Maximum Likelihood T")

#Comparing Tree Topologies for COI and HLA
# Here, we calculate a distance metric between the COI and HLA trees

cat("Comparing COI and RAG1 tree topologies...\n")
tree_distance_coi <- RF.dist(coi_nj_tree, rag1_nj_tree) # 16 and the maximum value for unrooted trees is (2×(13−1)=24) which indicates a relatively large discrepancy between trees
print(tree_distance_coi)
cat("Robinson-Foulds distance between COI and HLA trees: ", tree_distance_coi, "\n")

#I would also do this calculation for the maximum likelihood trees 

tree_distance_ml <- RF.dist(coi_ml_tree$tree, rag1_ml_tree$tree) # 16 and the maximum value for unrooted trees is (2×(13−1)=24) which indicates a relatively large discrepancy between trees
print(tree_distance_ml)
cat("Robinson-Foulds distance between COI and RAG1 trees: ", tree_distance_ml, "\n")

# I would add in boot strapping here so that you can determine which tree best represents the data on top of determining that they are different 

bs_coi <- bootstrap.pml(coi_ml_tree, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

bs_rag1 <- bootstrap.pml(rag1_ml_tree, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

plotBS(midpoint(coi_ml_tree$tree), bs_coi, type = "phylogram", cex = 0.6)

plotBS(midpoint(rag1_ml_tree$tree), bs_rag1, type = "phylogram", cex = 0.6)

# Visualizing Trees Side-by-Side for Comparison
par(mfrow = c(2, 2))  # Set up 2x2 plotting area

# Plot 1: Neighbor-Joining Tree for COI
plot(coi_nj_tree, main = "Neighbor-Joining Tree for COI Sequences")

# Plot 2: Maximum Likelihood Tree for COI
plot(coi_ml_tree$tree, main = "Maximum Likelihood Tree for COI Sequences")

# Plot 3: Neighbor-Joining Tree for HLA
plot(rag1_nj_tree, main = "Neighbor-Joining Tree for HLA Sequences")

# Plot 4: Maximum Likelihood Tree for HLA
plot(rag1_ml_tree$tree, main = "Maximum Likelihood Tree for HLA Sequences")

# Save alignments and tree plots for future reference and analysis
writeXStringSet(dfCOI.alignment, file = "COI_aligned_sequences.fasta")
writeXStringSet(dfRAG1.alignment, file = "RAG1_aligned_sequences.fasta")
saveRDS(coi_nj_tree, file = "COI_NJ_Tree.rds")
saveRDS(coi_ml_tree, file = "COI_ML_Tree.rds")
saveRDS(rag1_nj_tree, file = "RAG1_NJ_Tree.rds")
saveRDS(rag1_ml_tree, file = "RAG1_ML_Tree.rds")

cat("Analysis complete. Results saved.\n")

HLAgene <- read.FASTA("HLAortho.FASTA")
