- Script Enhancements:  
  - Added code for downloading sequence data for both COI and HLA genes at the end of the script, enabling reproducibility and seamless integration into future scripts.  
  - Improved data quality through the addition of filtering checks.  

- Data Filtering Improvements:  
  - Implemented functions to filter NCBI data by removing low-quality or unnecessary sequences.  
  - Added functionality to take random subsets of overrepresented species, ensuring balanced representation.  
  - Included filtering checks for sequence length, removing sequences that are too long or too short.  

- Phylogenetic Analysis Improvements:  
  - Tested the best model for determining the maximum likelihood tree.  
  - Added bootstrapping to compare the HLA and COI genes for better data representation in phylogeny.  
  - Produced phylogenetic trees with bootstrap values to evaluate statistical support.  
  - Introduced a tanglegram for direct visual comparison of the dendrograms generated from HLA and COI genes.  

These changes collectively enhanced the script’s reproducibility, data quality, and analytical rigor.  
