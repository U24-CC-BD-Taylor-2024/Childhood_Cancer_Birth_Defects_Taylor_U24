# R Corbett 2025
#
# Obtain HGVSg IDs from MTP variant table

# Load libraries
library(tidyverse)

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "preprocessing",
                          "allele-registry")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

# Define file paths
variant_file <- file.path(input_dir, 
                       "variant-level-snv-consensus-annotated-mut-freq.tsv.gz")

chr_key_file <- file.path(input_dir,
                          "chr_to_refseq.tsv")

# Wrangle data
variant_df <- read_tsv(variant_file) %>%
  # retain only unique variants
  distinct(Variant_ID_hg38, .keep_all = TRUE)

chr_key <- read_tsv(chr_key_file)

# write function to replace chr portion HGVSg ID with NC accession numbers
convert_hgvs <- function(hgvs_vector, ref_table) {
  # Iterate over each chromosome and replace it with RefSeq ID
  for (i in 1:nrow(ref_table)) {
    
    if (i == 1){
    
      new_hgvsg <- gsub(paste0("^chr", ref_table$Chromosome[i], ":"), 
                          paste0(ref_table$GRCh38[i], ":"), 
                          hgvs_vector)
      
    } else {
      
      new_hgvsg <- gsub(paste0("^chr", ref_table$Chromosome[i], ":"), 
                        paste0(ref_table$GRCh38[i], ":"), 
                        new_hgvsg)
      
    }
    
  }
  
  return(new_hgvsg)

}

# define cohorts included in variant df
cohorts <- c("GMKF", "PBTA", "TARGET")

# loop through cohorts, convert chr IDs to Refseq chr accession numbers in HGVSg, and write to output
for (cohort in cohorts){
  
  # subset variant df for cohort
  cohort_df <- variant_df %>%
    dplyr::filter(grepl(cohort, Dataset))
  
  # convert variant HGVSg ids
  cohort_df <- cohort_df %>%
    dplyr::mutate(HGVSg_clingen = convert_hgvs(cohort_df$HGVSg, chr_key))
  
  # write variants to output
  write_tsv(as.data.frame(unique(cohort_df$HGVSg_clingen)),
            file.path(results_dir, glue::glue("{cohort}-variants.tsv")),
            col_names = FALSE)
  
}

sessionInfo()
