rm(list = ls()) #clear data
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape")
#import mtags genus level data
deepargs_in <- readr::read_csv("data/deeparg_subtype_16S_norm.csv") %>%
  janitor::clean_names()

#import sample metadata (based on Serina's code)
sample_metadata <- readxl::read_excel("data/metadata/EcoImpact_Exp1_Exp2_DNA_samples_LC_2_metadata.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(sample_perc = dplyr::case_when(grepl("BT", sample_name) ~ paste0("BT_", stringr::word(sample_name, 2, sep = "_")),
                                               TRUE ~ stringr::word(sample_name, 1, sep = "_"))) %>%
  dplyr::mutate(timefix = dplyr::case_when(grepl("Week4", time) ~ "D28",  # make times consistent
                                           grepl("Week3", time) ~ "D21",
                                           grepl("Week2", time) ~ "D14",
                                           grepl("Week1", time) ~ "D07",
                                           grepl("Day", time) ~ gsub("Day", "D", time),
                                           TRUE ~ time)) %>%
  dplyr::select(sample_code, experiment, sample_type, sample_perc, timefix) %>%
  dplyr::filter(sample_code %in% deepargs_in$sample) %>%
  dplyr::arrange(sample_code)
deepargs <- deepargs_in %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "sample")
#keep exp2 only
sample_metadata <- sample_metadata %>%
  dplyr::filter(experiment == "2")
deepargs <- deepargs %>%
  dplyr::filter(rownames(deepargs) %in% sample_metadata$sample_code)  %>% 
  dplyr::select_if(colSums(.) != 0)
#perform hellinger transformation if desired
hellinger <- FALSE
if(hellinger == TRUE){
  deepargs <- deepargs %>%
    vegan::decostand(method = "hellinger")
}
#perform square root transformation if desired
square_root <- TRUE
if(square_root == TRUE){
  deepargs <- deepargs %>%
    dplyr::mutate(across(.cols = everything(), sqrt))
}



#compute distance matrix
dist_method <- "bray"
dist_mat <- deepargs %>%
  vegan::vegdist(method = dist_method)

#compute pcoa
pcoa <- ape::pcoa(dist_mat)
write_rds(pcoa, "data/deepargs_pcoa_exp2.rds")