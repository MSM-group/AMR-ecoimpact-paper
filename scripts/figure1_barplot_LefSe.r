#clear_environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR", "RColorBrewer", "gplots")
#import mtags phylum level data
mtags_in <- readr::read_tsv("data/mtags/all.phylum.tsv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(phylum_id = paste0("phylum_", as.character(dplyr::row_number())), .before = number_taxpath)
#create phylum table
mtags_phyla <- mtags_in %>%
  select(-number_taxpath) %>%
  tidyr::pivot_longer(-phylum_id, names_to = "sample", values_to = "count") %>%
  dplyr::mutate(sample = stringr::str_remove(sample, pattern = "_bins") %>%
                  stringr::str_to_upper()) %>%
  tidyr::pivot_wider(names_from = "phylum_id", values_from = "count")  %>%
  tibble::column_to_rownames(var = "sample")
#create taxonomy tabel
mtags_tax <- mtags_in %>%
  dplyr::select(phylum_id, number_taxpath) %>%
  dplyr::mutate(root = stringr::word(number_taxpath, sep = ";", start = 1, end = 1),
                domain = stringr::word(number_taxpath, sep = ";", start = 2, end = 2),
                phylum = stringr::word(number_taxpath, sep = ";", start = 3, end = 3)) %>%
  dplyr::select(-number_taxpath) %>%
  dplyr::mutate(dplyr::across(.cols = !ends_with("_id"), stringr::word, sep = "__", start = 2, end = 2))
mtags_tax[nrow(mtags_tax)-1, ncol(mtags_tax)] <- dplyr::pull(mtags_in, number_taxpath)[nrow(mtags_tax)-1]
mtags_tax[nrow(mtags_tax), ncol(mtags_tax)] <- dplyr::pull(mtags_in, number_taxpath)[nrow(mtags_tax)]
#import sample metadata 
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
  dplyr::filter(sample_code %in% rownames(mtags_phyla))
#keep only prokaryotes
mtags_tax <- mtags_tax %>%
  dplyr::filter(domain != "Eukaryota")
mtags_phyla <- mtags_phyla %>%
  dplyr::select(all_of(mtags_tax$phylum_id))
#keep exp2 only
sample_metadata <- sample_metadata %>%
  dplyr::filter(experiment == "2")
mtags_phyla <- mtags_phyla %>%
  dplyr::filter(rownames(mtags_phyla) %in% sample_metadata$sample_code)  %>% 
  dplyr::select_if(colSums(.) != 0)
mtags_tax <- mtags_tax %>%
  dplyr::filter(phylum_id %in% colnames(mtags_phyla))
#keep biofilm samples only
sample_metadata <- sample_metadata %>%
  dplyr::filter(sample_type == "Biofilm")
mtags_phyla <- mtags_phyla %>%
  dplyr::filter(rownames(mtags_phyla) %in% sample_metadata$sample_code)  %>% 
  dplyr::select_if(colSums(.) != 0)
mtags_tax <- mtags_tax %>%
  dplyr::filter(phylum_id %in% colnames(mtags_phyla))
#remove unaligned reads if desired
remove_unaligned <- FALSE
if(remove_unaligned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(phylum != "Unaligned")
  mtags_phyla <- mtags_phyla %>%
    dplyr::select(dplyr::any_of(mtags_tax$phylum_id))
}
#remove unassigned reads if desired
remove_unassigned <- FALSE
if(remove_unassigned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(phylum != "Unassigned")
  mtags_phyla <- mtags_phyla %>%
    dplyr::select(dplyr::any_of(mtags_tax$phylum_id))
}
#TMM normalization
normalize <- TRUE
#Format as DGElist for normalization
cts <- mtags_phyla %>%
  t()
grps <- sample_metadata %>%
  dplyr::select(sample_code, sample_perc) %>%
  arrange(desc(sample_code)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("sample_code") %>%
  t()

d <- edgeR::DGEList(counts = cts, group = grps)
#normalize using tmm method
d1 <- edgeR::calcNormFactors(d, method = "TMM")
d2 <- edgeR::cpm(d1)
if(normalize == TRUE){
  mtags_phyla <- d2 %>%
    t() %>%
    as.data.frame()
}
d2

# read in the phylum
phyla_pal <- read_csv("data/phylum_palette.csv")
colnames(phyla_pal) <- c("phylum", "color")
phyla_pal
#bin low abundance taxa as other
mtags_phyla_deciles  <- mtags_phyla %>%
  rownames_to_column(var = "sample") %>%
  tidyr::pivot_longer(-sample, names_to = "phylum_id", values_to = "count") %>%
  dplyr::left_join(sample_metadata, by = c("sample" = "sample_code")) %>%
  dplyr::group_by(phylum_id, sample_perc) %>%
  dplyr::summarise(count_sum = sum(count)) %>%
  dplyr::group_by(sample_perc) %>%
  dplyr::mutate(decile = dplyr::ntile(count_sum, 1)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(phylum_id) %>%
  dplyr::summarise(max_decile = max(decile))
upper_decile <- mtags_phyla_deciles %>%
  dplyr::filter(max_decile == 1) %>%
  dplyr::pull(phylum_id)
mtags_phyla_upper_decile <- mtags_phyla %>%
  dplyr::select(all_of(upper_decile))
mtags_phyla_other <- mtags_phyla %>%
  dplyr::select(!all_of(upper_decile)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(other = sum(c_across(where(is.numeric)))) %>%
  dplyr::select(other)
mtags_phyla <- mtags_phyla_upper_decile %>%
  dplyr::mutate(other = dplyr::pull(mtags_phyla_other, other))
mtags_phyla

#get mean abundance
mtags_phyla_mean <- mtags_phyla %>%
  rownames_to_column(var = "sample") %>%
  dplyr::left_join(sample_metadata, by = c("sample" = "sample_code")) %>%
  dplyr::select(-c(timefix, experiment, sample_type)) %>%
  tidyr::pivot_longer(-c(sample, sample_perc), names_to = "phylum_id", values_to = "count") %>%
  dplyr::group_by(sample_perc, phylum_id) %>%
  dplyr::summarise(mean_count = mean(count)) %>%
  dplyr::group_by(sample_perc) %>%
  dplyr::mutate(rel = mean_count/sum(mean_count)) %>%
  dplyr::left_join(mtags_tax, by = "phylum_id") %>%
  dplyr::select(-c(root, domain)) %>%
  #dplyr::filter(phylum %in% c(phyla_pal$phylum[!is.na(phyla_pal$phylum)])) %>%
  tidyr::replace_na(list(phylum = "other")) %>%
  dplyr::filter(phylum %in% c(phyla_pal$phylum[!is.na(phyla_pal$phylum)])) 
  
#color palette
n_colors <- length(levels(mtags_phyla_mean$phylum))-1

#get_palette <- colorRampPalette(brewer.pal(8, "Set1"))
set.seed(121)
#pal <- c(sample(phyla_pal$color, n_colors, replace = FALSE), gplots::col2hex("gray80"))
pal <- phyla_pal$color
mtags_phyla_mean$phylum

mtags_phyla_mean$phylum <- gsub("Firmicutes", "Bacillota", mtags_phyla_mean$phylum)
mtags_phyla_mean$phylum <- gsub("Proteobacteria", "Pseudomonadota", mtags_phyla_mean$phylum)
mtags_phyla_mean$phylum <- gsub("Chloroflexi", "Chloroflexota", mtags_phyla_mean$phylum)

mtags_phyla_mean <- mtags_phyla_mean %>%
  dplyr::mutate(sample_perc = case_when(sample_perc == "BT_CB" ~ "Stream water",
                                        sample_perc == "BT_WW" ~ "Wastewater",
                                        sample_perc == "BT_UF" ~ "Wastewater UF",
                                        sample_perc == "WW00" ~ "0% WW",
                                        #sample_perc == "WW10" ~ "10% WW",
                                        sample_perc == "WW30" ~ "30% WW",
                                        sample_perc == "WW80" ~ "80% WW",
                                        sample_perc == "WW30UF" ~ "30% WW UF",
                                        sample_perc == "WW80UF" ~ "80% WW UF") %>%
                  factor(levels = c("0% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF")))
levs <- sort(unique(mtags_phyla_mean$phylum))
mtags_phyla_mean$phylum <- factor(mtags_phyla_mean$phylum, levels = levs)

# make plot
plot <- ggplot2::ggplot(mtags_phyla_mean, aes(x = sample_perc, y = rel, fill = phylum)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") +
  #ggplot2::facet_grid(cols = vars(sample_type), scales = "free_x", space = "free") +
  #ggplot2::scale_x_discrete(labels = c("River water", "Wastewater", "Wastewater UF", "0% WW", "10% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF")) +
  ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%")) +
  ggplot2::scale_fill_manual(values = pal) +
  ggpubr::theme_pubr() +
  ggplot2::labs(x = "Biofilm treatment", y = "Relative abundance (%)", fill = "Phylum") +
  ggplot2::theme(legend.position = "right",
                 legend.text = element_text(size = 28, face = "italic"),
                 legend.title = element_text(size = 32),
                 axis.title = element_text(size = 32),
                 axis.text.x = element_text(size = 28, angle = 45, vjust = 1, hjust = 1),
                 axis.text.y = element_text(size = 28),
                 strip.text = element_text(size = 28))
plot


ggsave("output/figures/figure1.jpg",
       dpi = 300,
       device = "jpeg",
       units = "cm",
       width = 75/2.75,
       height = (75/2)/1.6)

