#clear_environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR")
#color palette
colors <- c(rgb(100/255, 170/255, 112/255),
            rgb(119/255, 106/255, 105/255),
            rgb(210/255, 208/255, 213/255),
            rgb(147/255, 213/255, 231/255),
            rgb(239/255, 234/255, 183/255),
            rgb(234/255, 166/255, 123/255),
            rgb(212/255, 191/255, 153/255),
            rgb(202/255, 154/255, 129/255))
#import mtags genus level data
mtags_in <- readr::read_tsv("data/mtags/all.genus.tsv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(genus_id = paste0("genus_", as.character(dplyr::row_number())), .before = number_taxpath)
#create genus table
mtags_genera <- mtags_in %>%
  select(-number_taxpath) %>%
  tidyr::pivot_longer(-genus_id, names_to = "sample", values_to = "count") %>%
  dplyr::mutate(sample = stringr::str_remove(sample, pattern = "_bins") %>%
                  stringr::str_to_upper()) %>%
  tidyr::pivot_wider(names_from = "genus_id", values_from = "count")  %>%
  tibble::column_to_rownames(var = "sample")
#create taxonomy tabel
mtags_tax <- mtags_in %>%
  dplyr::select(genus_id, number_taxpath) %>%
  dplyr::mutate(root = stringr::word(number_taxpath, sep = ";", start = 1, end = 1),
                domain = stringr::word(number_taxpath, sep = ";", start = 2, end = 2),
                phylum = stringr::word(number_taxpath, sep = ";", start = 3, end = 3),
                class = stringr::word(number_taxpath, sep = ";", start = 4, end = 4),
                order = stringr::word(number_taxpath, sep = ";", start = 5, end = 5),
                family = stringr::word(number_taxpath, sep = ";", start = 6, end = 6),
                genus = stringr::word(number_taxpath, sep = ";", start = 7, end = 7)) %>%
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
  dplyr::filter(sample_code %in% rownames(mtags_genera))
#select only prokaryotes
mtags_tax <- mtags_tax %>%
  dplyr::filter(domain %in% c("Bacteria", "Archaea")|genus %in% c("Unassigned", "Unaligned"))
mtags_genera <- mtags_genera %>%
  dplyr::select(dplyr::any_of(mtags_tax$genus_id))
#keep exp2 only
sample_metadata <- sample_metadata %>%
  dplyr::filter(experiment == "2")
mtags_genera <- mtags_genera %>%
  dplyr::filter(rownames(mtags_genera) %in% sample_metadata$sample_code)  %>% 
  dplyr::select_if(colSums(.) != 0)
mtags_tax <- mtags_tax %>%
  dplyr::filter(genus_id %in% colnames(mtags_genera))
#remove unaligned reads if desired
remove_unaligned <- TRUE
if(remove_unaligned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(genus != "Unaligned")
  mtags_genera <- mtags_genera %>%
    dplyr::select(dplyr::any_of(mtags_tax$genus_id))
}
#remove unassigned reads if desired
remove_unassigned <- TRUE
if(remove_unassigned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(genus != "Unassigned")
  mtags_genera <- mtags_genera %>%
    dplyr::select(dplyr::any_of(mtags_tax$genus_id))
}
#TMM normalization
normalize <- TRUE
#Format as DGElist for normalization
cts <- mtags_genera %>%
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
  mtags_genera <- d2 %>%
    t() %>%
    as.data.frame()
}

#perform hellinger transformation if desired
hellinger <- FALSE
if(hellinger == TRUE){
  mtags_genera <- mtags_genera %>%
    vegan::decostand(method = "hellinger")
}
#perform square root transformation if desired
square_root <- TRUE
if(square_root == TRUE){
  mtags_genera <- mtags_genera %>%
    dplyr::mutate(across(.cols = everything(), sqrt))
}
#compute distance matrix
dist_method <- "bray"
dist_mat <- mtags_genera %>%
  vegan::vegdist(method = dist_method)
#compute pcoa
pcoa <- ape::pcoa(dist_mat)
write_rds(pcoa, "data/mtags_prokaryotes_exp2_pcoa_ape_output.rds")
#create data for plotting pcoa
pcoa_plot_dat <- tibble::tibble(sample = rownames(pcoa$vectors), axis1 = pcoa$vectors[,1], axis2 = pcoa$vectors[,2]) %>%
  dplyr::left_join(sample_metadata, by = c("sample" = "sample_code")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", "WW30", "WW80", "WW30UF", "WW80UF")),
                experiment = as.factor(experiment))
pcoa_plot_dat <- pcoa_plot_dat %>%
  dplyr::mutate(domain = rep("Prokaryota", nrow(pcoa_plot_dat)))
readr::write_rds(pcoa_plot_dat, "data/mtags_prokaryotes_exp2_pcoa.rds")
#plot
pcoa_plot <- ggplot2::ggplot(pcoa_plot_dat, aes(x = axis1, y = axis2, color = sample_perc)) +
  ggplot2::geom_point() +
  ggplot2::stat_ellipse() +
  ggplot2::labs(x = paste0("PCoA Axis 1 (", as.character(round(pcoa$values$Relative_eig[1], 4)*100), " %)"), y = paste0("PCoA Axis 2 (", as.character(round(pcoa$values$Relative_eig[2], 4)*100), " %)"), color = "Treatment") +
  ggpubr::theme_pubr() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::scale_color_manual(labels = c("Stream water", "Wastewater", "Wastewater UF", "0% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF"), values = colors)
#save plot
ggplot2::ggsave("output/figures/figureS3.jpg",
                pcoa_plot,
                device = "jpg",
                dpi = 300,
                units = "cm",
                width = 20,
                height = 12.5)
stat_test <- adonis2(dist_mat ~ sample_perc, sample_metadata, permutations = 999, by = "margin")
