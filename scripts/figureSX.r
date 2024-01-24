#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR")

##color palette
colors <- c(rgb(100/255, 170/255, 112/255),
            rgb(119/255, 106/255, 105/255),
            rgb(210/255, 208/255, 213/255),
            rgb(147/255, 213/255, 231/255),
            rgb(239/255, 234/255, 183/255),
            rgb(234/255, 166/255, 123/255),
            rgb(212/255, 191/255, 153/255),
            rgb(202/255, 154/255, 129/255))

# Read in the diamond results
diamond <- readr::read_csv("data/Exp2_functional_gene_counts.csv") %>%
  janitor::clean_names() %>%
  remove_rownames %>% 
  column_to_rownames(var="sample")

#perform hellinger transformation if desired
hellinger <- TRUE
if(hellinger == TRUE){
  diamond <- diamond %>%
    vegan::decostand(method = "hellinger")
}

#perform square root transformation if desired
square_root <- FALSE
if(square_root == TRUE){
  diamond <- diamond %>%
    dplyr::mutate(across(.cols = everything(), sqrt))
}

#compute distance matrix
dist_method <- "bray"

dist_mat <- diamond %>%
  vegan::vegdist(method = dist_method)

#compute pcoa
pcoa <- ape::pcoa(dist_mat)

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
  dplyr::filter(experiment == 2)

pcoa_plot_dat <- tibble::tibble(sample = rownames(pcoa$vectors), axis1 = pcoa$vectors[,1], axis2 = pcoa$vectors[,2]) %>%
  dplyr::left_join(sample_metadata, by = c("sample" = "sample_code")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", "WW30", "WW80", "WW30UF", "WW80UF")),
                experiment = as.factor(experiment))

kegg_pcoa_plot <- ggplot2::ggplot(pcoa_plot_dat, aes(x = axis1, y = axis2, color = sample_perc)) +
  ggplot2::geom_point() +
  ggplot2::stat_ellipse() +
  ggplot2::labs(x = paste0("PCoA Axis 1 (", as.character(round(pcoa$values$Relative_eig[1], 4)*100), " %)"), y = paste0("PCoA Axis 2 (", as.character(round(pcoa$values$Relative_eig[2], 4)*100), " %)"), color = "Treatment") +
  ggpubr::theme_pubr() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::scale_color_manual(labels = c("River water", "Wastewater", "Wastewater UF", "0% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF"), values = colors)

ggplot2::ggsave("output/figures/figureSx.png",
                kegg_pcoa_plot, 
                device= "png", units= c("mm"), height = 70, width = 120, dpi = 500)

stat_test <- adonis2(dist_mat ~ sample_perc, sample_metadata, permutations = 999, by = "margin")
stat_test


