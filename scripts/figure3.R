#clear_environment
rm(list = ls())
#color palette
colors <- c(rgb(100/255, 170/255, 112/255),
            rgb(119/255, 106/255, 105/255),
            rgb(210/255, 208/255, 213/255),
            rgb(147/255, 213/255, 231/255),
            rgb(239/255, 234/255, 183/255),
            rgb(234/255, 166/255, 123/255),
            rgb(212/255, 191/255, 153/255),
            rgb(202/255, 154/255, 129/255))
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR")
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
  dplyr::select(sample_code, experiment, sample_type, sample_perc, timefix) 
prokaryotes <- readr::read_rds("data/mtags_prokaryotes_exp2_pcoa_ape_output.rds")
deeparg <- readr::read_rds("data/deepargs_pcoa_exp2.rds")
crustypro <- vegan::procrustes(X = prokaryotes$vectors, Y = deeparg$vectors, symmetric = TRUE)
ctest <- tibble(sample = rownames(crustypro$X),
                yrda1 = crustypro$Yrot[,1],
                yrda2 = crustypro$Yrot[,2],
                xrda1 = crustypro$X[,1],
                xrda2 = crustypro$X[,2]) %>%
  left_join(sample_metadata, by = c("sample" = "sample_code"))  %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", "WW30", "WW80", "WW30UF", "WW80UF")),
                experiment = as.factor(experiment))
plot_gg <- ggplot(ctest) +
  #geom_point(aes(x=yrda1, y=yrda2, color = sample_perc), size = 5) +
  geom_point(aes(x=xrda1, y=xrda2, color = sample_perc), size = 5) +
  geom_segment(aes(x=xrda1,y=xrda2,xend=yrda1,yend=yrda2, color = sample_perc), arrow=arrow(length=unit(0.5,"cm")), linewidth = 1.5) +
  theme_pubr() +
  ggplot2::theme(legend.position = "right",
                 legend.text = element_text(size = 28),
                 legend.title = element_text(size = 32),
                 strip.text.x = element_text(size = 32, face = "italic"),
                 axis.title = element_text(size = 32),
                 axis.text = element_text(size = 28)) +
  ggplot2::scale_color_manual(labels = c("River water", "Wastewater", "Wastewater UF", "0% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF"), values = colors) +
  ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = "Treatment")
#save plot
ggsave("paper_figures/figure3.jpg",
       dpi = 300,
       device = "jpeg",
       units = "cm",
       width = 75/2,
       height = (75/2)/2.3)
#statistical test
stat_test <- protest(X = prokaryotes$vectors, Y = deeparg$vectors, scores = "sites", permutations = 999)
