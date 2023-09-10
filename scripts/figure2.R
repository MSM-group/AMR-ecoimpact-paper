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
#import deepargs data
deepargs_in <- readr::read_csv("data/deeparg_subtype_16S_norm.csv") %>%
  janitor::clean_names()

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
  dplyr::filter(sample_code %in% deepargs_in$sample)
deepargs <- deepargs_in %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "sample")

#square root transformation
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
write_rds(pcoa, "data/deepargs_pcoa.rds")

#compute PREMANOVA
permanova <- adonis2(dist_mat ~ sample_perc + experiment, sample_metadata, permutations = 999, by = "margin")

#create data for plotting pcoa
pcoa_plot_dat <- tibble::tibble(sample = rownames(pcoa$vectors), axis1 = pcoa$vectors[,1], axis2 = pcoa$vectors[,2]) %>%
  dplyr::left_join(sample_metadata, by = c("sample" = "sample_code")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", "WW30", "WW80", "WW30UF", "WW80UF")),
                experiment = as.factor(experiment))

pcoa_plot_dat2 <- pcoa_plot_dat %>%
  filter(experiment == 2) #%>%


#plot
pcoa_plot <- ggplot2::ggplot(pcoa_plot_dat2, aes(x = axis1, y = axis2, color = sample_perc)) +
  ggplot2::geom_point(aes(x=axis1,y=axis2,colour=sample_perc),alpha = 1, size=1.5) +
  #ggplot2::geom_line()+
  ggplot2::labs(x = paste0("PCoA Axis 1 (", as.character(round(pcoa$values$Relative_eig[1]*100, 1)), " %)"), y = paste0("PCoA Axis 2 (", as.character(round(pcoa$values$Relative_eig[2]*100, 1)), " %)"), color = "Treatment") +
  ggplot2::guides(color = guide_legend(order=2),
                  shape = guide_legend(order=1))+
  ggpubr::theme_pubr() +
  ggplot2::theme(legend.position= "right", 
                 legend.spacing= unit(7, 'mm'), 
                 legend.key.size = unit(0.5, 'mm'),
                 legend.text = element_text(size=7),
                 legend.title = element_text(size=9, face="plain"),
                 axis.title=element_text(size=8,face="bold"),
                 axis.text.x = element_text(face="plain", size=8, colour = "black", hjust = 1),
                 axis.text.y = element_text(face="plain", size=8, colour = "black")) +
  ggplot2::scale_colour_manual(values = colors,
                               labels = c("River water", "Waste water", "Waste water UF", "0% WW", "30% WW", "80% WW", "30% WW UF", "80% WW UF"))+
  ggplot2::stat_ellipse(linewidth = 0.3)

ggplot2::ggsave("paper_figures/figure2.png",
                pcoa_plot, device= "png", units= c("mm"), height = 70, width = 120, dpi = 500)

