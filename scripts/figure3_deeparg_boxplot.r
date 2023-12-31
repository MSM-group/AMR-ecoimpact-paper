rm(list = ls()) 
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR")
#color palette
colors <- c(rgb(100/255, 170/255, 112/255),
            rgb(119/255, 106/255, 105/255),
            rgb(210/255, 208/255, 213/255),
            rgb(147/255, 213/255, 231/255),
            rgb(254/255, 203/255, 103/255),
            rgb(239/255, 234/255, 183/255),
            rgb(234/255, 166/255, 123/255),
            rgb(212/255, 191/255, 153/255),
            rgb(202/255, 154/255, 129/255))
#DeepARGs main types (classes)
deepargs_dat <- readr::read_csv("data/deeparg_type_16S_norm.csv") %>%
  janitor::clean_names()
#import sample metadata 
sample_metadat <- readxl::read_excel("data/metadata/EcoImpact_Exp1_Exp2_DNA_samples_LC_2_metadata.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(sample_perc = dplyr::case_when(grepl("BT", sample_name) ~ paste0("BT_", stringr::word(sample_name, 2, sep = "_")),
                                               TRUE ~ stringr::word(sample_name, 1, sep = "_"))) %>%
  dplyr::mutate(timefix = dplyr::case_when(grepl("Week4", time) ~ "D28",  # make times consistent
                                           grepl("Week3", time) ~ "D21",
                                           grepl("Week2", time) ~ "D14",
                                           grepl("Week1", time) ~ "D07",
                                           grepl("Day", time) ~ gsub("Day", "D", time),
                                           TRUE ~ time)) %>%
  dplyr::select(sample_code, experiment, sample_type, sample_perc) %>%
  dplyr::filter(sample_code %in% deepargs_dat$sample)%>%
  dplyr::rename(., sample = sample_code)

data <- dplyr::left_join(deepargs_dat, sample_metadat, by= "sample")%>%
  pivot_longer(cols= 2:41, names_to = "args", values_to = "read_count")

#Create labels for p-values
p_adj_stars = function(p){
  case_when(p > 0.05 ~ "ns",
            p <= 0.05 & p > 0.01 ~ "*",
            p <= 0.01 & p > 0.001 ~ "**",
            p <= 0.001 & p > 0.0001 ~ "***",
            p <= 0.0001 ~ "****")
}

# P-value correction (of the Kruskal-Wallis test) using Benjamini-Hochberg correction 
stati_test <-  data %>%
  dplyr::group_by(args) %>%
  rstatix::kruskal_test(read_count ~ sample_perc) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  arrange(p.adj) %>%
  mutate(lab = paste0(method, ": ", p_adj_stars(p.adj)))

p_corrected <- dplyr::left_join(data, stati_test)

write.csv(p_corrected, "output/stat_tests/20230725_deepargs_main_types_kruskal_wallistest_corrected_P_BHmethod.csv")
sig_deepargs <- p_corrected %>%
  filter(p.adj < 0.05) %>%
  select("sample", "sample_perc", "args", "read_count", "p.adj", "lab", "sample_type", "experiment")
write.csv(sig_deepargs, "output/stat_tests/20230725_sig_deepargs_maintypes_corrected_BH_kruskal_wallis.csv")

#tidying and filtering data
plot_data<- read_csv("output/stat_tests/20230725_sig_deepargs_maintypes_corrected_BH_kruskal_wallis.csv") %>%
  dplyr::mutate(read_count_corrected = read_count+1) %>%
  filter(as.numeric(read_count) > 0.001)%>%
  filter(experiment == 2)
plot_data2 <- plot_data %>%
  filter(args %in% c("aminoglycoside", "beta_lactam", "fluoroquinolone", "fosfomycin", "fosmidomycin", "glycopeptide", 
                     "mls", "multidrug", "peptide", "Phenicol", "rifamycin", "sulfonamide", "tetracycline")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", "WW30", "WW80", "WW30UF", "WW80UF")))

#creating labels
x_axis_labels<- c("Aminoglycoside", "Beta-lactam", "Fluoroquinolone", "Fosfomycin", "Fosmidomycin", "Glycopeptide", 
                  "MLS", "Multidrug", "Peptide", "Rifamycin", "Sulfonamide", "Tetracycline")  


plot_data3 <- plot_data %>%
  filter(args %in% c("aminoglycoside", "beta_lactam", "multidrug", "sulfonamide","tetracycline")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("BT_CB", "BT_WW", "BT_UF", "WW00", 
                                                                  "WW30", "WW80", "WW30UF", "WW80UF"))) 
arg_types<- c("aminoglycoside" = "Aminoglycoside", 
              "beta_lactam"= "Beta-lactam", 
              "fluoroquinolone"= "Fluoroquinolone",
              "fosfomycin" = "Fosfomycin",
              "fosmidomycin" = "Fosmidomycin",
              "glycopeptide" = "Glycopeptide",
              "mls" = "MLS",
              "peptide" = "Peptide",
              "rifamycin" = "Rifamycin",
              "multidrug"= "Multidrug",
              "sulfonamide"= "Sulfonamide", 
              "tetracycline"= "Tetracycline")

#filtering out water samples
plot_data_no_water <- plot_data2 %>%
  filter(!(sample_perc %in% c("BT_CB", "BT_WW", "BT_UF")))

stati_test2 <-  plot_data_no_water %>%
  dplyr::group_by(args) %>%
  rstatix::kruskal_test(read_count ~ sample_perc) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  dplyr::arrange(p.adj)

p_corrected2 <- dplyr::left_join(plot_data_no_water, stati_test2)

plot_data_no_water2 <- plot_data_no_water %>%
  filter(args %in% c("aminoglycoside", "beta_lactam", "fluoroquinolone", "glycopeptide", 
                     "mls", "multidrug", "peptide", "Phenicol", "rifamycin", "sulfonamide", "tetracycline"))
write.csv(plot_data_no_water2, "output/stat_tests/20230725_most_sig_deepargs_maintypes_corrected_BH_kruskal_wallis.csv")  

#plot
plot1 <- ggplot(plot_data_no_water2, aes(fill=sample_perc, y=as.numeric(read_count_corrected), x= sample_perc))+
  geom_boxplot() +
  geom_jitter(pch=21, color = "black", aes(x = sample_perc)) +
  scale_colour_manual(values = c("WW00"="#93D5E7",
                                 "WW30"="#EFEAB7",
                                 "WW80"="#EAA67B",
                                 "WW30UF"="#D4BF99",
                                 "WW80UF"="#CA9A81"))+
  scale_fill_manual(name = NULL, values = c("WW00"="#93D5E7",
                                            "WW30"="#EFEAB7",
                                            "WW80"="#EAA67B",
                                            "WW30UF"="#D4BF99",
                                            "WW80UF"="#CA9A81"),
                    labels=c("WW00" = "0% WW",
                             "WW30" = "30% WW",
                             "WW80" = "80% WW",
                             "WW30UF" = "30% WW UF",
                             "WW80UF" = "80% WW UF"))+
  scale_x_discrete(labels=c("WW00" = "0% WW", 
                            "WW30" = "30% WW",
                            "WW80" = "80% WW",
                            "WW30UF" = "30% WW UF",
                            "WW80UF" = "80% WW UF"))+
  stat_pvalue_manual(plot_data_no_water2, label = "lab")
  theme_pubr(legend = c("right")) +
  facet_wrap(nrow= 4, ncol= 3, vars(args), scales = "free_y", labeller = as_labeller(arg_types)) +
  labs(x = "Condition", y = "16S normalized read counts (Corrected +1)")+
  theme(legend.position= "bottom",
        legend.spacing= unit(7, 'mm'), 
        legend.key.size = unit(7, 'mm'),
        legend.text = element_text(size=11),
        legend.title = element_text(size=14, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="plain", size=10, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face="plain", size=10, colour = "black"), strip.text = element_text(size=14, face="bold"))
ggsave("output/figures/figure3_deeparg_boxplot.png", 
       plot1, device= "png", units= c("mm"), height = 240, width = 190, dpi = 500)
