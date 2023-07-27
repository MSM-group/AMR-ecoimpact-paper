#clear_environment
rm(list = ls()) #clear data
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR")

#importing SARG data
sig_genes_types <- read.csv("data/20230119_sig_genes_types_biofilmsamples102orfs.csv") %>%
  janitor::clean_names()%>%
  mutate(p_value= signif(p_adj, digits = 2)) %>%
  mutate(p_value_corr= str_trim(p_value, side = c("both")))
round(p_adj$sig_genes_types, digits = 3)

sig_genes_types_mod <- sig_genes_types %>%
  mutate(contig_mod = recode(contig,
                             "k121_858589_5" = "Beta-lactam",
                             "k121_252136_5" = "Beta-lactam",
                             "k121_1269806_2" = "Beta-lactam",
                             "k121_895379_1" = "Other",
                             "k121_1259309_2" = "Other",
                             "k121_1519353_3" = "Other",
                             "k121_660679_2" = "Other",
                             "k121_68825_2" = "Aminoglycoside",
                             "k121_906279_2" = "Multidrug",
                             "k121_671215_2" = "Macrolide-Lincomide-Streptogramin")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(condition, c("WW00", "WW10", 
                                                                "WW30", "WW80", 
                                                                "WW30UF", "WW80UF")))
#creating labels
contig_type <- c("k121_858589_5" = "Beta-lactam1",
                 "k121_252136_5" = "Beta-lactam2",
                 "k121_1269806_2" = "Beta-lactam3",
                 "k121_895379_1" = "Other1",
                 "k121_1259309_2" = "Other2",
                 "k121_1519353_3" = "Other3",
                 "k121_660679_2" = "Other4",
                 "k121_68825_2" = "Aminoglycoside",
                 "k121_906279_2" = "Multidrug",
                 "k121_671215_2" = "MLS")

#filtering experiment 2
sig_genes_types_mod2 <- sig_genes_types_mod %>%
  filter(experiment == 2)%>%
  mutate(contig_mod = recode(contig,
                             "k121_858589_5" = "Beta-lactam 1",
                             "k121_252136_5" = "Beta-lactam 2",
                             "k121_1269806_2" = "Beta-lactam 3",
                             "k121_895379_1" = "Other1",
                             "k121_1259309_2" = "Other2",
                             "k121_1519353_3" = "Other3",
                             "k121_660679_2" = "Other4",
                             "k121_68825_2" = "Aminoglycoside",
                             "k121_906279_2" = "Multidrug",
                             "k121_671215_2" = "Macrolide-Lincomide-Streptogramin")) %>%
  dplyr::mutate(sample_perc = forcats::fct_relevel(sample_perc, c("WW00", "WW30", 
                                                                  "WW80", "WW30UF",
                                                                  "WW80UF"))) %>%
  dplyr::mutate(corrected_reads= normalized_reads+1)
sig_genes_types_mod3 <- sig_genes_types_mod2 %>%
  filter(contig_mod %in% c("Beta-lactam 1", "Beta-lactam 2", "Beta-lactam 3", "Aminoglycoside", "Multidrug", 
                           "Macrolide-Lincomide-Streptogramin")) %>%
  dplyr::mutate(contig_mod = forcats::fct_relevel(contig_mod, c("Beta-lactam 1", 
                                                                "Beta-lactam 2", 
                                                                "Beta-lactam 3", 
                                                                "Aminoglycoside", 
                                                                "Multidrug", 
                                                                "Macrolide-Lincomide-Streptogramin")))

#plot

plot1<- ggplot(sig_genes_types_mod3, aes(x= sample_perc, y= as.numeric(corrected_reads), 
                                               fill = sample_perc))+
  geom_boxplot() +
  geom_jitter(aes(x = sample_perc, as.numeric(corrected_reads), fill = sample_perc, color= sample_perc)) +
  scale_colour_manual(guide = "none", values = c("WW00"="#93D5E7",
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
  theme_pubr() +
  facet_wrap(vars(contig_mod), scales = "free_y", nrow = 2, ncol = 3) +
  labs(x = "Condition", y = "Normalized average coverage (corrected +1)")+
  theme(legend.position= "bottom",
        legend.spacing= unit(5, 'mm'), 
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, face="bold"),
        axis.title=element_text(size=9,face="bold"),
        axis.text.x = element_text(face="plain", size=8, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face="plain", size=8, colour = "black"), 
        strip.text = element_text(size=8, face="bold"))

ggsave("paper_figures/figure5.png", plot1, 
       device= "png", units= c("mm"), height = 120, width = 190, dpi = 500)
