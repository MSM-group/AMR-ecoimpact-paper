#clear_environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr", "scales",
               "janitor", "ape", "edgeR", "phyloseq", "mia", "lefser", 
               "phyloseqCompanion", "microeco", "file2meco", "ggtree", "RColorBrewer")

#import mtags genus level data
mtags_in <- readr::read_tsv("data/mtags/all.genus.tsv") %>%
  janitor::clean_names() 
genus_tab <- mtags_in %>%
  dplyr::select(-number_taxpath) %>%
  phyloseq::otu_table(taxa_are_rows = TRUE)
colnames(genus_tab) <- colnames(genus_tab) %>%
  stringr::str_remove(pattern = "_bins") %>%
  stringr::str_to_upper()
rownames(genus_tab) <- rownames(genus_tab) %>%
  stringr::str_replace(pattern = "sp", replacement = "Genus")
tax_tab <- mtags_in %>%
  dplyr::select(number_taxpath)  %>%
  dplyr::mutate(root = stringr::word(number_taxpath, sep = ";", start = 1, end = 1),
                domain = stringr::word(number_taxpath, sep = ";", start = 2, end = 2),
                phylum = stringr::word(number_taxpath, sep = ";", start = 3, end = 3),
                class = stringr::word(number_taxpath, sep = ";", start = 4, end = 4),
                order = stringr::word(number_taxpath, sep = ";", start = 5, end = 5),
                family = stringr::word(number_taxpath, sep = ";", start = 6, end = 6),
                genus = stringr::word(number_taxpath, sep = ";", start = 7, end = 7)) %>%
  dplyr::select(-number_taxpath) %>%
  dplyr::mutate(dplyr::across(.cols = !ends_with("_id"), stringr::word, sep = "__", start = 2, end = 2)) %>%
  dplyr::select(-root) %>%
  as.matrix() %>%
  phyloseq::tax_table()
colnames(tax_tab) <- colnames(tax_tab) %>%
  stringr::str_to_title()
rownames(tax_tab) <- rownames(tax_tab) %>%
  stringr::str_replace(pattern = "sp", replacement = "Genus")

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
  dplyr::filter(sample_code %in% colnames(genus_tab)) %>%
  dplyr::arrange(sample_code) %>%
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$sample_code
sample_metadata <- sample_metadata %>%
  dplyr::select(-sample_code) %>%
  dplyr::mutate(sample_group = dplyr::case_when(stringr::str_detect(sample_perc, "UF") ~ "30-80 % WW UF",
                                                sample_perc == "WW00" ~ "0% WW",
                                                TRUE ~ "30-80 % WW"))
#Convert to microeco format and subset to only Exp2 water samples
sample_data <- phyloseq::sample_data(sample_metadata)
dat <- phyloseq::merge_phyloseq(genus_tab, tax_tab, sample_data) %>%
  phyloseq::subset_samples(experiment == "2") %>%
  phyloseq::subset_samples(sample_type == "Biofilm") %>%
  phyloseq::subset_taxa(Domain %in% c("Bacteria", "Archaea"))
dat_meco <- dat %>% file2meco::phyloseq2meco()
#Perform LefSE analysis
t1 <- microeco::trans_diff$new(dataset = dat_meco, method = "lefse", group = "sample_group", p_adjust_method = "fdr", lefse_subgroup = "sample_perc", lefse_min_subsam = 4)
df3 <- t1$res_diff %>%
  dplyr::mutate(phylum = word(word(Taxa, sep = "p__", 2), sep = "\\|", 1)) 

#Create cladogram
use_labels <- t1$res_abund %>%
  janitor::clean_names() %>%
  dplyr::filter(stringr::str_count(taxa, "__") == 2) %>%
  dplyr::arrange(desc(mean)) %>%
  dplyr::select(taxa) %>%
  dplyr::distinct() %>%
  dplyr::mutate(taxa = stringr::str_remove(taxa, "d__[:alpha_:]+\\|")) %>%
  unlist()
clades <- t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5, group_order = c("0% WW", "30-80 % WW", "30-80 % WW UF")) +
  geom_cladelab(node = 127, label = "p__Cyanobacteria", size = 2)
phyla_to_keep <- unique(clades$data$label[grepl("p__", clades$data$label)])
clades

ggsave("output/figures/supplemental_LefSE.jpg",
       clades,
       dpi = 300,
       device = "jpeg",
       units = "cm",
       width = 50,
       height = 75/1.6)

taxa <- t1$abund_table %>%
  rownames()
taxa
taxa <- tibble(taxpath = taxa) %>%
  dplyr::mutate(level = stringr::str_count(taxpath, "\\|") + 1,
                taxon = stringr::word(taxpath, sep = "\\|", start = level, end = level),
                phylum = dplyr::case_when(level == 1 ~ NA,
                                          TRUE ~ stringr::word(taxpath, sep = "\\|", start = 2, end = 2) %>%
                                            stringr::str_remove(pattern = "p__")))
clades_dat <- clades$data %>%
  dplyr::left_join(taxa, by = c("node_label" = "taxon"))
annot <- t1$res_diff %>% 
  janitor::clean_names() %>%
  dplyr::mutate(label = strsplit(taxa, split = "|", fixed = TRUE) %>% 
                  purrr::map_chr(utils::tail, n = 1))
annot <- clades_dat %>%
  dplyr::left_join(annot, by = c("node_label" = "label")) %>%
  dplyr::mutate(group = replace_na(group, "none")) #%>%
  #dplyr::filter(group != "none")

annot_highlight <- annot
  
# Set custom color palette
n_colors <- clades_dat$phylum %>%
  unique() %>%
  length()
pal <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)
pal2 <- gsub("#E41A1C", "gray40", pal)
pal3 <- gsub("#629363", "#E41A1C", pal2)
pal4 <- gsub("#F2E631", "royalblue3", pal3)
pal5 <- gsub("#BD6253", "black", pal4)
pal6 <- gsub("#BF862B", "lightslateblue", pal5)
pal7 <- gsub("#C4625D", "forestgreen", pal6)
pal8 <- gsub("#46A169", "navyblue", pal7)
pal9 <- gsub("#CE8BAE", "yellow4", pal8)
pal10 <- gsub("#815375", "deeppink", pal9)

ww80.roots <- c(200, 150, 145, 146, 160, 164, 191, 195, 199, 177)
ww80.colors <- rep("#CA9A81", length(ww80.roots))
ww30.roots <- c(95, 98, 97, 104, 106, 102, 103, 95, 187, 188, 100, 
                103, 108, 109, 110, 112, 124, 125, 126,  136, 137, 
                138, 139, 140, 141,101, 105)
ww30.colors <- rep("#EFEAB7",length(ww30.roots))
ww00.roots <- c(189, 180, 190)
ww00.colors <- rep("#93D5E7", length(ww00.roots))
group.roots <- c(ww80.roots, ww30.roots, ww00.roots)
group.colors <- c(ww80.colors, ww30.colors, ww00.colors)  
intersect(ww80.roots, ww30.roots)

# Manually fix taxa names to conform with Oren & Garrity nomenclature (2021)
g <- data.frame(gnode = group.roots, gfill = group.colors)
write_csv(data.frame(bind_cols(unique(clades_dat$phylum), pal10)), file = "data/phylum_palette.csv")

clades_dat$phylum <- gsub("Firmicutes", "Bacillota", clades_dat$phylum)
annot$phylum <- gsub("Firmicutes", "Bacillota", annot$phylum)
clades_dat$phylum <- gsub("Proteobacteria", "Pseudomonadota", clades_dat$phylum)
annot$phylum <- gsub("Proteobacteria", "Pseudomonadota", annot$phylum)
clades_dat$phylum <- gsub("Chloroflexi", "Chloroflexota", clades_dat$phylum)
annot$phylum <- gsub("Chloroflexi", "Chloroflexota", annot$phylum)
clades_dat$phylum <- gsub("Cyanobacteria", "Cyanobacteriota", clades_dat$phylum)
annot$phylum <- gsub("Cyanobacteria", "Cyanobacteriota", annot$phylum)

pdf("output/figures/figure1_LEfSe_with_legend.pdf", width = 6, height = 6)
ggtree::ggtree(clades_dat, layout = "circular", aes(color = phylum)) +
  ggtree::geom_point2(data = annot, aes(size = I(node_size), fill = group), shape = 21) +
  ggplot2::labs(fill = "Enrichment", color = "Phylum") +
  ggplot2::scale_fill_manual(values = c("#93D5E7", "#EFEAB7", "#CA9A81", "white"), guide = guide_legend(label.theme = element_text(size = 18))) +
  ggplot2::scale_color_manual(values = pal10,
                              guide = guide_legend(label.theme = element_text(face = "italic", size = 14))) +
  ggplot2::theme(legend.title = element_text(size = 16)) +
                # legend.key.size = unit(3, "cm")) +
                # legend.key.width = unit(1,"cm"),
                # legend.key.height= unit(1, 'cm')) +
  lapply(seq(nrow(g)), function(i) {
    geom_hilight(
      node = g$gnode[i],
      alpha = 0.5,
      type = "roundrect",
      fill = g$gfill[i],
      to.bottom = TRUE
    )
  })
dev.off()

ggsave("output/figures/figure1_LEfSe_legend.jpg",
       dpi = 300,
       device = "jpeg",
       units = "cm",
       width = 14,
       height = 14)

       