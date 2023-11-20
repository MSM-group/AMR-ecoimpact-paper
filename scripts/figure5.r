#clear_environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape", "edgeR", "DESeq2")
#import mtags genus level data
mtags_in <- readr::read_tsv("data/mtags/all.family.tsv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(family_id = paste0("family_", as.character(dplyr::row_number())), .before = number_taxpath)

#create genus table
mtags_families <- mtags_in %>%
  dplyr::select(-number_taxpath) %>%
  tidyr::pivot_longer(-family_id, names_to = "sample", values_to = "count") %>%
  dplyr::mutate(sample = stringr::str_remove(sample, pattern = "_bins") %>%
                  stringr::str_to_upper()) %>%
  tidyr::pivot_wider(names_from = "family_id", values_from = "count")  %>%
  tibble::column_to_rownames(var = "sample")

#create taxonomy tabel
get_family <- function(tp){
  case_when(tp %in% c("Unassigned", "Unaligned") ~ tp,
            TRUE ~ stringr::word(tp, sep = ";", start = 3, end = 3))
}
modify_names <- function(nm){
  dplyr::case_when(str_detect(nm, "__") ~ stringr::word(nm, sep = "__", start = 6, end = 6),
                   TRUE ~ nm)
}
mtags_tax <- mtags_in %>%
  dplyr::select(family_id, number_taxpath) %>%
  dplyr::rename(taxpath = number_taxpath) %>%
  dplyr::mutate(root = stringr::word(taxpath, sep = "__", 2),
                domain = stringr::word(taxpath, sep = "__", 3),
                phylum = stringr::word(taxpath, sep = "__", 4),
                class = stringr::word(taxpath, sep = "__", 5),
                order = stringr::word(taxpath, sep = "__", 6),
                family = stringr::word(taxpath, sep = "__", 7))

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
  dplyr::filter(sample_code %in% rownames(mtags_families))

#remove unaligned reads if desired
remove_unaligned <- FALSE
if(remove_unaligned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(family != "Unaligned")
  mtags_families <- mtags_families %>%
    dplyr::select(dplyr::all_of(mtags_tax$family_id))
}
#remove unassigned reads if desired
remove_unassigned <- FALSE
if(remove_unassigned == TRUE){
  mtags_tax <- mtags_tax %>%
    dplyr::filter(family != "Unassigned")
  mtags_families <- mtags_families %>%
    dplyr::select(dplyr::all_of(mtags_tax$family_id))
}
#bring into DESeq2 format
DESeq2_coldata <- sample_metadata %>%
  dplyr::arrange(sample_code) %>%
  filter(stringr::str_detect(sample_perc, "WW30")) %>%
  tibble::column_to_rownames(var = "sample_code") %>%
  as.matrix()

DESeq2_counts <- mtags_families %>%
  dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
  as.matrix() %>%
  t() 
DESeq2_counts <- DESeq2_counts[,colnames(DESeq2_counts) %in% rownames(DESeq2_coldata)]
colnames(DESeq2_counts) == rownames(DESeq2_coldata)
#Fit DESeq2 model with treatment as explanatory model
dds <- DESeq2::DESeqDataSetFromMatrix(countData = DESeq2_counts,
                                      colData = DESeq2_coldata,
                                      design = ~ sample_perc)
#Remove families w/ <10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#Relevel factors to set 30% WW UF as reference
dds$sample_perc <- factor(dds$sample_perc, levels = c("WW30", "WW30UF"))
dds$sample_perc
#Run differential abundance analysis
dds <- DESeq(dds)
resu <- results(dds, contrast = c("sample_perc", "WW30UF", "WW30"))
resu
#convert results to tibble format for further analysis and add taxonomic data
res <- resu %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "family_id") %>%
  tibble::as_tibble() %>%
  dplyr::left_join(mtags_tax, by = "family_id")
res_sig <- res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)
#Make volcano plot
is_sig <- function(log2fc, p){
  dplyr::case_when(log2fc > 2 & p < 0.05 ~ "up",
                   log2fc < -2 & p < 0.05 ~ "down",
                   TRUE ~ "no")
}
make_label <- function(phyl){
  dplyr::case_when(grepl("Plancto", class) ~ class,
                   TRUE ~ "")
}

res <- res %>%
  drop_na(padj) %>%
  dplyr::mutate(minusLog10padj = -log10(padj),
                significance = is_sig(log2FoldChange, padj),
                lab = make_label(phylum)) 

sigdat <- res %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(grepl("Bacteria", domain)) %>%
  dplyr::filter(abs(log2FoldChange) > 2) %>%
  dplyr::filter(abs(minusLog10padj) > 10)

volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = minusLog10padj)) +
  ggplot2::geom_point(aes(color = significance), alpha = 0.5, stroke = 0) +
  ggplot2::geom_text(data = subset(res, log2FoldChange > 2.4), aes(label = gsub(";class", "", lab), hjust = 0, vjust = 0.5, fontface = "italic", size = 3)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  ggplot2::geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  ggplot2::geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  ggplot2::scale_color_manual(values = c("red", "gray", "blue")) +
  ggplot2::xlim(-10, 10) +
  ggpubr::theme_pubr() +
  ggplot2::theme(legend.position = "none") +
  labs(x = expression(log[2](FC)), y = expression(-log[10](p*"-"*adj)))
volcano_plot

ggsave("output/figures/figure5a.jpg",
       volcano_plot,
       device = "jpeg",
       dpi = 300,
       units = "cm",
       width = 14,
       height = 10)
