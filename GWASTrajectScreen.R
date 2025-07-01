# set of functions to merge pseudotime analysis with GWAS screens 



library(tidyverse)
library(stringr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggrepel)
library(biomaRt)
library(ggplot2)
library(data.table)



run_gwas_overlap <- function(gene_list, gwas_df, padding = 250000, mart) {
  
  gene_map <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                    filters = "hgnc_symbol",
                    values = gene_list,
                    mart = mart)
  
  if (nrow(gene_map) == 0) {
    message("No gene mappings found.")
    return(NULL)
  }
  
  exons <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", 
                                "exon_chrom_start", "exon_chrom_end"),
                 filters = "ensembl_gene_id",
                 values = gene_map$ensembl_gene_id,
                 mart = mart)
  
  gene_coords <- exons %>%
    filter(chromosome_name %in% as.character(1:22)) %>%
    mutate(chromosome = as.integer(chromosome_name)) %>%
    left_join(gene_map, by = "ensembl_gene_id") %>%
    group_by(hgnc_symbol, chromosome) %>%
    summarise(start = min(exon_chrom_start, exon_chrom_end) - padding,
              end = max(exon_chrom_start, exon_chrom_end) + padding,
              .groups = "drop")
  
  


    # 6. Query SNPs near those genes
    deg_hits <- map_dfr(seq_len(nrow(gene_coords)), function(i) {
      row <- gene_coords[i, ]
      gene <- row$hgnc_symbol
      chr <- row$chromosome
      start <- row$start
      end <- row$end
    
    
    
      snps_in_gene <- gwas %>%
        filter(chromosome == chr,
               base_pair_location >= start,
               base_pair_location <= end)
      
      if (nrow(snps_in_gene) == 0) return(NULL)
      
      snps_in_gene %>%
        mutate(matched_gene = gene)
    })
    
 
    if (nrow(deg_hits) > 0) {
      deg_hits <- as.data.table(deg_hits)
      deg_hits[, fdr := p.adjust(p_value, method = "fdr")]
      significant_hits <- deg_hits[fdr < 0.05, .(chromosome, base_pair_location, matched_gene, p_value, fdr, beta, standard_error)]
      
      message("This run: ", length(unique(significant_hits$matched_gene)), 
              " unique matched genes with significant SNPs")
      return(significant_hits)
    } else {
      message("No significant hits.")
      return(NULL)
    }
    
}


DAMvsHomeo.markers <- read_csv("/Volumes/PC60/PDGWASStuff/DamvsHomeMarkers.csv", show_col_types = FALSE)






gwas <- fread("/Volumes/PC60/PDGWASStuff/GCST009325.tsv")
head(gwas)




# 1. Set up Ensembl grch38 connection
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")




# DEG gene names
geneNames=DAMvsHomeo.markers$...1
deg_genes <- (geneNames[DAMvsHomeo.markers$p_val_adj < 0.05 ])


# Running actual overlap analysis
UKBioHits <- run_gwas_overlap(deg_genes, gwas, mart = mart)

cat("Real: ", length(unique(UKBioHits$matched_gene)), "unique matched genes with significant SNPs\n")


write.csv(significant_hits, "/Volumes/PC60/PDGWASStuff/GWASCompleteMatchUKBioBankData500KMBWindow_Final.csv")





#settin seeds
set.seed(42)

# Pre-allocate vector
perm_results <- numeric(50)
perm_length <- numeric(50)

# Loop to generate permutations, this takes like 3 hours
for (i in 1:50) {
  cat("Running permutation", i, "...\n")
  
  random_genes <- sample(setdiff(geneNames, deg_genes), size = length(deg_genes))
  result <- run_gwas_overlap(random_genes, gwas, mart = mart)
  
  if (is.null(result)) {
    perm_results[i] <- 0
  } else {
    perm_results[i] <- length(unique(result$matched_gene))
    perm_length[i]<- length(result$matched_gene)
  }
  
}



## now making manhattan plot

library(ggrepel)

fdr_logp_cutoff <- -log10(0.05)i


# Step 1: One SNP per unique gene (with smallest p-value)
label_points <- significant_hits %>%
  group_by(matched_gene) %>%
  slice_min(order_by = p_value, n = 1) %>%
  ungroup() %>%
  mutate(logp = -log10(p_value))


# Start with your DEG-associated GWAS hits (assumes this came from earlier filtering)
gwas_for_plot <- deg_hits %>%
  filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value)) %>%
  transmute(
    CHR = as.numeric(chromosome),
    BP = as.numeric(base_pair_location),
    P = as.numeric(p_value),
    matched_gene = matched_gene,
    beta = beta,
    sebeta = standard_error,
    logp = -log10(P)
  )



# Step 2: Merge with plotting data to get CHR and BP
label_points_plot <- label_points %>%
  rename(BP = base_pair_location, CHR = `chromosome`) %>%
  mutate(CHR = as.numeric(CHR))

# Step 3: Plot with labels
ggplot(gwas_for_plot, aes(x = BP, y = logp)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = fdr_logp_cutoff, linetype = "dashed", color = "red") +
  geom_text_repel(
    data = label_points_plot,
    aes(label = matched_gene),
    size = 2.2,
    min.segment.length = 0,
    box.padding = 0.25,
    max.overlaps = 50
  ) +
  facet_wrap(~ CHR, scales = "free_x", nrow = 1, strip.position = "bottom") +
  labs(x = "Chromosome", y = "-log10(p-value)", title = "Manhattan Plot of NADEL et al GWAS filtered for DAM Trajectory Genes") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0.1, "lines")
  )
