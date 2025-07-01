

#First just testing what known GWAS hits are sig different in microglia cluster

DAMvsHomeo.markers <- read_csv("/Volumes/PC60/PDGWASStuff/DamvsHomeMarkers.csv", show_col_types = FALSE)


gwas_data <- read_tsv("/Volumes/PC60/PDGWASStuff/GWASPDAssociation.tsv", show_col_types = FALSE)


# Get DEG table with rownames as a column
deg_df <- DAMvsHomeo.markers  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column("gene")

# Function to split and trim gene names
split_and_trim <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(str_split(x, "[-, ]+")) %>% str_trim() %>% unique()
}

# Match GWAS entries to DEG list and extract padj
gwas_matches <- gwas_data %>%
  rowwise() %>%
  mutate(
    mapped_genes = list(split_and_trim(MAPPED_GENE)),
    reported_genes = list(split_and_trim(`REPORTED GENE(S)`)),
    matched_genes_mapped = list(intersect(mapped_genes, deg_df$gene)),
    matched_genes_reported = list(intersect(reported_genes, deg_df$gene)),
    matched_any = length(matched_genes_mapped) > 0 | length(matched_genes_reported) > 0
  ) %>%
  ungroup() %>%
  filter(matched_any) %>%
  mutate(
    matched_gene_names = map2_chr(matched_genes_mapped, matched_genes_reported, ~paste(unique(c(.x, .y)), collapse = ", ")),
    MicroglialTrajectoryPvalAdj = map2_chr(matched_genes_mapped, matched_genes_reported, ~{
      matched <- unique(c(.x, .y))
      padj_vals <- deg_df %>% filter(gene %in% matched) %>% pull(p_val_adj)
      paste(format(padj_vals, digits = 3, scientific = TRUE), collapse = ", ")
    })
  ) %>%
  select(SNPS, CHR_ID, CHR_POS, MAPPED_GENE, `REPORTED GENE(S)`,
         matched_gene_names, MicroglialTrajectoryPvalAdj,
         INTERGENIC, STUDY, `DISEASE/TRAIT`, `FIRST AUTHOR`, `DATE ADDED TO CATALOG`, JOURNAL)

# View result
print(gwas_matches)

write_csv(gwas_matches, "/Volumes/PC60/PDGWASStuff/Transition_vs_Homeostatic_Matches_allPDGWAS.csv")
#######################
