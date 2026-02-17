# Build GMT database for cotton
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)

# Load G.hirsutum annotation file
gh_annot <- read.delim(
  "data/gh_annot/Ghirsutum_578_v3.1.annotation_info.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = ""
)

head(gh_annot)

# Extract the ortholog mapping
orthologs <- gh_annot %>%
  dplyr::select(
    gh_gene_id = locusName,
    at_besthit = Best.hit.arabi.name
  ) %>%
  dplyr::filter(!is.na(at_besthit), at_besthit != "")

# Convert gmt_at into long format (gmt_at object generated from arabidopsis_GMT.R script)
pathway_genes_at <- gmt_at %>%
  mutate(gene_list = str_split(genes, "\t")) %>%
  unnest(gene_list) %>%
  rename(gene_id = gene_list) %>%
  dplyr::select(pathway_id, pathway_name, gene_id)

head(pathway_genes_at)

# Join G.hirsutum with AGI
pathway_genes_gh <- pathway_genes_at %>%
  inner_join(orthologs, by = c("gene_id" = "at_besthit")) %>%
  dplyr::select(pathway_id, pathway_name, gh_gene_id)

# Build the G.hirsutum GMT
gmt_gh <- pathway_genes_gh %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(genes = paste(unique(gh_gene_id), collapse = "\t"), .groups = "drop") %>%
  mutate(gmt_line = paste(pathway_id, pathway_name, genes, sep = "\t"))

writeLines(gmt_gh$gmt_line, "output/gossypium_hirsutum_pathways.gmt")

# inspect
library(GSEABase)

gmt_gh_obj <- getGmt("output/gossypium_hirsutum_pathways.gmt")
length(gmt_gh_obj)
gmt_gh_obj[[1]]
geneIds(gmt_gh_obj[[1]])

# How to verify which pathways were lost
at_ids  <- sapply(gmt_obj, setName)
gh_ids  <- sapply(gmt_gh_obj, setName)

lost <- setdiff(at_ids, gh_ids)
lost



