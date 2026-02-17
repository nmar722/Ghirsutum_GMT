# Built GMT database for Arabidopsis

library(tidyverse)
library(data.table)

# Inspect the data
readLines("data/aracyc/genes.dat", n = 50)

# Define read_dat_blocks() function
read_dat_blocks <- function(file) {
  lines <- readLines(file)
  # split into blocks using the "//" separator
  blocks <- split(lines, cumsum(lines == "//"))
  # remove empty blocks (sometimes the last one)
  blocks <- blocks[sapply(blocks, function(x) any(x != "//"))]
  blocks
}

# 1. Parse genes.dat
parse_aracyc_genes <- function(file) {
  blocks <- read_dat_blocks(file)
  
  tibble(
    tair_id = map_chr(blocks, ~{
      line <- .x[grepl("^UNIQUE-ID - ", .x)]
      if (length(line) == 0) NA else sub("^UNIQUE-ID - ", "", line)
    }),
    symbol = map_chr(blocks, ~{
      line <- .x[grepl("^COMMON-NAME - ", .x)]
      if (length(line) == 0) NA else sub("^COMMON-NAME - ", "", line)
    }),
    products = map(blocks, ~{
      lines <- .x[grepl("^PRODUCT - ", .x)]
      sub("^PRODUCT - ", "", lines)
    })
  )
}

genes <- parse_aracyc_genes("data/aracyc/genes.dat")
head(genes)

# 2. Parse proteins.dat
parse_aracyc_proteins <- function(file) {
  blocks <- read_dat_blocks(file)
  
  tibble(
    protein_id = map_chr(blocks, ~{
      line <- .x[grepl("^UNIQUE-ID - ", .x)]
      if (length(line) == 0) return(NA)
      sub("^UNIQUE-ID - ", "", line[1])   # take FIRST match
    }),
    
    gene_id = map_chr(blocks, ~{
      line <- .x[grepl("^GENE - ", .x)]
      if (length(line) == 0) return(NA)
      sub("^GENE - ", "", line[1])        # take FIRST match
    }),
    
    species = map_chr(blocks, ~{
      line <- .x[grepl("^SPECIES - ", .x)]
      if (length(line) == 0) return(NA)
      sub("^SPECIES - ", "", line[1])     # take FIRST match
    })
  )
}

proteins <- parse_aracyc_proteins("data/aracyc/proteins.dat")
head(proteins)

# 3. Parse pathways
parse_aracyc_pathways <- function(file) {
  blocks <- read_dat_blocks(file)
  
  tibble(
    pathway_id = map_chr(blocks, ~{
      line <- .x[grepl("^UNIQUE-ID - ", .x)]
      if (length(line) == 0) return(NA_character_)
      sub("^UNIQUE-ID - ", "", line[1])
    }),
    
    pathway_name = map_chr(blocks, ~{
      line <- .x[grepl("^COMMON-NAME - ", .x)]
      if (length(line) == 0) return(NA_character_)
      sub("^COMMON-NAME - ", "", line[1])
    }),
    
    reactions = map(blocks, ~{
      lines <- .x[grepl("^REACTION-LIST - ", .x)]
      if (length(lines) == 0) return(character(0))
      sub("^REACTION-LIST - ", "", lines)
    })
  )
}

pathways <- parse_aracyc_pathways("data/aracyc/pathways.dat")
head(pathways)

# 4. Parse enxrxns.dat
parse_aracyc_enzrxns <- function(file) {
  blocks <- read_dat_blocks(file)
  
  tibble(
    enzrxn_id = map_chr(blocks, ~{
      line <- .x[grepl("^UNIQUE-ID - ", .x)]
      if (length(line) == 0) return(NA_character_)
      sub("^UNIQUE-ID - ", "", line[1])
    }),
    
    reaction_id = map_chr(blocks, ~{
      line <- .x[grepl("^REACTION - ", .x)]
      if (length(line) == 0) return(NA_character_)
      sub("^REACTION - ", "", line[1])
    }),
    
    enzymes = map(blocks, ~{
      lines <- .x[grepl("^ENZYME - ", .x)]
      if (length(lines) == 0) return(character(0))
      sub("^ENZYME - ", "", lines)
    })
  )
}

enzrxns <- parse_aracyc_enzrxns("data/aracyc/enzrxns.dat")
head(enzrxns)

# 5. Parse reactions.dat
parse_aracyc_reactions <- function(file) {
  blocks <- read_dat_blocks(file)
  
  tibble(
    reaction_id = map_chr(blocks, ~{
      line <- .x[grepl("^UNIQUE-ID - ", .x)]
      if (length(line) == 0) return(NA_character_)
      sub("^UNIQUE-ID - ", "", line[1])
    }),
    
    enzrxn_ids = map(blocks, ~{
      lines <- .x[grepl("^ENZYMATIC-REACTION - ", .x)]
      if (length(lines) == 0) return(character(0))
      sub("^ENZYMATIC-REACTION - ", "", lines)
    })
  )
}

reactions <- parse_aracyc_reactions("data/aracyc/reactions.dat")
head(reactions)

# Join pipeline
## 1. Pathway → Reaction
pw_rxn <- pathways %>%
  select(pathway_id, pathway_name, reactions) %>%
  unnest(reactions) %>%
  rename(reaction_id = reactions)

## 2. Reaction → Enzymatic-Reaction
rxn_to_enzrxn <- reactions %>%
  unnest(enzrxn_ids) %>%
  rename(enzrxn_id = enzrxn_ids)

## 3. Clean enzrxns to avoid duplicate reaction_id
enzrxns_clean <- enzrxns %>%
  select(enzrxn_id, enzymes)

## 4. Enzymatic-Reaction → Enzyme (protein)
rxn_to_protein <- rxn_to_enzrxn %>%
  left_join(enzrxns_clean, by = "enzrxn_id") %>%
  unnest(enzymes) %>%
  rename(protein_id = enzymes)

## 5. Protein → Gene
rxn_to_gene <- rxn_to_protein %>%
  left_join(proteins, by = "protein_id") %>%
  select(reaction_id, gene_id)

## 6. Pathway → Gene
pathway_genes <- pw_rxn %>%
  left_join(rxn_to_gene, by = "reaction_id") %>%
  select(pathway_id, pathway_name, gene_id) %>%
  distinct()

head(pathway_genes)

# ------------ Building Arabidopsis GMT file
# 1. Inspect pathway_genes
# Count of NA
sum(is.na(pathway_genes$gene_id))

# Count NA per pathway
pathway_genes %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(n_NA = sum(is.na(gene_id)), .groups = "drop") %>%
  filter(n_NA > 0)

# Count pathways that become empty after removing NA
pathway_genes %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(n_genes = sum(!is.na(gene_id)), .groups = "drop") %>%
  filter(n_genes == 0)

# Remove NA rows
pathway_genes_clean <- pathway_genes %>% 
  filter(!is.na(gene_id))
# Keep only Arabidopsis real names
pathway_genes_at <- pathway_genes_clean %>%
  filter(str_detect(gene_id, "^AT[1-5CM]G\\d{5}$"))

sum(grepl("^AT[1-5CM]G[0-9]{5}$", pathway_genes$gene_id))

# Identify and document empty pathways
all_pathways <- pathway_genes %>%
  distinct(pathway_id, pathway_name) #pathway originally present

nonempty_pathways <- pathway_genes_at %>%
  distinct(pathway_id)               #survived AGI filtering

empty_after_filter <- all_pathways %>%
  anti_join(nonempty_pathways, by = "pathway_id") #pathway that become empty

write.csv(empty_after_filter, "output/empty_pathways.csv",row.names = FALSE)

## Define valid pathways
valid_pathways <- pathway_genes_at

write.csv(valid_pathways, 
          "output/pathways_genes_arabidopsis_valid.csv", row.names = FALSE)
head(valid_pathways)

# Clean valid_pathways HTML formatting
# Remove HTML tags like <i>...</i>
library(stringr)

valid_pathways_clean <- valid_pathways %>%
  mutate(pathway_name = str_replace_all(pathway_name, "<[^>]+>", ""))

# Decode HTML entities (&alpha;, &beta;, etc.)
library(xml2)

decode_html <- function(x) {
  xml2::xml_text(xml2::read_html(paste0("<x>", x, "</x>")))
}

valid_pathways_clean <- valid_pathways_clean %>%
  mutate(pathway_name = vapply(pathway_name, decode_html, FUN.VALUE = character(1)))

head(valid_pathways_clean$pathway_name)

# ----------- Build the GMT file
library(dplyr)
library(stringr)

gmt_at <- valid_pathways_clean %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(genes = paste(unique(gene_id), collapse = "\t"), .groups = "drop") %>%
  mutate(gmt_line = paste(pathway_id, pathway_name, genes, sep = "\t"))

# Save GMT file
write.csv(gmt_at, "output/gmt_at.csv", row.names = FALSE)
writeLines(gmt_at$gmt_line, "output/arabidopsis_pathways.gmt")

#-------------------
# Install GSEABase
#BiocManager::install("GSEABase", ask = FALSE)
library(GSEABase)

gmt_obj <- getGmt("output/arabidopsis_pathways.gmt")
length(gmt_obj)
gmt_obj[[1]]
geneIds(gmt_obj[[1]])

library(tibble)
library(purrr)

gmt_df <- tibble(
  pathway_id   = names(gmt_obj),
  pathway_name = map_chr(gmt_obj, setName),
  genes        = map(gmt_obj, geneIds)
)





