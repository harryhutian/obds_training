library(tidyverse)
library(biomaRt)
library(pheatmap)

sample_table <- read_tsv("data/tidyverse/obds_sampletable.tsv")
count_table <- read_tsv("data/tidyverse/obds_countstable.tsv.gz")

processed_countable <- count_table %>%
  pivot_longer(-Geneid, names_to = "sample", values_to = "count")

# Join with gene info to get mgi_symbol
listMarts()
ensembl <- useMart("ensembl") # connecting to specific database
datasets <- listDatasets(ensembl)  # list dataset
head(datasets)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listFilters(ensembl)
attributes[1:5,]


# matching
ensembl_gene_id
mgi_symbol

#packageVersion('biomaRt')

matching <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
      values = unique(processed_countable$Geneid),
      mart = ensembl
      )
# joining tables to match geneid
processed_countable <- processed_countable %>%
  left_join(matching, by = c("Geneid" = "ensembl_gene_id"))

# tidy metadata file:
# seperate
processed_sample_table <- sample_table %>%
  separate(sample_title, c("genotype", "knockout","cell type","replicates"), sep = "_") %>%
  unite("Genotype", genotype,knockout, sep = "_" )%>%
  dplyr::select(-library_layout, -read_count)

# joining 2 tables
processed_joined <- processed_countable %>%
  left_join(processed_sample_table, by = c("sample" = "Sample_accession"))

# add a new column with count per million
calculated <- processed_joined %>%
  group_by(sample) %>%
  mutate(total_count_per_sample = sum(count)) %>%
  mutate(total_count_in_million = total_count_per_sample/1000000) %>%
  mutate(CPM = count/total_count_in_million) %>%
  mutate(log_CPM = log2(CPM+1))
