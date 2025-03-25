# Load required libraries
library(TreeSummarizedExperiment)
library(S4Vectors)
library(mia)

data("data_survival", package = "coda4microbiome")
data_survival <- x
assay_matrix <- as.matrix(data_survival)
assay_matrix <- t(assay_matrix)
taxa_names <- rownames(assay_matrix)

# Function to Parse Taxonomy Labels
parse_taxonomy <- function(taxa_names) {
  taxonomy_df <- data.frame(TaxonID = taxa_names, stringsAsFactors = FALSE)
  
  # Add structured taxonomy placeholders
  taxonomy_df$Kingdom <- NA  
  taxonomy_df$Phylum <- NA
  taxonomy_df$Class <- NA
  taxonomy_df$Order <- NA
  taxonomy_df$Family <- NA
  taxonomy_df$Genus <- NA
  
  # Fill known taxonomic ranks
  for (i in seq_along(taxa_names)) {
    taxon <- taxa_names[i]
    if (grepl("^g__", taxon)) {
      taxonomy_df$Genus[i] <- sub("^g__", "", taxon)
    } else if (grepl("^f__", taxon)) {
      taxonomy_df$Family[i] <- sub("^f__", "", taxon)
    } else if (grepl("^o__", taxon)) {
      taxonomy_df$Order[i] <- sub("^o__", "", taxon)
    }
  }
  
  return(taxonomy_df)
}

Event <- as.numeric(Event)  
Event_time <- as.numeric(Event_time) 

rowData <- DataFrame(parse_taxonomy(taxa_names), row.names = taxa_names)
colData <- DataFrame(SampleID = colnames(assay_matrix), row.names = colnames(assay_matrix), Event, Event_time)

data_survival <- TreeSummarizedExperiment(
  assays = list(counts = assay_matrix),
  rowData = rowData,
  colData = colData
)

# Convert logical columns to character
owData(data_survival)$Kingdom <- as.character(rowData(data_survival)$Kingdom)
rowData(data_survival)$Phylum <- as.character(rowData(data_survival)$Phylum)
rowData(data_survival)$Class <- as.character(rowData(data_survival)$Class)

save(data_survival, file = "~/miaTime/data/data_survival.rda")
