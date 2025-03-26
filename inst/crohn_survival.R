# Load required libraries
library(TreeSummarizedExperiment)
library(S4Vectors)


data("data_survival", package = "coda4microbiome")
crohn_survival <- x
assay_matrix <- as.matrix(crohn_survival)
assay_matrix <- t(assay_matrix)
taxa_names <- rownames(assay_matrix)
str(x)
# Function to Parse Taxonomy Labels
parse_taxonomy <- function(taxa_names) {
  # Create the taxonomy data frame
  taxonomy_df <- data.frame(TaxonID = taxa_names, stringsAsFactors = FALSE)
  
  # Add structured taxonomy placeholders
  taxonomy_df$Kingdom <- NA_character_
  taxonomy_df$Phylum <- NA_character_
  taxonomy_df$Class <- NA_character_
  taxonomy_df$Order <- NA_character_
  taxonomy_df$Family <- NA_character_
  taxonomy_df$Genus <- NA_character_
  taxonomy_df$Species <- NA_character_
  
  # Return the modified taxonomy dataframe
  return(taxonomy_df)
}


# Fill known taxonomic ranks


# Create an example taxonomy dataframe with columns for Genus, Family, and Order
taxonomy_df <- data.frame(Genus = rep(NA, length(taxa_names)),
                          Family = rep(NA, length(taxa_names)),
                          Order = rep(NA, length(taxa_names)),
                          stringsAsFactors = FALSE)

# Now your loop should work as expected
taxonomy_df$Genus[grepl("^g__", taxa_names)] <- sub("^g__", "", taxa_names[grepl("^g__", taxa_names)])
taxonomy_df$Family[grepl("^f__", taxa_names)] <- sub("^f__", "", taxa_names[grepl("^f__", taxa_names)])

Event <- as.numeric(Event)  
Event_time <- as.numeric(Event_time) 

rowData <- DataFrame(parse_taxonomy(taxa_names), row.names = taxa_names)
colData <- DataFrame(SampleID = colnames(assay_matrix), Event, Event_time, row.names = colnames(assay_matrix))

crohn_survival <- TreeSummarizedExperiment(
     assays = SimpleList(counts = assay_matrix),
     rowData = rowData,
     colData = colData
   )


save(crohn_survival, file = "data/crohn_survival.rda")







