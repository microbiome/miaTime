# Load required libraries
library(TreeSummarizedExperiment)
library(S4Vectors)

# Get data from coda4microbiome package
data("data_survival", package = "coda4microbiome")
crohn_survival <- x
assay_matrix <- as.matrix(crohn_survival)
assay_matrix <- t(assay_matrix)
taxa_names <- rownames(assay_matrix)

# Prepare rowData
rd <- lapply(mia:::.taxonomy_rank_prefixes, function(rank){
    tax <- rep(NA_character_, length(taxa_names))
    prefix <- paste0("^", rank, "__")
    found <- grepl(prefix, taxa_names)
    tax[found] <- gsub(prefix, "", taxa_names[found])
    return(tax)
})
rd <- do.call(cbind, rd) |> DataFrame()
empty <- sapply(rd, function(x) all(is.na(x)))
rd <- rd[, !empty]
rd[["taxonomy_id"]] <- taxa_names

# Prepare colData
event <- as.numeric(Event)
event_time <- as.numeric(Event_time)
diagnosis <- ifelse(event == 1, "CD", "control")

cd <- DataFrame(
    sample_id = colnames(assay_matrix),
    event = as.integer(event),
    event_time,
    diagnosis,
    row.names = colnames(assay_matrix)
)

# Construct TreeSE
crohn_survival <- TreeSummarizedExperiment(
     assays = SimpleList(counts = assay_matrix),
     rowData = rd,
     colData = cd
   )

save(crohn_survival, file = "data/crohn_survival.rda")
