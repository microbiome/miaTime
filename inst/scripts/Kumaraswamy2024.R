# Load required libraries
library(TreeSummarizedExperiment)
library(Cairo)
library(dplyr)
library(readxl)

# Function to load Excel data
read_data <- function (f) {
    x <- read_excel(f)
    rownams <- unname(unlist(x[,1]))
    x <- x[, -1]
    x <- as.matrix(x)
    rownames(x) <- rownams
    return(x)
}
# Get data from : https://zenodo.org/records/14424024
# Abundance profiles
gen <- read_data(file.path("..", "data", "Genus_hitchip.xlsx"))
phy <- read_data(file.path("..", "data", "Phylum_hitchip.xlsx"))
oli <- read_data(file.path("..", "data", "Oligo_hitchip.xlsx"))

# Metadata
md <- read_data(file.path("..", "data", "modified_file.xlsx"))
rownames(md) <- unname(md[, "sample"])
md <- as.data.frame(md)
md[14:61] <- lapply(md[14:61], as.logical)
# Group-A: never consumed Hawaijar and Dahi (n=20, control)
# Group-B: consume Hawaijar and Dahi (n=21)
# Group-C: consume Hawaijar, not Dahi (n=23)
# Group-D: consume Dahi, not Hawaijar (n=14)
md[, "timepoint"] <- as.numeric(unlist(md[, "timepoint"]))
md[, "season"] <- factor(unlist(md[, "season"]),
                         evels=c("summer", "autumn", "winter"))
factors <- c("age", "sex", "bmi", "clan", "nature_of_birth",
             "marital_status", "residence", "subject", "group")
for(f in factors) {
  md[, f] <- factor(unlist(md[, f]), levels=sort(unique(md[, f])))
}

# Create tse data object
tse <-TreeSummarizedExperiment(
    assays=SimpleList(signal=gen), colData=DataFrame(md))
# Add altExps
altExp(tse, "phylum")  <- TreeSummarizedExperiment(
    assays=SimpleList(signal=phy))

altExp(tse, "oligo")   <- TreeSummarizedExperiment(
    assays=SimpleList(signal=oli))
# There is one NA, replace it with min value
assay(altExp(tse, "oligo"), "signal")[is.na(assay(altExp(tse, "oligo"), "signal"))] <- min(assay(altExp(tse, "oligo"), "signal"), na.rm=TRUE)
# Round the oligo assay to 8 decimal places
assay(altExp(tse, "oligo")) <- round(assay(altExp(tse, "oligo")), 8)
# -------------------------------------------

# Total load in LOG10_16S _RNA_gene copies_per_g
# tabs 6 and 8 have different sample names
tabs <- list()
for (i in 1:11) {
  tabs[[i]] <- read_excel(file.path("..", "data", "AbsoluteloadTaxaspecificqPCRdata.xlsx"), sheet = i)
}
tabs <- tabs[-c(6,8)]
d <- Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="sample"), tabs)
d <- data.frame(d)
rownams <- unname(unlist(d[, "sample"]))
d <- d[, -1]
d[d %in% c("missing data", "NA")] <- NA
d <- apply(d, 2, as.numeric)
rownames(d) <- rownams
altExp(tse, "total_loads")  <- TreeSummarizedExperiment(assays=SimpleList(signal=t(d)))

# 'Fecal metabolite profile_LC-HRMS Data.xlsx'
x <- read_excel(file.path("..", "data", "Fecal\ metabolite\ profile_LC-HRMS Data.xlsx"), sheet = 1)
colnams <- as.character(x[3,])
x <- x[-c(1,2,3),]
colnames(x) <- colnams
xr <- x[, 1:5]
rownames(xr) <- paste0("feature_", 1:nrow(xr))
xd <- apply(as.matrix(x[, 6:ncol(x)]), 2, as.numeric)
M <- matrix(NA, nrow=nrow(xd), ncol=ncol(tse))
colnames(M) <- colnames(tse)
# Match samples
M[, colnames(xd)] <- xd
rownames(M) <- rownames(xr)
altExp(tse, "metabolites")  <- TreeSummarizedExperiment(
    assays=SimpleList(signal=M), rowData=xr)

# 'SCFA data-HPLC.xlsx'
x <- read_excel(file.path("..", "data", "SCFA\ data-HPLC.xlsx"))
colnams <- unname(unlist(x[1,]))
x <- x[-1, ]
colnames(x) <- colnams
rownams <- x$sample
x <- x[,-1]
x <- as.matrix(x)
x <- apply(x,2,as.numeric)
scfa <- t(x)
colnames(scfa) <- rownams
M <- matrix(NA, nrow=nrow(scfa), ncol=ncol(tse))
colnames(M) <- colnames(tse)
M[, colnames(scfa)] <- scfa
rownames(M) <- colnams[-1]
altExp(tse, "scfa")  <- TreeSummarizedExperiment(assays=SimpleList(signal=M))

# -----------------------------------------------------------------------------
save(tse, file = "Kumaraswamy2024.rda", compress = "xz")
