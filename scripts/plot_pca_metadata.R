#!/usr/bin/env Rscript

# Overlay metadata as coloured points onto PCA

library(tidyverse)
library(ggplot2)
library(glue)
library(googleCloudStorageR)
library(gargle)
library(magrittr)
library(viridis)

# Google cloud setup/token authorisation (to get files from GCP)
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- gargle::token_fetch(scopes = scope)
googleCloudStorageR::gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-bioheart-test")

scores_file <- "gs://cpg-bioheart-test-analysis/tenk10k/externalid_scores.csv"
metadata_file <- "gs://cpg-bioheart-test-analysis/tenk10k/metadata_tenk10k.csv"
# Copy in scores and metadata files
system(glue(
    "gsutil cp {scores_file} externalid_scores.csv"
))
system(glue(
    "gsutil cp {metadata_file} metadata_tenk10k.csv"
))
# Read in files once copied
scores <- read.csv("externalid_scores.csv")
metadata <- read.csv("metadata_tenk10k.csv")

# tidy data --------------------------------------------------------------
metadata[metadata == ""] <- NA
# remove characters before or after digits and turn to numeric values
str_strip_to_numeric <- function(character_vector) {
    as.numeric(sub("^\\D+", "", character_vector) %>% sub("(\\d)[^0-9]+$", "\\1", .))
}
character_variables <- c("Contamination..S.", "X..Aligned", "Bases...30X", "Error.rate", 
"X..Mapped", "X..Proper.Pairs", "Median.Coverage", "Insert.Size", "nt_bnp", "bnp")
for (i in character_variables){
    metadata[,i] <- str_strip_to_numeric(metadata[,i])
}

# match bioheart ID in metadata with scores
metadata <- metadata[match(scores$sample_id, metadata$sample_id), ]

# plot data --------------------------------------------------

gcs_image_outdir <- "gs://cpg-bioheart-test-web/tenk10k/"

# assign covariates to plot
sample_ids <- c("sample_id", "Sample.Name", "internal_id", "external_id", "record_id", "bioheart_id")
descriptive_columns <- colnames(select(metadata,contains("desc")))
subtract <- c(sample_ids, descriptive_columns)

# get index of unwanted variables
subtract <- which(colnames(metadata) %in% subtract)
covariate.names <- colnames(metadata)[-subtract]
all.covars.df <- metadata[,covariate.names]

# plot first ten PCs
pcs_to_plot <- seq_along(1:10)

plot.pca <- function(scores_df, metadata_df, variable_name){
    metadata_plot <- paste0("metadata_",variable_name, ".pdf")
    pdf(metadata_plot, width = 8, height = 7)
    for (i in pcs_to_plot){
        pca_axis1 = paste0("PC",i)
        pca_axis2 = paste0("PC",i+1)
        df <- data.frame(PC1 = scores_df[,pca_axis1], PC2 = scores_df[,pca_axis2], covariate = metadata_df[,variable_name])
        if (is.character(metadata_df[,variable_name])){
            p <- df %>% 
            ggplot(aes(x = PC1, y = PC2)) + geom_point(aes(fill = covariate), alpha = 0.9, shape = 21, size = 3) + 
            theme_bw() + ggtitle(variable_name) + theme(legend.title=element_blank()) +
            xlab(pca_axis1) + ylab(pca_axis2) + scale_fill_viridis_d(option = "plasma")
        } else {
            p <- df %>% 
            ggplot(aes(x = PC1, y = PC2)) + geom_point(aes(fill = covariate), alpha = 0.9, shape = 21, size = 3) + 
            theme_bw() + ggtitle(variable_name) + theme(legend.title=element_blank()) +
            xlab(pca_axis1) + ylab(pca_axis2) + scale_fill_viridis_c(option = "plasma")
        }
        print(p)
        # Copy pdf to system
        system(glue("gsutil cp {metadata_plot} {gcs_image_outdir}"))
    }
    dev.off()
}

for (name in colnames(all.covars.df)){
    plot.pca(scores_df = scores, metadata_df = metadata, variable_name = name)
}

# Perform ANOVA test ----------------------------------------------------------

pc.assoc <- function(pca.data, metadata_table){
    all.pcs <- data.frame()
    for (i in 1:ncol(pca.data)){
        all.assoc <- vector()
        for (j in 1:ncol(metadata_table)){
            test.assoc <- anova(lm(pca.data[,i] ~ metadata_table[,j]))[1,5]
            all.assoc <- c(all.assoc, test.assoc)
        }
        single.pc <- c(i, all.assoc)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", colnames(metadata_table))
    print (all.pcs)
}

# make all character variables as factors for the ANOVA
all.covars.df[sapply(all.covars.df, is.character)] <- lapply(all.covars.df[sapply(all.covars.df, is.character)], as.factor)
# drop variables that have less than one level
all.covars.df <- all.covars.df[, sapply(all.covars.df, function(col) length(unique(col))) > 2]
# select only PC scores to plot
pc_scores <- select(scores,contains("PC"))
all.pcs <- pc.assoc(pc_scores, all.covars.df)
