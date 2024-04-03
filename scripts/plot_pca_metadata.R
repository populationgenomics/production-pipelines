#!/usr/bin/env Rscript

# Overlay metadata as coloured points onto PCA

library(tidyverse)
library(ggplot2)
library(glue)
library(googleCloudStorageR)
library(gargle)
library(magrittr)

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
# tidy data
metadata[metadata == ""] <- NA

# match bioheart ID in metadata with scores
metadata = metadata[match(scores$external_id, metadata$bioheart_id), ]

# plot data
gcs_image_outdir <- "gs://cpg-bioheart-test-web/tenk10k/"
pcs_to_plot=seq_along(1:10)

plot.pca <- function(scores_df, metadata_df, variable_name){
    for (i in pcs_to_plot){
        pca_axis1=paste0("PC",i)
        pca_axis2=paste0("PC",i+1)
        df <- data.frame(PC1 = scores_df[,pca_axis1], PC2 = scores_df[,pca_axis2], covariate = metadata_df[,variable_name])
        p <- df %>% 
        ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(fill=covariate), alpha=0.9, shape=21, size=3) + 
        theme_bw() + ggtitle(variable_name) + theme(legend.title=element_blank()) +
        xlab(pca_axis1) + ylab(pca_axis2)
        # Save plot
        metadata_plot <- paste0("metadata_",pca_axis1,"_",variable_name, ".pdf")
        pdf(metadata_plot,  width = 8, height = 7)
        print(p)
        dev.off()
        # Copy pdf to system
        system(glue("gsutil cp {metadata_plot} {gcs_image_outdir}"))
    }
}

# Plot PCA
for (name in colnames(metadata)){
    plot.pca(scores_df=scores, metadata_df=metadata, variable_name=name)
}
