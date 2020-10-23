#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(deseq2)
library(data.table)
library(ggplot2)

BiocParallel::register(
  BiocParallel::MulticoreParam(workers = snakemake@threads))

# dev
BiocParallel::register(
  BiocParallel::MulticoreParam(workers = 8))
dds_file <- "output/030_deseq/dds.Rds"
alpha <- 0.1
lfc_threshold <- log(1.5, 2)

# set up
dds <- readRDS(dds_file)

# filter genes with low expression
row_means <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
keep_genes <- (row_means > 5 | row_max > 10)

dds_filtered <- dds[keep_genes,]
design(dds_filtered) <- ~ caste

# do each within-hive comparison
coldt <- as.data.table(colData(dds_filtered), keep.rownames = TRUE)

GetNestResults <- function(my_group) {
  my_dds <- dds_filtered[, my_group$rn]
  colData(my_dds) <- droplevels(colData(my_dds))
  my_dds <- DESeq(my_dds)
  my_res <- data.table(results(my_dds,
                               tidy = TRUE,
                               alpha = alpha,
                               lfcThreshold = lfc_threshold))
  return(my_res)
}

de_res_by_nest <- coldt[, GetNestResults(.SD), by = nest]
de_res_by_nest[which.min(padj)]
de_res_by_nest[order(padj)]


plotCounts(dds, "gene-HZH66_001581", intgroup = c("nest", "caste"))
 
