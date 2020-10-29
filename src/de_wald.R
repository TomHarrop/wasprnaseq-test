#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(DESeq2)
library(data.table)
library(ggplot2)
library(scales)

#############
# FUNCTIONS #
#############

GetMaValues <- function(my_rn,
                        dds_object){
  res <- lfcShrink(dds_object,
                   coef = my_rn,
                   type = "apeglm",
                   lfcThreshold = lfc_threshold,
                   parallel = TRUE)
  
  ma_dt <- data.table(
    plotMA(res,
           alpha = alpha,
           returnData = TRUE),
    keep.rownames = TRUE)
  return(ma_dt)
}


GetCoefContrasts <- function(my_coef) {
  my_levs <- levels(unlist(coldt[, ..my_coef]))
  comb_dt <- data.table(t(combn(my_levs, 2, simplify = TRUE)),
                        coef = my_coef)
  my_contrast_res <- comb_dt[, GetContrastResults(.SD), by = 1:nrow(comb_dt)]
  my_contrast_res[, nrow := NULL]
  return(my_contrast_res)
}

GetContrastResults <- function(my_contrast){
  res <- results(dds_filtered,
                 contrast = my_contrast[, c(coef, V1, V2)],
                 lfcThreshold = lfc_threshold,
                 alpha = alpha,
                 tidy = TRUE)
  my_res <- as.data.table(res)
  my_res[, contrast := my_contrast[, paste(coef, V1, V2, sep = "_")]]
  return(my_res)
}

# nice, but didn't use them
PlotAllMaPlots <- function(my_rn) {
  res <- lfcShrink(dds_filtered,
                   coef = my_rn,
                   type = "apeglm",
                   lfcThreshold = lfc_threshold)
  
  ma_dt <- data.table(
    plotMA(res,
           alpha = alpha,
           returnData = TRUE),
    keep.rownames = TRUE)
  
  gp <- ggplot(ma_dt, aes(x = mean,
                          y = lfc,
                          colour = isDE)) +
    ylab(expression("Log"[2] ~ "fold change")) +
    xlab("Mean of normalized counts") +
    ggtitle(my_rn) +
    scale_colour_viridis_d(guide = FALSE, direction = -1) +
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      breaks = trans_breaks("log10", function(x) 10^x)) +
    geom_point(shape = 16, alpha = 0.5)
  return(gp)
}

###########
# GLOBALS #
###########

BiocParallel::register(
  BiocParallel::MulticoreParam(workers = snakemake@threads))
dds_file <- snakemake@input[["dds"]]
ma_file <- snakemake@output[["ma"]]
wald_file <- snakemake@output[["res"]]
alpha <- snakemake@params[["alpha"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]

# dev
# BiocParallel::register(
#   BiocParallel::MulticoreParam(workers = 8))
# dds_file <- "test/star_dds.Rds"
# alpha <- 0.1
# lfc_threshold <- log(1.5, 2)
# wald_file <- "test/star_res.csv"


########
# MAIN #
########

# set up
dds <- readRDS(dds_file)

# filter genes with low expression
# row_means <- rowMeans(counts(dds))
# row_max <- rowMax(counts(dds))
# keep_genes <- (row_means > 5 | row_max > 10)

# at least 3 samples with a count of 10 or higher
# chucks out a handful more genes
keep_genes <- rowSums(counts(dds) >= 10) >= 3

dds_filtered <- dds[keep_genes,]

# run deseq
dds_filtered <- DESeq(dds_filtered,
                      parallel = TRUE)

# results for annotating (not used)
# res_dt <- data.table(results(
#   dds_filtered,
#   name = "caste_f_vs_d",
#   tidy = TRUE,
#   alpha = alpha,
#   lfcThreshold = lfc_threshold))

# generate results for all pairwise contrasts
coldt <- as.data.table(colData(dds_filtered),
                       keep.rownames = TRUE)
coefs <- c("nest", "caste")
contrast_res <- rbindlist(lapply(coefs, GetCoefContrasts))

fwrite(contrast_res,
       wald_file)

# run once per results_names
# first get the missing comparison
dds_relevel <- copy(dds_filtered)
dds_relevel$nest <- relevel(dds_relevel$nest,
                            levels(dds_relevel$nest)[[
                              length(levels(dds_relevel$nest))
                            ]])
dds_relevel <- DESeq(dds_relevel,
                     parallel = TRUE)

# get the results from each dds run
rn1 <- resultsNames(dds_filtered)[-1]
names(rn1) <- rn1
rn2 <- resultsNames(dds_relevel)[-1]
names(rn2) <- rn2

# work out which resultsNames in dds_relevel are already in dds_filtered
rn1_names <- sapply(rn1, function(x) sort(unlist(strsplit(x, "_"))))
rn1_tr <- apply(rn1_names, 2, paste, collapse = "_")

rn2_names <- sapply(rn2, function(x) sort(unlist(strsplit(x, "_"))))
rn2_tr <- apply(rn2_names, 2, paste, collapse = "_")
rn2_subset <- rn2[names(rn2_tr[!rn2_tr %in% rn1_tr])]

# get the MA values for each resultsName
set1 <- rbindlist(lapply(rn1, GetMaValues, dds_object = dds_filtered),
                  idcol = "contrast")
set2 <- rbindlist(lapply(rn2_subset, GetMaValues, dds_object = dds_relevel),
                  idcol = "contrast")

all_ma <- rbindlist(list(set1, set2))

# generate a plot
gp <- ggplot(all_ma, aes(x = mean,
                   y = lfc,
                   colour = isDE)) +
  ylab(expression("Log"[2] ~ "fold change")) +
  xlab("Mean of normalized counts") +
  facet_wrap(~contrast) +
  scale_colour_viridis_d(guide = FALSE, direction = -1) +
  scale_x_continuous(
    trans = "log10",
    labels = trans_format("log10", math_format(10^.x)),
    breaks = trans_breaks("log10", function(x) 10^x)) +
  geom_hline(yintercept = c(lfc_threshold, -lfc_threshold),
             linetype = 2) +
  geom_point(shape = 16, alpha = 0.5, size = 1)

ggsave(ma_file,
       gp,
       width = 16,
       height = 9,
       units = "cm")

sessionInfo()
