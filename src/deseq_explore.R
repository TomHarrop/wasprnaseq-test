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

BiocParallel::register(
  BiocParallel::MulticoreParam(workers = 8))
dds_file <- "output/030_deseq/dds.Rds"

dds <- readRDS(dds_file)

# filter genes with low expression
row_means <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
keep_genes <- (row_means > 5 | row_max > 10)

dds_filtered <- dds[keep_genes,]
# run DESeq2
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)

vst <- varianceStabilizingTransformation(dds_filtered, blind = FALSE)

vst_counts <- assay(vst)

pc <- prcomp(t(vst_counts), center = TRUE)

pct_var <- 100 * (pc$sdev^2 / sum(pc$sdev^2))

pct_dt <- data.table(pc = paste0("PC", 1:length(pct_var)),
                     pct_var)
pc_long <- melt(data.table(pc$x, keep.rownames = TRUE),
                id.vars = "rn", variable.name = "pc", value.name = "score")

cd_dt <- as.data.table(colData(vst), keep.rownames = TRUE)

pc_pct <- merge(pc_long, pct_dt)

pc_plot <- merge(pc_pct, cd_dt, by = "rn")
n_pcs <- 4

pd <- pc_plot[pc %in% paste0("PC", 1:n_pcs)]
pd[, pc_label := paste0(pc, " (", round(pct_var, 1), "%)")]
pd[, pc_label := factor(pc_label,
                        levels = gtools::mixedsort(unique(pc_label)))]
pd[, nest := factor(as.character(nest),
                    levels = gtools::mixedsort(unique(as.character(nest))))]

gp <- ggplot(pd, aes(x = nest, y = score, colour = caste)) +
  facet_wrap(~pc_label) +
  scale_colour_viridis_d() + 
  geom_point(position = position_jitter(width = 0.2),
             shape = 16,
             alpha = 0.8)

ggsave("quick_pca.pdf", gp, width = 16, height = 9, units = "cm")

# any DE genes?
resultsNames(dds_filtered)
res_dt <- data.table(results(dds_filtered,
                             name = "caste_f_vs_d",
                             tidy = TRUE,
                             lfcThreshold = log(1.5, 2),
                             alpha = 0.05))

res_dt[order(abs(log2FoldChange), decreasing = TRUE)]

plotCounts(dds_filtered, "gene-HZH66_009572", intgroup = c("nest", "caste"))

