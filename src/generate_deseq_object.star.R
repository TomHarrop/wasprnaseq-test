#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)

star_files <- snakemake@input[["quant_files"]]
dds_file <- snakemake@output[["dds"]]

names(star_files) <- sub(".ReadsPerGene.out.tab", "", basename(star_files))

star_list <- lapply(star_files, fread)
star_counts_raw <- rbindlist(star_list, idcol = "sample")
star_counts <- star_counts_raw[!startsWith(V1, "N_")]

# V3 and V4 are the strand-specific cols
# column 1:  gene ID
# column 2:  counts for unstranded RNA-seq
# column 3:  counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4:  counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
# illumina truseq stranded == column 4
x <- star_counts[, lapply(.SD, sum), .SDcols = c("V2", "V3", "V4"),
            by = sample]
y <- x[,.SD / V2, .SDcols = c("V3", "V4"), by = sample]

# generate count matrix
# V4 is the strand-specific column for illumina TruSeq
sample_order <- star_counts[, gtools::mixedsort(unique(sample))]
count_dt <- dcast(star_counts, V1 ~ sample, value.var = "V4")
setcolorder(count_dt, sample_order)

# generate col_data
col_data <- data.table(
  names = sample_order
)
col_data[, c("caste", "indiv") := tstrsplit(names, "_")]

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(data.frame(count_dt, row.names = "V1")),
  colData = data.frame(col_data, row.names = "names"),
  design = ~ caste)

# write output
saveRDS(dds, dds_file)

# log
sessionInfo()
