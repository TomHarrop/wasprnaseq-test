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
library(tximeta)

transcript_file <- snakemake@input[["mrna"]]
gff_file <- snakemake@input[["gff"]]
file_list <- snakemake@input[["quant_files"]]
index_dir <- snakemake@params[["index"]]
dds_file <- snakemake@output[["dds"]]

# dev
# transcript_file <- "output/000_ref/vvulg.mrna.fa"                 # mrna
# gff_file <- "data/ref/GCA_014466185.1_ASM1446618v1_genomic.gff"   # gff
# index_dir <- "output/005_index"
# file_list <- list.files("output/020_salmon",
#                         full.names = TRUE,
#                         recursive = TRUE,
#                         pattern = "quant.sf")

# generate the txome
genome_name <- sub("_genomic.gff", "", basename(gff_file))
tmp <- tempdir()

json_file <- file.path(tmp, paste0(genome_name, ".json"))

setTximetaBFC(tmp)
makeLinkedTxome(indexDir = index_dir,
                source = "ncbi",       # should be RefSeq but NCBI hasn't posted genome info
                organism = "Vespula vulgaris",
                release = "1",
                genome = genome_name,
                fasta = transcript_file,
                gtf = gff_file,
                jsonFile = json_file)

# generate col_data
names(file_list) <- gsub(".*/", "", dirname(file_list))
col_data <- data.table(
  files = file_list,
  names = names(file_list))
col_data[, splitname := gsub("^(n[[:digit:]]+)([d|f])([[:digit:]]+)",
                             "\\1_\\2_\\3",
                             names)]
col_data[, c("nest", "caste", "indiv") := tstrsplit(splitname, "_")]
col_data[, splitname := NULL]

# read the salmon info
se <- tximeta(col_data)

# mapping info (don't run)
# metadata(se)[["quantInfo"]]$percent_mapped
# metadata(se)[["quantInfo"]]
# metadata(se)[["quantInfo"]]$num_decoy_fragments / metadata(se)[["quantInfo"]]$num_processed * 100

# summarize counts to gene-level
gse <- summarizeToGene(se)

# generate deseq object
dds <- DESeqDataSet(gse,
                    design = ~ nest + caste )

# write output
saveRDS(dds, dds_file)

# log
sessionInfo()
