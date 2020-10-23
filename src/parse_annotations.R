#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(rtracklayer)

LookupParentGene <- function(cds_parent){
  mc[ID == cds_parent, unique(Parent)]}

gff_file <- snakemake@input[["gff"]]

# read the gff
gff <- import.gff(gff_file)
mc <- as.data.table(mcols(gff))

# get CDS
cds <- mc[as.character(type) == "CDS"]

# get the parent transcript for each CDS
cds[, parent_char := as.character(Parent)]
cds[, note_list := unlist(Note), by = ID]

# lookup the parent gene for each CDS
cds[, parent_gene := LookupParentGene(parent_char),
    by = parent_char]

cds_out <- unique(cds[, .(gene = parent_gene,
        protein_name = Name,
        locus_tag = locus_tag,
        transcript_id = orig_transcript_id,
        note = note_list)])

fwrite(cds_out, snakemake@output[["annot"]])

sessionInfo()
