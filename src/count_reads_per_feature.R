#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")


library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(Rsamtools)
library(rtracklayer)
library(systemPipeR)

#############
# FUNCTIONS #
#############

FirstElement <- function(y){unlist(y)[1][1]}

ReadCountWrapper <- function(features,
                             bamfile_list) {
  my_counts <- GenomicAlignments::summarizeOverlaps(
    features = features,
    reads = bamfile_list,
    mode = "Union",
    singleEnd = FALSE,
    ignore.strand = TRUE,
    fragments = FALSE)
  my_rd <- rowData(my_counts)
  my_unlisted_rd_list <- lapply(seq_len(ncol(my_rd)), function(i)
    sapply(my_rd[,i], FirstElement))
  names(my_unlisted_rd_list) <- names(my_rd)
  my_unlisted_rd <- data.table(do.call(cbind, my_unlisted_rd_list))
  data.table(my_unlisted_rd,
             data.table(assay(my_counts)))
}

###########
# GLOBALS #
###########

bamfile_list <- snakemake@input[["bam_files"]]
gff_file <- snakemake@input[["gff"]]
summary_file <- snakemake@output[["summary"]]
fc_file <- snakemake@output[["feature_counts"]]
cpus <- snakemake@threads[[1]]

# cpus <- 8
# gff_file = "data/ref/GCA_014466185.1_ASM1446618v1_genomic.gff"
# bamfile_list <- list.files("output/025_star/pass2",
#            pattern = "Aligned.sortedByCoord.out.bam",
#            full.names = TRUE)[1:2]
# plot1_file <- "test/counts_per_category.pdf"
# plot2_file <- "test/intron_exon_counts.pdf"

########
# MAIN #
########

# set up multiprocessing
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load bamfiles
names(bamfile_list) <- gsub("\\..*$", "", basename(bamfile_list))
bamfiles <- BamFileList(bamfile_list)

# prepare annotations
gff <- rtracklayer::import.gff(gff_file)

# found the nuclear rRNA genes
ribo_gff <- subset(gff, type == "rRNA")
ribo_parents <- unique(unlist(ribo_gff$Parent))

# get all nuclear annotations
chr_mask <- !startsWith(as.character(seqnames(gff)), "JACSEA")
nuclear_gff <- subset(gff, chr_mask) 
nuclear_txdb_with_rrna <- GenomicFeatures::makeTxDbFromGRanges(nuclear_gff)

# remove rRNA genes from nuclear annotations
# there aren't any, so skip
# chuck <- apply(sapply(ribo_parents, function(x)
#   grepl(x, nuclear_gff$ID)), 1, any)
# nuc_no_rRNA <- nuclear_gff[!chuck]
nuc_no_rRNA <- copy(nuclear_gff)

# generate nuclear-only non-rRNA txdb
nuclear_txdb <- GenomicFeatures::makeTxDbFromGRanges(nuc_no_rRNA)
nuc_exons <- GenomicFeatures::exonicParts(nuclear_txdb,
                                          linked.to.single.gene.only = TRUE)
nuc_introns <- GenomicFeatures::intronicParts(nuclear_txdb,
                                              linked.to.single.gene.only = TRUE)

# mito regions... not defined, skip
# mito_gff <- subset(gff, seqnames == "ChrM")
# mito_txdb <- GenomicFeatures::makeTxDbFromGRanges(mito_gff)
# mito_exons <- GenomicFeatures::exonicParts(mito_txdb,
#                                            linked.to.single.gene.only = TRUE)


# intergenic
# https://support.bioconductor.org/p/73648/
gen_features <- genFeatures(nuclear_txdb_with_rrna,
                            reduce_ranges = TRUE)

mrna <- disjoin(gen_features$mRNA_red)
mrna_gaps <- disjoin(gaps(mrna))
intergenic <- disjoin(gen_features$intergenic)
mrna_plus1kb <- disjoin(mrna + 500)
mrna_plus1kb_gaps <- disjoin(gaps(mrna_plus1kb))

# list of things to count
feature_list <- list(exons = nuc_exons,
                     introns = nuc_introns,
                     mrna = mrna,
                     mrna_gaps <- mrna_gaps,
                     intergenic = intergenic,
                     mrna_plus1kb = mrna_plus1kb,
                     mrna_plus1kb_gaps = mrna_plus1kb_gaps)

feature_count_list <- lapply(feature_list,
                             ReadCountWrapper,
                             bamfile_list = bamfiles)

feature_counts <- rbindlist(feature_count_list,
                            idcol = "feature",
                            fill = TRUE,
                            use.names = TRUE)

# PLOT PER CATEGORY
# feature_counts_long <- melt(feature_counts,
#                             id.vars = c("feature",
#                                         "gene_id",
#                                         "exon_name",
#                                         "feature_by"),
#                             measure.vars = names(bamfile_list),
#                             variable.name = "sample_name",
#                             value.name = "counts")

# pd <- feature_counts_long[, .(counts = sum(counts)),
#                           by = .(feature, sample_name)]
# pd[, counts_per_sample := sum(counts), by = sample_name]
# pd[, percentage_of_counts := counts* 100 / counts_per_sample]

# write output
fwrite(feature_counts, fc_file)
# fwrite(feature_counts_long, fc_file)
# fwrite(pd, summary_file)

# write session info
sessionInfo()
