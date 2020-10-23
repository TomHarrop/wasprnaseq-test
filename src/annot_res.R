#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)

annot_file <- snakemake@input[["annot"]]
res_file <- snakemake@input[["res"]]

annot <- fread(annot_file)
res <- fread(res_file)

res_annot <- merge(res,
                   annot,
                   by.x = "row",
                   by.y = "gene",
                   all.x = TRUE,
                   all.y = FALSE)
res_annot[, abs_l2fc := abs(log2FoldChange)]

sortnames <- c(
  ifelse("contrast" %in% names(res_annot),
         "contrast",
         NA),
  "abs_l2fc")

sortorder <- c(1, rep(-1, length(sortnames[!is.na(sortnames)]) - 1))

setorderv(res_annot,
          sortnames[!is.na(sortnames)],
          na.last = TRUE,
          order = sortorder)

res_annot[, abs_l2fc := NULL]

fwrite(res_annot, snakemake@output[["res_annot"]])

sessionInfo()
