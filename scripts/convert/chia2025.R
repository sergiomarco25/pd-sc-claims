# Convert the chia2025 Seurat object (GSE303823_rna_annotated.rds) to Matrix Market + TSV:
# RNA assay raw counts (counts slot / counts layers; NOT SCT, NOT data) and meta.data. Works without SeuratObject by
# reading the S4 slots directly via attr() (slotNames/@ need the class definition) (Seurat v3/v4 Assay: @counts; Seurat v5 Assay5: @layers$counts*).
# Usage: Rscript scripts/convert/chia2025.R <in.rds[.gz]> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(Matrix))
t0 <- Sys.time()
# GEO file for dehestani2024 is gzip-compressed twice (.rds.gz around a gzip .rds): fall back to two layers
obj <- tryCatch(readRDS(infile), error = function(e) {
  cat("readRDS failed (", conditionMessage(e), "); trying double-gzip\n")
  con <- gzcon(gzfile(infile, "rb")); on.exit(close(con)); readRDS(con) })
cat("class:", class(obj), " read in", format(Sys.time() - t0), "\n")
S <- function(x, n) attr(x, n, exact = TRUE)
sn <- function(x) setdiff(names(attributes(x)), "class")
cat("slots:", paste(sn(obj), collapse = ","), "\n")
cat("assays:", paste(names(S(obj, "assays")), collapse = ","), "\n")
a <- S(obj, "assays")[["RNA"]]
cat("RNA assay class:", class(a), " slots:", paste(sn(a), collapse = ","), "\n")
dn <- function(x) { d <- S(x, "dimnames"); m <- matrix(as.logical(unclass(asS4(x, FALSE, complete = FALSE))), nrow = length(d[[1]]), dimnames = d); m }
if ("layers" %in% sn(a)) {
  ln <- grep("^counts", names(S(a, "layers")), value = TRUE)
  cat("counts layers:", paste(ln, collapse = ","), "\n")
  feat <- dn(S(a, "features")); cells <- dn(S(a, "cells"))
  mats <- lapply(ln, function(l) {
    m <- S(a, "layers")[[l]]
    rownames(m) <- rownames(feat)[feat[, l]]
    colnames(m) <- rownames(cells)[cells[, l]]
    m
  })
  if (length(mats) > 1) {
    g <- Reduce(union, lapply(mats, rownames))
    mats <- lapply(mats, function(m) { z <- Matrix(0, length(setdiff(g, rownames(m))), ncol(m), sparse = TRUE);
      rownames(z) <- setdiff(g, rownames(m)); colnames(z) <- colnames(m); rbind(m, z)[g, ] })
  }
  counts <- do.call(cbind, mats)
} else {
  counts <- S(a, "counts")
}
counts <- as(counts, "CsparseMatrix")
meta <- S(obj, "meta.data")
cat("counts:", dim(counts), " meta:", dim(meta), " sum:", sum(counts@x), "\n")
cat("non-integer values in counts:", sum(counts@x != round(counts@x)), "\n")
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop = FALSE]
writeMM(counts, file.path(out, "counts.mtx"))
writeLines(rownames(counts), file.path(out, "genes.tsv"))
writeLines(colnames(counts), file.path(out, "barcodes.tsv"))
write.table(meta, file.path(out, "meta.tsv"), sep = "\t", quote = TRUE, col.names = NA)
# also record reductions/other assays present, for the notes
cat("reductions:", paste(names(S(obj, "reductions")), collapse = ","), "\n")
system(paste("gzip -f", shQuote(file.path(out, "counts.mtx"))))
cat("done in", format(Sys.time() - t0), "\n")
