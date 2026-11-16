# Extract counts and metadata from the Wei et al. 2023 ventral-midbrain Seurat object (opioid overdose vs controls)
# without Seurat installed: S4 slots are read as attributes. Writes per donor x cell type pseudobulk counts.
# Usage: Rscript scripts/46a_wei2023_extract.R
suppressMessages(library(Matrix))
pdrep_raw <- Sys.getenv("PDREP_RAW")  # set by slurm/run_r.sbatch; else read config/project.toml
if (pdrep_raw == "") pdrep_raw <- system2("python3", c("-I", "src/pdrep/config.py", "RAW"), stdout = TRUE)
raw <- file.path(pdrep_raw, "wei2023")
x <- readRDS(gzcon(file(file.path(raw, "GSE240457_SeuratObject-GEO.rds.gz"), "rb")))
a <- attributes(x)
meta <- a$meta.data
cat("cells:", nrow(meta), "\n"); print(colnames(meta)); print(head(meta, 3))
rna <- attributes(a$assays$RNA)
cnt <- if (!is.null(rna$counts)) rna$counts else rna$layers$counts
if (is.null(rownames(cnt))) stop("counts without gene names")
cat("counts:", dim(cnt), class(cnt), "\n")
write.table(meta, file.path(raw, "wei2023_meta.tsv"), sep = "\t", quote = FALSE)
saveRDS(cnt, file.path(raw, "wei2023_counts.rds"))
# pseudobulk by every donor-like x cell-type-like column pair present
dcol <- intersect(c("sampleID", "Donor", "donor", "donor_id", "DonorID", "Sample", "sample"), colnames(meta))
ccol <- intersect(c("celltype", "CellType", "cell_type", "Celltype", "predicted.celltype", "seurat_clusters"), colnames(meta))
cat("donor columns:", dcol, " cell-type columns:", ccol, "\n")
key <- paste(meta[[dcol[1]]], meta[[ccol[1]]], sep = "|")
f <- factor(key)
ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
pb <- as.matrix(cnt %*% ind)
colnames(pb) <- levels(f)
write.csv(pb, gzfile(file.path(raw, "wei2023_pseudobulk.csv.gz")))
cat("pseudobulk:", dim(pb), "\n")
