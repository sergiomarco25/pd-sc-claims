# Pseudobulk of the Wei et al. 2023 ventral-midbrain counts per donor (sampleID, demultiplexed singlets) x cell type.
# Usage: Rscript scripts/46c_wei2023_pseudobulk.R   (after scripts/46a_wei2023_extract.R has saved counts and metadata)
suppressMessages(library(Matrix))
pdrep_raw <- Sys.getenv("PDREP_RAW")  # set by slurm/run_r.sbatch; else read config/project.toml
if (pdrep_raw == "") pdrep_raw <- system2("python3", c("-I", "src/pdrep/config.py", "RAW"), stdout = TRUE)
raw <- file.path(pdrep_raw, "wei2023")
cnt <- readRDS(file.path(raw, "wei2023_counts.rds"))
meta <- read.delim(file.path(raw, "wei2023_meta.tsv"), check.names = FALSE)
meta <- meta[colnames(cnt), ]
f <- factor(paste(meta$sampleID, meta$celltype, sep = "|"))
ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
pb <- as.matrix(cnt %*% ind)
colnames(pb) <- levels(f)
write.csv(pb, gzfile(file.path(raw, "wei2023_pseudobulk.csv.gz")))
cat("pseudobulk:", dim(pb), "\n")
