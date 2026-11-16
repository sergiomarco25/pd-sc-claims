#!/usr/bin/env Rscript
# DE cache: edgeR quasi-likelihood pseudobulk DE for every region x cell class x contrast of one dataset that has
# >= 3 case and >= 3 reference units (manifest.tsv from 12_export_pseudobulk.py), plus a continuous Lewy-body stage
# fit ("braak") across all units with a stage value. Model as revision/scripts/05d_edger_pseudobulk.R:
#   ~ group + age + sex + pmi   (covariates used only when known for every unit in the fit and not constant)
# Usage: Rscript 12b_de_cache.R <dataset_dir>   (dataset_dir = results/de_cache/<dataset>)
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
dsdir <- commandArgs(trailingOnly = TRUE)[1]
indir <- file.path(dsdir, "input"); outdir <- file.path(dsdir, "de"); dir.create(outdir, showWarnings = FALSE)
units <- fread(file.path(indir, "units.tsv"), colClasses = "character")
man <- fread(file.path(indir, "manifest.tsv"))
scl <- function(v) as.numeric(scale(v))
braak_of <- function(u) {
  b <- suppressWarnings(as.numeric(u$braak_lb))
  if (all(is.na(b)) && "unified_lb_stage" %in% names(u)) b <- suppressWarnings(as.numeric(u$unified_lb_stage))
  if (any(!is.na(b))) b[is.na(b) & u$disease == "control"] <- 0  # as pdrep.contrasts.braak_value
  b
}

fit_one <- function(counts, u, term) {
  covs <- c()
  if (all(!is.na(suppressWarnings(as.numeric(u$age)))) && sd(as.numeric(u$age)) > 0) { u$age_s <- scl(as.numeric(u$age)); covs <- c(covs, "age_s") }
  if (all(u$sex %in% c("M", "F")) && length(unique(u$sex)) > 1) { u$sex <- factor(u$sex); covs <- c(covs, "sex") }
  if (all(!is.na(suppressWarnings(as.numeric(u$pmi)))) && sd(as.numeric(u$pmi)) > 0) { u$pmi_s <- scl(as.numeric(u$pmi)); covs <- c(covs, "pmi_s") }
  # keep covariates only while the design stays estimable with >= 2 residual df
  repeat {
    design <- model.matrix(as.formula(paste("~", paste(c(term, covs), collapse = " + "))), data = u)
    if (qr(design)$rank == ncol(design) && nrow(design) - ncol(design) >= 2) break
    if (!length(covs)) return(NULL)
    covs <- head(covs, -1)
  }
  y <- DGEList(counts); keep <- filterByExpr(y, design)
  if (sum(keep) < 50) return(NULL)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE)
  coef <- grep(paste0("^", term), colnames(design), value = TRUE)[1]
  tt <- as.data.table(topTags(glmQLFTest(fit, coef = coef), n = Inf)$table, keep.rownames = "gene")
  setnames(tt, c("PValue", "FDR"), c("p", "padj"))
  tt[, covariates := paste(covs, collapse = "+")]
  tt
}

todo <- unique(man[run == TRUE, .(region, cell_class, file)])
all_files <- unique(man[, .(region, cell_class, file)])
log <- list()
for (i in seq_len(nrow(all_files))) {
  f <- all_files[i]
  counts <- fread(cmd = paste("zcat", shQuote(file.path(indir, f$file))))
  gn <- counts[[1]]; counts <- as.matrix(counts[, -1, with = FALSE]); rownames(counts) <- gn
  u <- units[match(colnames(counts), units$unit)]
  stem <- sub("\\.csv\\.gz$", "", f$file)
  for (con in man[file == f$file & run == TRUE, contrast]) {
    g <- u[[con]]; sel <- g %in% c("case", "ref")
    uu <- u[sel]; uu$group <- factor(g[sel], levels = c("ref", "case"))
    tt <- tryCatch(fit_one(counts[, sel, drop = FALSE], uu, "group"), error = function(e) { message(stem, " ", con, ": ", conditionMessage(e)); NULL })
    if (!is.null(tt)) {
      fwrite(tt, file.path(outdir, paste0(stem, "__", con, ".tsv.gz")), sep = "\t")
      log[[length(log) + 1]] <- data.table(file = f$file, region = f$region, cell_class = f$cell_class, contrast = con,
                                           n_case = sum(uu$group == "case"), n_ref = sum(uu$group == "ref"),
                                           n_genes = nrow(tt), n_fdr05 = sum(tt$padj < 0.05), covariates = tt$covariates[1])
    }
  }
  b <- braak_of(u); sel <- !is.na(b)
  if (sum(sel) >= 6 && length(unique(b[sel])) >= 3) {
    uu <- u[sel]; uu$stage <- b[sel]
    tt <- tryCatch(fit_one(counts[, sel, drop = FALSE], uu, "stage"), error = function(e) NULL)
    if (!is.null(tt)) {
      fwrite(tt, file.path(outdir, paste0(stem, "__braak.tsv.gz")), sep = "\t")
      log[[length(log) + 1]] <- data.table(file = f$file, region = f$region, cell_class = f$cell_class, contrast = "braak",
                                           n_case = sum(sel), n_ref = NA, n_genes = nrow(tt), n_fdr05 = sum(tt$padj < 0.05),
                                           covariates = tt$covariates[1])
    }
  }
}
res <- rbindlist(log)
if (!nrow(res)) res <- data.table(file = character(), region = character(), cell_class = character(), contrast = character(),
                                  n_case = integer(), n_ref = integer(), n_genes = integer(), n_fdr05 = integer(), covariates = character())
fwrite(res, file.path(dsdir, "de_log.tsv"), sep = "\t")
message(basename(dsdir), ": ", nrow(res), " fits written")
if (nrow(res)) print(res[, .(fits = .N, median_fdr05 = median(n_fdr05)), by = contrast])
