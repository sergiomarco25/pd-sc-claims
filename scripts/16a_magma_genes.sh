#!/bin/bash
# PD GWAS gene-level association with MAGMA v1.10 (Nalls 2019, GCST009325, GRCh37; 1000G EUR reference).
# Sumstats carry no rsIDs: variants are matched to the reference .bim by chr, position and alleles.
# Output: results/gwas/nalls2019.genes.out / .genes.raw (gene-level Z, input to gene-set tests per dataset).
set -euo pipefail
cd "$(dirname "$0")/.."
eval "$(python3 -I src/pdrep/config.py --env)"   # paths from config/project.toml
G=$PDREP_RAW/gwas
O=$PDREP_RESULTS/gwas
mkdir -p $O $G/magma_bin/x $G/magma_ref/x
[ -x $G/magma_bin/x/magma ] || unzip -oq $G/magma_bin/magma_static.zip -d $G/magma_bin/x
[ -s $G/magma_ref/x/NCBI37.3.gene.loc ] || unzip -oq $G/magma_ref/gene_loc_b37.zip -d $G/magma_ref/x
[ -s $G/magma_ref/x/g1000_eur.bim ] || unzip -oq $G/magma_ref/g1000_eur.zip -d $G/magma_ref/x
M=$G/magma_bin/x/magma; R=$G/magma_ref/x
"$PDREP_PYTHON" -I scripts/16b_magma_pfile.py $G/nalls2019/GCST009325.tsv $R/g1000_eur.bim $O/nalls2019.pval
$M --annotate window=35,10 --snp-loc $R/g1000_eur.bim --gene-loc $R/NCBI37.3.gene.loc --out $O/annot
$M --bfile $R/g1000_eur --pval $O/nalls2019.pval ncol=N --gene-annot $O/annot.genes.annot --out $O/nalls2019
