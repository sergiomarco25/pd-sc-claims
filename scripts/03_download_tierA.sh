#!/bin/bash
# Download the open tier-A datasets that were not yet on disk (R1, 2026-10-21).
# Each dataset goes to its own folder under $RAW/<paper_id>/. Already-complete files are skipped (wget -c).
# For large GEO series only the per-sample count matrices are fetched (ATAC fragment files skipped, ~46 GB).
# Usage: bash scripts/03_download_tierA.sh [paper_id ...]   (no argument = all)
set -u
eval "$(python3 -I "$(dirname "$0")/../src/pdrep/config.py" --env)"   # paths from config/project.toml
RAW=$PDREP_RAW
GEO=https://ftp.ncbi.nlm.nih.gov/geo
log(){ echo "$(date '+%F %T') $*"; }
fetch(){ # url outdir
  local f; f=$(basename "${1%%\?*}"); mkdir -p "$2"
  if wget -q -c -t 5 -O "$2/$f" "$1"; then log "OK $2/$f"; else log "FAIL $1"; fi; }
series_file(){ local d; d=$(echo $1 | sed 's/...$/nnn/'); fetch "$GEO/series/$d/$1/suppl/$2" "$3"; }
soft(){ local d; d=$(echo $1 | sed 's/...$/nnn/'); fetch "$GEO/series/$d/$1/soft/$1_family.soft.gz" "$2"; }
gsm_files(){ # acc regex outdir : per-sample files whose name matches regex, from the series filelist
  local d; d=$(echo $1 | sed 's/...$/nnn/')
  curl -s "$GEO/series/$d/$1/suppl/filelist.txt" | awk -F'\t' 'NR>1 && $1=="File"{print $2}' | grep -E "$2" |
  while read -r f; do local gsm=${f%%_*}; local gd; gd=$(echo $gsm | sed 's/...$/nnn/')
    fetch "$GEO/samples/$gd/$gsm/suppl/$f" "$3"; done; }

ma2025(){ local o=$RAW/ma2025; soft GSE253462 $o; series_file GSE253462 GSE253462_RAW.tar $o
  soft GSE253975 $o/visium; series_file GSE253975 GSE253975_RAW.tar $o/visium; }
adams2024(){ local o=$RAW/adams2024; soft GSE193688 $o; gsm_files GSE193688 'filtered_feature_bc_matrix\.h5$|atac_peaks\.bed\.gz$' $o; }
chia2025(){ local o=$RAW/chia2025; soft GSE303823 $o; series_file GSE303823 GSE303823_rna_annotated.rds $o; }
dehestani2024(){ local o=$RAW/dehestani2024; soft GSE272760 $o; series_file GSE272760 GSE272760_all_new_SCT_harmony_res0.1_v1.rds.gz $o; }
feleke2021(){ local o=$RAW/feleke2021; soft GSE178146 $o; series_file GSE178146 GSE178146_RAW.tar $o; }
glock2026(){ local o=$RAW/glock2026; soft GSE295998 $o; series_file GSE295998 GSE295998_UmapCoordinatesAndClusterAssignments.txt.gz $o
  gsm_files GSE295998 'cellbender_filtered\.h5$' $o; }
jin2024(){ local o=$RAW/jin2024; soft GSE235914 $o; series_file GSE235914 GSE235914_RAW.tar $o; }
lee2023(){ local o=$RAW/lee2023; soft GSE148434 $o; series_file GSE148434 GSE148434_RNA.log2TPM.combat.qq.txt.gz $o
  gsm_files GSE148434 'snRNA\.molecule_info\.h5$|snATAC\.filtered_peak_bc_matrix\.h5$|snATAC\.singlecell\.csv\.gz$' $o; }
volpato2025(){ local o=$RAW/volpato2025; soft GSE329625 $o
  series_file GSE329625 GSE329625_full_metadata.csv.gz $o; series_file GSE329625 GSE329625_processed_gene_matrix.txt.gz $o; }
psychad(){ # CELLxGENE PsychAD MSSM cohort (contains the PD donors); subset later
  local o=$RAW/psychad_mssm; mkdir -p $o
  curl -s https://api.cellxgene.cziscience.com/curation/v1/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f > $o/collection.json
  local url; url=$(python3 -c "import json,sys; c=json.load(open(sys.argv[1]))
print([a['url'] for d in c['datasets'] if d['dataset_id']=='37a17b78-4864-4a42-b67b-31c00962795a' for a in d['assets'] if a['filetype']=='H5AD'][0])" $o/collection.json)
  fetch "$url" $o; }

ALL="ma2025 adams2024 chia2025 dehestani2024 feleke2021 glock2026 jin2024 lee2023 volpato2025 psychad"
for ds in ${@:-$ALL}; do log "== $ds"; $ds; done
log DONE
