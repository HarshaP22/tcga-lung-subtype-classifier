#!/usr/bin/env bash
# =============================================================================
# One-command data downloader for the LUAD-vs-LUSC overlapping-feature-selection
# reproduction. Downloads ONLY public reference datasets into ./data/.
#
#   Usage:   bash download_data.sh
#
# - Re-runnable: files already present (non-empty) are skipped.
# - Self-healing: tries several candidate URLs per file (e.g. STAR vs HTSeq
#   counts, gencode v36 vs v22) and keeps the first that works.
# - Touches nothing outside this project's data/ folders. No git, no accounts.
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW="$SCRIPT_DIR/data/raw"
EXT="$SCRIPT_DIR/data/external"
GS="$RAW/genesets"
mkdir -p "$RAW" "$EXT" "$GS"

XENA="https://gdc-hub.s3.us-east-1.amazonaws.com/download"
GEO="https://ftp.ncbi.nlm.nih.gov/geo"
ENR="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName"
# Xena's official Ensembl->symbol probeMap lives in this public repo (not an account):
GHXENA="https://raw.githubusercontent.com/ucscXena/xena-GDC-ETL/master/xena_gdc_etl/resources"
ACC="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

OK=(); MISS=()

# fetch <dest> <url1> [url2 ...] : try each url until one succeeds
fetch () {
  local dest="$1"; shift
  local name; name="$(basename "$dest")"
  if [ -s "$dest" ]; then echo "  ✓ already have $name"; OK+=("$name"); return 0; fi
  local url
  for url in "$@"; do
    echo "  … $name  <-  $url"
    if curl -fL --retry 3 --retry-delay 2 -C - \
            --connect-timeout 30 -A "Mozilla/5.0 (luad-lusc-dl)" \
            -o "$dest.part" "$url" 2>/dev/null && [ -s "$dest.part" ]; then
      mv "$dest.part" "$dest"
      # integrity check for gzip files
      case "$dest" in
        *.gz) if ! gzip -t "$dest" 2>/dev/null; then
                echo "    ! downloaded but not a valid gzip, trying next candidate"
                rm -f "$dest"; continue; fi ;;
      esac
      echo "    -> saved ($(du -h "$dest" | cut -f1))"
      OK+=("$name"); return 0
    fi
    rm -f "$dest.part"
  done
  echo "    ✗ could not fetch $name from any candidate URL"
  MISS+=("$name"); return 1
}

echo "==============================================================="
echo " A. CORE: TCGA expression (needed to start)"
echo "==============================================================="
fetch "$RAW/TCGA-LUAD.star_counts.tsv.gz"  "$XENA/TCGA-LUAD.star_counts.tsv.gz"  "$XENA/TCGA-LUAD.htseq_counts.tsv.gz"
fetch "$RAW/TCGA-LUSC.star_counts.tsv.gz"  "$XENA/TCGA-LUSC.star_counts.tsv.gz"  "$XENA/TCGA-LUSC.htseq_counts.tsv.gz"
fetch "$RAW/gencode.v22.annotation.gene.probeMap" \
      "$GHXENA/gencode.v22.annotation.gene.probeMap" \
      "$XENA/gencode.v36.annotation.gene.probeMap" \
      "$XENA/gencode.v22.annotation.gene.probeMap"

echo "==============================================================="
echo " B. EXTENSIONS: phenotype + survival"
echo "==============================================================="
fetch "$RAW/TCGA-LUAD.clinical.tsv.gz" "$XENA/TCGA-LUAD.clinical.tsv.gz" "$XENA/TCGA-LUAD.GDC_phenotype.tsv.gz"
fetch "$RAW/TCGA-LUSC.clinical.tsv.gz" "$XENA/TCGA-LUSC.clinical.tsv.gz" "$XENA/TCGA-LUSC.GDC_phenotype.tsv.gz"
fetch "$RAW/TCGA-LUAD.survival.tsv.gz"      "$XENA/TCGA-LUAD.survival.tsv.gz"
fetch "$RAW/TCGA-LUSC.survival.tsv.gz"      "$XENA/TCGA-LUSC.survival.tsv.gz"

echo "==============================================================="
echo " C. EXTERNAL VALIDATION: GSE28582 (GEO)"
echo "==============================================================="
fetch "$EXT/GSE28582-GPL570_series_matrix.txt.gz" \
      "$GEO/series/GSE28nnn/GSE28582/matrix/GSE28582-GPL570_series_matrix.txt.gz"
# GEO has no curated .annot for GPL570, so pull the platform's full annotation
# table (plain text SOFT: includes an "ID" and "Gene Symbol" column).
fetch "$RAW/GPL570.annot.txt" \
      "$ACC?acc=GPL570&targ=self&form=text&view=full"

echo "==============================================================="
echo " D. PATHWAY LIBRARIES: GO / KEGG / Reactome (.gmt)"
echo "==============================================================="
fetch "$GS/GO_Biological_Process_2021.gmt" "$ENR=GO_Biological_Process_2021"
fetch "$GS/KEGG_2021_Human.gmt"            "$ENR=KEGG_2021_Human"
fetch "$GS/Reactome_2022.gmt"              "$ENR=Reactome_2022"

echo
echo "==============================================================="
echo " SUMMARY"
echo "==============================================================="
echo " downloaded/present: ${#OK[@]} file(s)"
if [ "${#MISS[@]}" -gt 0 ]; then
  echo " MISSING (${#MISS[@]}): ${MISS[*]}"
  echo
  echo " For any missing TCGA file, use the always-current UCSC Xena page:"
  echo "   LUAD: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Adenocarcinoma%20(LUAD)&host=https%3A%2F%2Fgdc.xenahubs.net"
  echo "   LUSC: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&host=https%3A%2F%2Fgdc.xenahubs.net"
else
  echo " All files present. Next: make reproduce"
fi
echo
echo " Files landed in:"
echo "   $RAW"
echo "   $EXT"
echo "   $GS"
