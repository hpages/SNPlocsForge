#!/bin/sh
#
# To run this script in "batch mode":
#
#   cd /home/hpages/SNPlocsForge/forged/SNPlocs.Hsapiens.dbSNP155.GRCh38/inst/extdata
#   /home/hpages/SNPlocsForge/SNPlocsForge/inst/scripts/dbSNP155/build_GRCh38_OnDiskLongTable.sh >build_OnDiskLongTable.log 2>&1 &
#

set -e  # Exit immediately if a simple command exits with a non-zero status

## Settings for rex3:
PKGNAME="SNPlocs.Hsapiens.dbSNP155.GRCh38"
ASSEMBLY="GRCh38.p13"
SELECTED_SNVS_DIR="/home/hpages/SNPlocsForge/downloads/dbSNP155/GRCh38_snvs"
Rscript="/home/hpages/R/R-4.2.0/bin/Rscript"

SEQNAMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
#SEQNAMES="22 MT"

Rexpr="library(SNPlocsForge)"
Rexpr="$Rexpr; system.time(build_OnDiskLongTable('$SELECTED_SNVS_DIR', '$SEQNAMES', assembly='$ASSEMBLY'))"
Rexpr="$Rexpr; sessionInfo()"
$Rscript -e "$Rexpr"

