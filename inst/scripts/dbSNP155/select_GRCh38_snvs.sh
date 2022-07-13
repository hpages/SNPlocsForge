#!/bin/sh
#
# To run this script in "batch mode":
#
#   /home/hpages/SNPlocsForge/SNPlocsForge/inst/scripts/dbSNP155/select_GRCh38_snvs.sh >select_GRCh38_snvs.log 2>&1 &
#

set -e  # Exit immediately if a simple command exits with a non-zero status

## Settings for rex3:
ASSEMBLY="GRCh38.p13"
DUMP_DIR="/home/hpages/SNPlocsForge/downloads/dbSNP155/snvs_dump"
OUT_DIR="/home/hpages/SNPlocsForge/downloads/dbSNP155/GRCh38_snvs"
Rscript="/home/hpages/R/R-4.2.0/bin/Rscript"

SEQNAMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
#SEQNAMES="22 MT"

select_snvs()
{
	chr="$1"
	dump_subdir="$DUMP_DIR/$chr"
	out_subdir="$OUT_DIR/$chr"
	mkdir $out_subdir
	Rexpr="library(SNPlocsForge)"
	Rexpr="$Rexpr; system.time(select_snvs('$dump_subdir', '$out_subdir', assembly='$ASSEMBLY'))"
	#Rexpr="$Rexpr; sessionInfo()"
	$Rscript -e "$Rexpr"
}

for seqname in $SEQNAMES; do
	select_snvs "chr${seqname}"
done

