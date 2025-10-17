#!/bin/sh
#
# To run this script in "batch mode":
#
#   cd /home/hpages/snplocs_forge/downloads/dbSNP157
#   mkdir GRCh38_snvs
#   cd GRCh38_snvs
#   /home/hpages/snplocs_forge/SNPlocsForge/inst/scripts/dbSNP157/select_GRCh38_snvs.sh >select_GRCh38_snvs.log 2>&1 &
#

set -e  # Exit immediately if a simple command exits with a non-zero status

## Settings for rex3:
ASSEMBLY="GRCh38.p14"
DUMP_DIR="/home/hpages/snplocs_forge/downloads/dbSNP157/snvs_dump"
OUT_DIR="/home/hpages/snplocs_forge/downloads/dbSNP157/GRCh38_snvs"
Rscript="/home/hpages/R/R-4.5.1/bin/Rscript"

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

