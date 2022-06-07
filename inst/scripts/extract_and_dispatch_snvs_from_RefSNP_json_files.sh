#!/bin/sh
#
# Run this script in batch mode with:
#
#   path/to/extract_and_dispatch_snvs_from_RefSNP_json_files.sh >extract_and_dispatch_snvs_from_RefSNP_json_files.log 2>&1 &
#

set -e  # exit immediately if a simple command returns a non-zero status

# Settings for dbSNP155 on rex3:
JSON_DIR="/home/hpages/SNPlocsForge/downloads/dbSNP155"
DUMP_DIR="$JSON_DIR/snvs_dump"
JSON_FILES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12"
JSON_PREFIX="refsnp-"
JSON_SUFFIX=".json.bz2"
Rscript="/home/hpages/R/R-4.2.0/bin/Rscript"

mkdir "$DUMP_DIR"

cd "$JSON_DIR"

for file in $JSON_FILES; do
	filename="${JSON_PREFIX}${file}${JSON_SUFFIX}"
	echo "Start processing ${filename}."
	logfile="$DUMP_DIR/${seqname}.log"
	Rexpr="library(SNPlocsForge)"
	Rexpr="$Rexpr;system.time(extract_and_dispatch_snvs_from_RefSNP_json('$filename', '$DUMP_DIR'))"
	$Rscript -e "$Rexpr" >$logfile 2>&1 &
done

echo "DONE."
echo ""

