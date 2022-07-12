#!/bin/bash
#

set -e  # exit immediately if a simple command returns a non-zero status

# Settings for dbSNP155 on rex3 (80 logical cpus), will extract 58G of data:
JSON_DIR="/home/hpages/SNPlocsForge/downloads/dbSNP155"
DUMP_DIR="$JSON_DIR/snvs_dump"
JSON_PREFIX="refsnp-"
JSON_SUFFIX=".json.bz2"
Rscript="/home/hpages/R/R-4.2.0/bin/Rscript"

cd "$JSON_DIR"

start_extraction()
{
	chr="$1"
	nworkers="$2"
	chunksize="$3"
	jsonfile="${JSON_PREFIX}${chr}${JSON_SUFFIX}"
	dump_subdir="$DUMP_DIR/$chr"
	mkdir $dump_subdir
	logfile="$DUMP_DIR/${chr}.log"
	echo "Start processing ${jsonfile}."
	Rexpr="library(SNPlocsForge)"
	Rexpr="$Rexpr; BPPARAM <- if ($nworkers == 1) NULL else MulticoreParam($nworkers)"
	Rexpr="$Rexpr; system.time(extract_snvs_from_RefSNP_json('$jsonfile', '$dump_subdir', $chunksize, BPPARAM))"
	$Rscript -e "$Rexpr" >$logfile 2>&1
}

# Takes about 100 hours on rex3 (80 logical cpus):
start_extraction chr1  9 15000 &
start_extraction chr2  9 15000 &
start_extraction chr3  7 12500 &
start_extraction chr4  6 12000 &
start_extraction chr5  5 11000 &
start_extraction chr6  5 11000 &
start_extraction chr7  5 11000 &
start_extraction chr8  5 11000 &
start_extraction chr9  1  5000 &
start_extraction chr10 1  5000 &
start_extraction chr11 1  5000 &
start_extraction chr12 1  5000 &
start_extraction chr13 1  5000 &
start_extraction chr14 1  5000 &
start_extraction chr15 1  5000 &
start_extraction chr16 1  5000 &
start_extraction chr17 1  5000 &
start_extraction chr18 1  5000 &
start_extraction chr19 1  5000 &
start_extraction chr20 1  5000 &
start_extraction chr21 1  5000 &
start_extraction chr22 1  5000 &
start_extraction chrX  1  5000 &
start_extraction chrY  1  5000 &
start_extraction chrMT 1  5000 &

