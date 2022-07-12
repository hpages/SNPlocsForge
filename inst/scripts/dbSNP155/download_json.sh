#!/bin/bash
# 
# To run this script in "batch mode":
#
#   ./download_json.sh >download_json.log 2>&1 &
#

set -e  # exit immediately if a simple command exits with a non-zero status

SOURCE_DATA_URL="ftp.ncbi.nih.gov/snp/archive/b155/JSON"  # do NOT include https:// or ftp:// here!

# This will download about 380 GB of data!
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT"
#CHROMOSOMES="chrMT chrY"
JSON_PREFIX="refsnp-"
JSON_SUFFIX=".json.bz2"

if [ -e .listing ]; then
        echo "The current directory already contains a .listing file, which was probably generated "
        echo "by a previous run of this script. Please remove it (and eventually clean the current "
        echo "directory) before you run this script again. Alternatively, run this script in "
        echo "another directory."
        echo ""
        exit -1
fi

echo "List of files available at $SOURCE_DATA_URL:"
curl --silent "ftp://$SOURCE_DATA_URL/" >.listing
cat .listing
echo ""

for chr in $CHROMOSOMES; do
        filename="${JSON_PREFIX}${chr}${JSON_SUFFIX}"
        url="https://$SOURCE_DATA_URL/$filename"  # https:// is faster than ftp://
        echo "Downloading $url ..."
        curl -O "$url"
        echo ""
done

echo "Downloading CHECKSUMS to check MD5 sums (we didn't download all the files listed in CHECKSUMS so expect some 'No such file or directory' errors below) ..."
rm CHECKSUMS
url="https://$SOURCE_DATA_URL/CHECKSUMS"
curl -O "$url"
echo ""
md5sum --check CHECKSUMS
echo ""

echo "DONE."
echo ""

echo "A listing of all the files available at $SOURCE_DATA_URL"
echo "at the time of the download was dumped in the .listing file."

