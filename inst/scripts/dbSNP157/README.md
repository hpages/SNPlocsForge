To forge the SNPlocs.Hsapiens.dbSNP157.GRCh38 and SNPlocs.Hsapiens.dbSNP157.GRCh37 packages, run the following scripts in this order:

1. `download_json_files.sh`

2. `extract_snvs_from_RefSNP_json_files.sh`

3. `select_GRCh38_snvs.sh` and `select_GRCh37_snvs.sh` (can be run in any order or simultaneously)

4. `build_GRCh38_OnDiskLongTable.sh` and `build_GRCh37_OnDiskLongTable.sh` (can be run in any order or simultaneously)

See comment at the beginning of the scripts for the details of how to run each of them.

