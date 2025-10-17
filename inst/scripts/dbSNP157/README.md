To build the **SNPlocs.Hsapiens.dbSNP157.GRCh38** and **SNPlocs.Hsapiens.dbSNP157.GRCh37** packages, proceed as follow:

1. Run script `download_json_files.sh`.

2. Run script `extract_snvs_from_RefSNP_json_files.sh`.

3. Run scripts `select_GRCh38_snvs.sh` and `select_GRCh37_snvs.sh` (can be run simultaneously).

4. Run scripts `build_GRCh38_OnDiskLongTable.sh` and `build_GRCh37_OnDiskLongTable.sh` (can be run simultaneously).

5. Edit files `DESCRIPTION` and `man/package.Rd` in newly forged packages by replacing placeholders `@TOTAL_SNPS@` (found in `DESCRIPTION` _and_ `man/package.Rd`) and `@NB_OF_SNP_LOCI_ON_CHR22@` (found in `man/package.Rd`) with corresponding values.

See comment at the beginning of each shell script for the details of how to run them.

