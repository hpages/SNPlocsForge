### =========================================================================
### select_snvs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### I/O helpers
###

.load_snvs <- function(file)
{
    COL2CLASS <- c(rsid="integer",
                   is_ptlp="logical",
                   pos0="integer",
                   deleted_sequence="character",
                   inserted_sequences="character")
    read.table(file, sep="\t", quote="", col.names=names(COL2CLASS),
               na.strings="?", colClasses=unname(COL2CLASS),
               stringsAsFactors=FALSE)
}

### Used by .load_selected_snvs_from_multiple_files() defined in
### file build_OnDiskLongTable.R
load_snvs_from_multiple_files <- function(filepaths)
{
    lapply(setNames(filepaths, basename(filepaths)),
        function(filepath) {
            cat("- loading SNVs from ", filepath, " ... ", sep="")
            snvs <- .load_snvs(filepath)
            cat("ok [", nrow(snvs), " SNVs loaded]\n", sep="")
            snvs
        })
}

.save_snvs <- function(snvs, file)
{
    write.table(snvs, file=file,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

.save_selected_snvs_to_multiple_files <- function(selected_snv_dfs, out_dir,
                                                  original_snv_dfs)
{
    stopifnot(is.list(selected_snv_dfs), is.list(original_snv_dfs))
    filenames <- names(selected_snv_dfs)
    stopifnot(identical(filenames, names(original_snv_dfs)))

    for (i in seq_along(selected_snv_dfs)) {
        snvs <- selected_snv_dfs[[i]]
        num_snvs <- nrow(snvs)
        out_file <- file.path(out_dir, filenames[[i]])
        if (num_snvs == 0L) {
            cat("- nothing to save to ", out_file, "\n", sep="")
        } else {
            cat("- saving ", num_snvs, "/", nrow(original_snv_dfs[[i]]),
                " SNVs to ", out_file, " ... ", sep="")
            .save_snvs(snvs, out_file)
            cat("ok\n")
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### select_snvs()
###

### Select the snv files that belong to the specified assembly.
.select_snv_files <- function(dump_dir, assembly)
{
    stopifnot(isSingleString(dump_dir))

    chrominfo <- getChromInfoFromNCBI(assembly, assembled.molecules.only=TRUE)
    target_files <- paste0(chrominfo[ , "RefSeqAccn"], ".tab")

    snv_files <- dir(dump_dir, pattern="\\.tab$", full.names=TRUE)
    snv_files[basename(snv_files) %in% target_files]
}

### Must be run after extract_snvs_from_RefSNP_json().
### Usage:
###   dump_dir <- "~/SNPlocsForge/downloads/dbSNP155/snvs_dump/chr22"
###   out_dir  <- "~/SNPlocsForge/downloads/dbSNP155/GRCh37_snvs/chr22"
###   select_snvs(dump_dir, out_dir, assembly="GRCh37.p13")
select_snvs <- function(dump_dir, out_dir, assembly="GRCh38.p13")
{
    stopifnot(isSingleString(out_dir))

    snv_files <- .select_snv_files(dump_dir, assembly)

    original_snv_dfs <- load_snvs_from_multiple_files(snv_files)
    nrows <- vapply(original_snv_dfs, nrow, integer(1))
    original_snv_dfs <- original_snv_dfs[order(nrows, decreasing=TRUE)]

    cat("- keeping one placement per RefSNP id ... ")
    all_rsids <- IntegerList(lapply(original_snv_dfs,
                                    function(snvs) snvs[ , "rsid"]))
    rsids <- unlist(all_rsids, use.names=FALSE)
    keep_idx <- relist(!duplicated(rsids), all_rsids)
    selected_snv_dfs <- mapply(function(snvs, idx) snvs[idx, , drop=FALSE],
                               original_snv_dfs, keep_idx, SIMPLIFY=FALSE)
    cat("ok\n")

    .save_selected_snvs_to_multiple_files(selected_snv_dfs, out_dir,
                                          original_snv_dfs)
}

