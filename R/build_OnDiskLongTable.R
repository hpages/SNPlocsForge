### =========================================================================
### build_OnDiskLongTable()
### -------------------------------------------------------------------------


.collect_snv_files <- function(seqname, chrominfo, all_snv_files)
{
    m <- match(seqname, chrominfo[ , "SequenceName"])
    if (is.na(m))
        stop(wmsg("unknown chromosome: ", seqname))
    fname <- paste0(chrominfo[m, "RefSeqAccn"], ".tab")
    all_snv_files[which(basename(all_snv_files) == fname)]
}

### Return SNVs in a 3-col data.frame:
###   1. rsid: integer vector (RefSNP id without "rs" prefix)
###   2. pos: integer vector (one-based position)
###   3. alleles: raw vector (alleles as an IUPAC letter turned into
###      byte value).
.cook_snvs <- function(snvs)
{
    rsid <- snvs[ , "rsid"]
    stopifnot(is.integer(rsid))

    pos <- snvs[ , "pos0"] + 1L

    alleles <- gsub(",", "", snvs[ , "inserted_sequences"], fixed=TRUE)
    alleles <- BSgenome:::encode_letters_as_bytes(mergeIUPACLetters(alleles))

    ans <- data.frame(rsid=rsid,
                      pos=pos,
                      alleles=alleles,
                      stringsAsFactors=FALSE)
    ans <- ans[order(ans[ , "pos"]), , drop=FALSE]
    row.names(ans) <- NULL
    ans
}

### Believe it or not, but some SNPs in dbSNP are reported to be at a position
### that is beyond the end of the chromosome. For example, rs553244808 (from
### dbSNP 151) was reported to be at position 143544518 on chromosome 14 in
### GRCh38.p7, even though the length of this chromosome is 107043718. We
### drop these SNPs.
.drop_out_of_bounds_snvs <- function(cooked_snvs, seqlength)
{
    pos <- cooked_snvs[ , "pos"]
    is_out_of_bounds <- pos < 1L | pos > seqlength
    nb_out_of_bounds <- sum(is_out_of_bounds)
    if (nb_out_of_bounds != 0L) {
        cat("  DROP ", nb_out_of_bounds, " OUT OF BOUNDS SNVS! ... ", sep="")
        keep_idx <- which(!is_out_of_bounds)
        cooked_snvs <- cooked_snvs[keep_idx, , drop=FALSE]
        cat("OK\n")
    }
    cooked_snvs
}

.build_spatial_index <- function(pos, batchsize, seqname, seqinfo)
{
    stopifnot(is.integer(pos))
    if (is.unsorted(pos))
        stop(wmsg("'pos' must be sorted"))
    chunks <- breakInChunks(length(pos), chunksize=batchsize)
    spatial_ranges <- range(relist(pos, chunks))
    GRanges(seqname, IRanges(spatial_ranges[ , 1L],
                             spatial_ranges[ , 2L]),
            batchsize=width(chunks),
            seqinfo=seqinfo)
}

### 'seqnames' must be a single string (e.g. "20 21 22")
build_OnDiskLongTable <- function(dump_dir, seqnames, assembly="GRCh38.p13",
                                  batchsize=200000L, rowids_nchunk=6L)
{
    stopifnot(isSingleString(dump_dir),
              isSingleString(seqnames),
              isSingleString(assembly),
              isSingleNumber(batchsize),
              isSingleNumber(rowids_nchunk),
              rowids_nchunk >= 1L,
              rowids_nchunk <= 9L)

    COLNAMES <- c("SequenceName", "SequenceLength", "circular",
                  "GenBankAccn", "RefSeqAccn")
    ci <- read.table("seqinfo.txt", col.names=COLNAMES, stringsAsFactors=FALSE)
    ci2 <- getChromInfoFromNCBI(assembly, assembled.molecules.only=TRUE)
    stopifnot(identical(ci, ci2[ , colnames(ci)]))
    seqinfo <- Seqinfo(ci[ , "SequenceName"], ci[ , "SequenceLength"],
                       ci[ , "circular"], assembly)

    all_snv_files <- dir(dump_dir, pattern="\\.tab$", recursive=TRUE,
                         full.names=TRUE)

    cat("\n")
    cat("***************** START build_OnDiskLongTable() ******************\n")

    seqnames <- strsplit(seqnames, " ", fixed=TRUE)[[1L]]
    rowids <- vector("list", length=length(seqnames))
    names(rowids) <- seqnames

    append <- FALSE
    for (seqname in seqnames) {
        cat("\n")
        cat("Processing SNVs for chromosome ", seqname, ":\n", sep="")

        filepaths <- .collect_snv_files(seqname, ci, all_snv_files)
        snvs <- do.call(rbind, load_snvs_from_multiple_files(filepaths))
        cat("- total number of SNVs on chromosome ", seqname, ": ",
            nrow(snvs), "\n", sep="")

        cat("- cooking the SNVs ... ", sep="")
        cooked_snvs <- .cook_snvs(snvs)
        cat("ok\n")

        seqlength <- seqlengths(seqinfo)[[seqname]]
        cooked_snvs <- .drop_out_of_bounds_snvs(cooked_snvs, seqlength)

        rowids[[seqname]] <- cooked_snvs[ , 1L]

        df <- cooked_snvs[ , -1L, drop=FALSE]
        spatial_index <- .build_spatial_index(df[ , "pos"], batchsize,
                                              seqname, seqinfo)
        fmt <- paste0("- %s ", nrow(df), " cooked SNVs ",
                      "%s OnDiskLongTable directory structure")
        if (!append) {
            msg <- sprintf(fmt, "writing", "as")
        } else {
            msg <- sprintf(fmt, "appending", "to")
        }
        cat(msg, " ... ", sep="")
        writeOnDiskLongTable(df, spatial_index=spatial_index,
                                 append=append)
        append <- TRUE
        cat("ok\n")
    }

    cat("\n")

    rowids <- unlist(rowids, recursive=FALSE, use.names=FALSE)

    cat("Adding RefSNP ids to OnDiskLongTable directory structure ... ")
    writeOnDiskLongTableRowids(rowids, nchunk=rowids_nchunk)
    cat("OK\n")

    cat("\n")
    cat("****************** END build_OnDiskLongTable() *******************\n")
    cat("Total number of SNVs written to disk: ", length(rowids), "\n", sep="")
    cat("\n")
}

