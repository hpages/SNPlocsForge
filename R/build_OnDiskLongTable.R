###

.load_raw_snps <- function(filepath)
{
    COL2CLASS <- c(rsid="integer",
                   variant_type="character",
                   seqname="character",
                   pos0="integer",
                   deleted_sequence="character",
                   inserted_sequences="character")
    read.table(filepath, sep="\t", quote="", col.names=names(COL2CLASS),
               na.strings="?", colClasses=unname(COL2CLASS),
               stringsAsFactors=FALSE)
}

.bad_snps <- function(raw_snps, bad_idx)
{
    bad_snps <- raw_snps[bad_idx, "rsid"]
    nbad <- length(bad_snps)
    if (nbad > 6L)
        bad_snps <- c(head(bad_snps, n=5L), "...")
    bad_snps <- paste(bad_snps, collapse=", ")
    if (nbad > 6L)
        bad_snps <- paste0(bad_snps, " [", nbad - 5L, " more]")
    bad_snps
}

### Return SNPs of type "snv" in a 3-col data.frame:
###   1. rsid: integer vector (RefSNP id without "rs" prefix)
###   2. pos: integer vector (one-based position)
###   3. alleles: raw vector (alleles as an IUPAC letter turned into
###      byte value).
.cook_raw_snps <- function(raw_snps, expected_seqname)
{
    stopifnot(isSingleString(expected_seqname))

    ## Keep SNPs located on expected seqname only.
    seqname <- raw_snps[ , "seqname"]
    keep_idx <- which(seqname == expected_seqname)
    raw_snps <- raw_snps[keep_idx, , drop=FALSE]

    ## Keep SNPs of type "snv" only.
    type <- raw_snps[ , "variant_type"]
    keep_idx <- which(type == "snv")
    raw_snps <- raw_snps[keep_idx, , drop=FALSE]

    rsid <- raw_snps[ , "rsid"]
    stopifnot(is.integer(rsid))

    pos <- raw_snps[ , "pos0"] + 1L

    alleles <- gsub("/", "", raw_snps[ , "inserted_sequences"], fixed=TRUE)
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
.drop_out_of_bounds_snps <- function(cooked_snps, seqlength)
{
    pos <- cooked_snps[ , "pos"]
    is_out_of_bounds <- pos < 1L | pos > seqlength
    nb_out_of_bounds <- sum(is_out_of_bounds)
    if (nb_out_of_bounds != 0L) {
        cat("  DROP ", nb_out_of_bounds, " OUT OF BOUNDS SNPS! ... ", sep="")
        keep_idx <- which(!is_out_of_bounds)
        cooked_snps <- cooked_snps[keep_idx, , drop=FALSE]
        cat("OK\n")
    }
    cooked_snps
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
build_OnDiskLongTable <- function(tmp_dir, seqnames, chr_prefix="chr",
                                  batchsize=200000L)
{
    cat("\n")
    cat("***************** START build_OnDiskLongTable() ******************\n")

    seqnames <- strsplit(seqnames, " ", fixed=TRUE)[[1L]]

    seqinfo <- BSgenome:::read_seqinfo_table("seqinfo.txt", genome="GRCh38.p7")

    rowids <- vector("list", length=length(seqnames))
    names(rowids) <- seqnames

    append <- FALSE
    for (seqname in seqnames) {
        cat("\n")
        cat("Processing SNPs for chromosome ", seqname, ":\n", sep="")

        filename <- paste0(chr_prefix, seqname, "_raw_snps.tab")
        cat("  Loading raw SNPs from ", filename, " ... ", sep="")
        filepath <- file.path(tmp_dir, filename)
        raw_snps <- .load_raw_snps(filepath)
        cat("OK [", nrow(raw_snps), " raw SNPs loaded]\n", sep="")

        cat("  Cooking raw SNPs ... ", sep="")
        cooked_snps <- .cook_raw_snps(raw_snps, seqname)
        cat("OK\n")

        seqlength <- seqlengths(seqinfo)[[seqname]]
        cooked_snps <- .drop_out_of_bounds_snps(cooked_snps, seqlength)

        rowids[[seqname]] <- cooked_snps[ , 1L]

        df <- cooked_snps[ , -1L, drop=FALSE]
        spatial_index <- .build_spatial_index(df[ , "pos"], batchsize,
                                              seqname, seqinfo)
        fmt <- paste0("%s ", nrow(df), " cooked SNPs ",
                      "%s OnDiskLongTable directory structure")
        if (!append) {
            msg <- sprintf(fmt, "Writing", "as")
        } else {
            msg <- sprintf(fmt, "Appending", "to")
        }
        cat("  ", msg, " ... ", sep="")
        writeOnDiskLongTable(df, spatial_index=spatial_index,
                                 append=append)
        append <- TRUE
        cat("OK\n")
    }

    cat("\n")

    rowids <- unlist(rowids, recursive=FALSE, use.names=FALSE)
    cat("Adding RefSNP ids to OnDiskLongTable directory structure ... ")
    ## Using compress="xz" reduces the size on disk by < 2% but makes further
    ## loading of the row ids (with readRDS()) 8x slower. Not worth it!
    #writeOnDiskLongTableRowids(rowids, compress="xz")
    writeOnDiskLongTableRowids(rowids)
    cat("OK\n")

    cat("\n")
    cat("****************** END build_OnDiskLongTable() *******************\n")
    cat("Total number of SNPs written to disk: ", length(rowids), "\n", sep="")
    cat("\n")
}

