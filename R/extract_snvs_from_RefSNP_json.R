### =========================================================================
### extract_snvs_from_RefSNP_json()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Open a compress file
###

### Vectorized!
.has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x))
    stopifnot(isSingleString(suffix), !is.na(suffix))
    nc <- nchar(x)
    substr(x, nc - nchar(suffix) + 1L, nc) == suffix
}

.open_local_file <- function(filepath, open="rb")
{
    if (!isSingleString(filepath))
        stop(wmsg("path to local file must be a single string"))
    if (.has_suffix(filepath, ".gz")) {
        con <- gzfile(filepath, open=open)
    } else if (.has_suffix(filepath, ".bz2")) {
        con <- bzfile(filepath, open=open)
    } else if (.has_suffix(filepath, ".xz")) {
        con <- xzfile(filepath, open=open)
    } else {
        con <- file(filepath, open=open)
    }
    con
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .json_apply()
###

.json_apply <- function(json_lines, FUN, ..., BPPARAM=NULL)
{
    FUN_WRAPPER <- function(json_line, FUN, ...) {
        ## rjson::fromJSON() and jsonlite::parse_json() are both fast.
        ## Beware of a caveat with the latter: using 'simplifyVector=TRUE'
        ## makes it very slow and do weird things!
        ## Very slow: it makes jsonlite::parse_json() about 20x slower
        ## on 'refsnp-chrMT.json'! This is because the simplification is
        ## implemented in pure R (via jsonlite:::simplify()).
        ## Weird things: when parsing a RefSNP JSON file (e.g.
        ## 'refsnp-chrMT.json') the "present_obs_movements" field gets
        ## transformed in a weird way that makes it hard to work with.
        snp <- rjson::fromJSON(json_line)
        FUN(snp, ...)
    }
    if (is.null(BPPARAM)) {
        lapply(json_lines, FUN_WRAPPER, FUN, ...)
    } else {
        json_lines <- as.list(json_lines)
        bplapply(json_lines, FUN_WRAPPER, FUN, ..., BPPARAM=BPPARAM)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers for parsing a RefSNP JSON files
###
### RefSNP JSON (partly) documented at:
###   https://github.com/USF-HII/snptk/wiki/dbSNPJson
###
### We are only interested in the following parts of a JSON record:
###   { -- named list of length 10 (dict) --
###     "refsnp_id": string
###     "present_obs_movements": array (1 elt per movement)
###       [
###         { -- named list of length 6 (dict) --
###           "allele_in_cur_release":
###             { -- SPDI object (named list of length 4) --
###               "seq_id": string
###               "position": integer
###               "deleted_sequence": string
###               "inserted_sequence": string
###             }
###         }
###       ]
###     "primary_snapshot_data":
###       { -- named list of length 6 (dict) --
###         "variant_type": string
###         "placements_with_allele": array (1 elt per placement)
###           [
###             { -- named list of length 4 (dict) --
###               "seq_id": string
###               "is_ptlp": bool
###               "placement_annot":
###                 { -- named list of length 5 (dict) --
###                   "seq_type": string
###                 }
###               "alleles": array (1 elt per allele)
###                 [
###                   { -- named list of length 2 (dict) --
###                     "allele":
###                       { -- named list of length 1 (dict) --
###                         "spdi":
###                           { -- SPDI object (named list of length 4) --
###                             "seq_id": string
###                             "position": integer
###                             "deleted_sequence": string
###                             "inserted_sequence": string
###                           }
###                       }
###                   }
###                 ]
###             }
###           ]
###       }
###   }

.EXPECTED_SPDI_FIELDS <- c(
    "seq_id",
    "position",
    "deleted_sequence",
    "inserted_sequence"
)

### Returns a 1x4 data frame.
.make_spdi_record <- function(spdi)
{
    stopifnot(is.list(spdi), identical(names(spdi), .EXPECTED_SPDI_FIELDS))
    as.data.frame(spdi)
}

### Returns a list of 1x4 data frames.
.extract_spdi_records_from_alleles <- function(alleles)
{
    spdi_records <- lapply(alleles,
        function(allele) {
            stopifnot(is.list(allele),
                      identical(names(allele), c("allele", "hgvs")))
            allele <- allele$allele
            stopifnot(is.list(allele))
            spdi <- allele$spdi
            if (is.null(spdi))
                return(NULL)
            .make_spdi_record(spdi)
        })
    spdi_records <- S4Vectors:::delete_NULLs(spdi_records)
    stopifnot(length(spdi_records) != 0L)
    spdi_records
}

### 'spdi_records' must be a list of 1x4 data frames as returned
### by .extract_spdi_records_from_alleles() above.
### Returns a 4-component **list** (NOT a data.frame).
.collapse_spdi_records <- function(spdi_records)
{
    stopifnot(is.list(spdi_records))
    df <- do.call(rbind, spdi_records)
    stopifnot(identical(colnames(df), .EXPECTED_SPDI_FIELDS))
    seq_id <- unique(df[ , "seq_id"])
    stopifnot(length(seq_id) == 1L)
    position <- unique(as.integer(df[ , "position"]))
    if (length(position) != 1L)
        print(spdi_records)
    stopifnot(length(position) == 1L)
    deleted_sequence <- unique(df[ , "deleted_sequence"])
    if (length(deleted_sequence) != 1L)
        print(spdi_records)
    stopifnot(length(deleted_sequence) == 1L)
    inserted_sequences <- unique(df[ , "inserted_sequence"])
    list(seq_id=seq_id,
         position=position,
         deleted_sequence=deleted_sequence,
         inserted_sequences=inserted_sequences)
}

.check_placement <- function(placement)
{
    expected_fields <- c("seq_id", "is_ptlp", "placement_annot", "alleles")
    stopifnot(is.list(placement), identical(names(placement), expected_fields))
}

.extract_placement_alleles <- function(placement)
{
    alleles <- placement$alleles
    stopifnot(is.list(alleles))
    alleles
}

.extract_primary_snapshot_data <- function(snp)
{
    EXPECTED_TOP_FIELDS <- c(
        "refsnp_id", "create_date", "last_update_date", "last_update_build_id",
        "dbsnp1_merges", "citations", "lost_obs_movements",
        "present_obs_movements", "primary_snapshot_data", "mane_select_ids"
    )
    stopifnot(is.list(snp), identical(names(snp), EXPECTED_TOP_FIELDS))

    ## Get variant type from "primary_snapshot_data".
    data <- snp$primary_snapshot_data
    EXPECTED_DATA_FIELDS <- c(
        "placements_with_allele", "allele_annotations", "support", "anchor",
        "variant_type", "ga4gh"
    )
    stopifnot(is.list(data), identical(names(data), EXPECTED_DATA_FIELDS))

    data
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### quick_preview_RefSNP_json()
###

### Returns a 1x6 data frame.
.summarize_placement <- function(placement)
{
    is_ptlp <- placement$is_ptlp
    stopifnot(isTRUEorFALSE(is_ptlp))
    alleles <- .extract_placement_alleles(placement)
    seq_type <- placement$placement_annot$seq_type
    stopifnot(isSingleString(seq_type))
    spdi_records <- .extract_spdi_records_from_alleles(alleles)
    alleles <- .collapse_spdi_records(spdi_records)
    alleles$inserted_sequences <- paste(alleles$inserted_sequences,
                                        collapse=",")
    summarized_placement <- c(list(is_ptlp=is_ptlp, seq_type=seq_type), alleles)
    as.data.frame(summarized_placement)
}

### Returns a 3-component list or a NULL (if the SNP type is not one
### of the types specified via 'variant_type').
.summarize_snp <- function(snp, variant_type=NULL, seq_type=NULL)
{
    data <- .extract_primary_snapshot_data(snp)
    variant_type0 <- data$variant_type
    stopifnot(isSingleString(variant_type0))
    if (!(is.null(variant_type) || variant_type0 %in% variant_type))
        return(NULL)

    summarized_placements <- lapply(data$placements_with_allele,
        function(placement) {
            .check_placement(placement)
            if (!is.null(seq_type)) {
                seq_type0 <- placement$placement_annot$seq_type
                if (!(seq_type0 %in% seq_type))
                    return(NULL)
            }
            .summarize_placement(placement)
        })
    summarized_placements <- S4Vectors:::delete_NULLs(summarized_placements)

    if (length(summarized_placements) != 0L) {
        ## Turn 'summarized_placements' into a 5-col data frame with 1 row
        ## per placement.
        summarized_placements <- do.call(rbind, summarized_placements)
    }

    list(refsnp_id=snp$refsnp_id,
         variant_type=variant_type0,
         placements=summarized_placements)
}

### Returns a named list with 1 list element per SNP.
quick_preview_RefSNP_json <-
    function(con, variant_type=NULL, seq_type=NULL, n=6)
{
    if (is.character(con)) {
        if (!isSingleString(con))
            stop(wmsg("'con' must be a single string or a connection"))
        con <- .open_local_file(con)
        on.exit(close(con))
    }
    if (!(is.null(variant_type) || is.character(variant_type)))
        stop(wmsg("'variant_type' must be NULL or a character vector"))
    if (!(is.null(seq_type) || is.character(seq_type)))
        stop(wmsg("'seq_type' must be NULL or a character vector"))
    json_lines <- readLines(con, n=n)
    .json_apply(json_lines,
                .summarize_snp, variant_type=variant_type, seq_type=seq_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_snvs_from_RefSNP_json()
###

### Returns a 6-col data frame or a NULL.
.summarize_snv <- function(snp)
{
    summarized_snp <- .summarize_snp(snp, variant_type="snv",
                                          seq_type="refseq_chromosome")
    if (is.null(summarized_snp))
        return(NULL)
    placements <- summarized_snp$placements
    stopifnot(is.data.frame(placements))
    keep_columns <- c("is_ptlp", "seq_id", "position",
                      "deleted_sequence", "inserted_sequences")
    cbind(refsnp_id=summarized_snp$refsnp_id,
          placements[ , keep_columns, drop=FALSE])
}

### Returns a SplitDataFrameList object or a NULL.
.group_summarized_snvs_by_seq_id <- function(summarized_snvs)
{
    stopifnot(is.list(summarized_snvs))
    DF <- as(do.call(rbind, summarized_snvs), "DataFrame")
    if (nrow(DF) == 0L)
        return(NULL)
    j <- match("seq_id", colnames(DF))
    stopifnot(!is.na(j))
    f <- DF[ , j]
    split(DF[ , -j, drop=FALSE], f)
}

.dump_summarized_snvs <- function(summarized_snvs, dump_dir)
{
    snvs_per_seq_id <- .group_summarized_snvs_by_seq_id(summarized_snvs)
    nb_snvs <- if (is.null(snvs_per_seq_id)) 0L
               else nrow(unlist(snvs_per_seq_id, use.names=FALSE))
    cat("  --> ", nb_snvs, " snvs extracted\n", sep="")
    if (nb_snvs == 0L)
        return(NULL)
    for (i in seq_along(snvs_per_seq_id)) {
        seq_id <- names(snvs_per_seq_id)[[i]]
        out_file <- file.path(dump_dir, paste0(seq_id, ".tab"))
        summarized_snvs <- as.data.frame(snvs_per_seq_id[[i]])
        cat("    - writing ", nrow(summarized_snvs),
            " snvs to ", out_file, " ... ", sep="")
        write.table(summarized_snvs, file=out_file, append=TRUE,
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        cat("ok\n")
    }
    lengths(snvs_per_seq_id)
}

### Optimal chunksize seems to somehwat depend on number of workers. Some
### timings for processing a RefSNP json file with 150000 lines (1 snp/line):
###
###   nb of workers   optimal chunksize   total time in sec.
###   -------------   -----------------   -----------------
###               8               10000             390.778
###
###              10                2500             420.672
###              10                5000             377.874
###              10               10000             373.494
###              10               15000             370.211
###              10               20000             370.749
###              10               25000             393.568
###
###              16                5000
###              16               10000             305.232
###              16               15000             330.971
###              16               20000             275.598
###              16               25000             348.339
###              16               30000             346.602
###
###              32               40000             269.055
###
extract_snvs_from_RefSNP_json <- function(con, dump_dir,
                                          chunksize=10000, BPPARAM=NULL)
{
    if (is.character(con)) {
        if (!isSingleString(con))
            stop(wmsg("'con' must be a single string or a connection"))
        con <- .open_local_file(con)
        on.exit(close(con))
    }

    if (!isSingleString(dump_dir))
        stop(wmsg("'dump_dir' must be a single string specifying the path ",
                  "to the directory where to dump the snvs"))
    if (!dir.exists(dump_dir))
        stop(wmsg("'dump_dir' must be the path to an existing directory"))

    if (!isSingleNumber(chunksize))
        stop(wmsg("'chunksize' must be a single integer"))
    if (!is.integer(chunksize))
        chunksize <- as.integer(chunksize)

    if (!(is.null(BPPARAM) || is(BPPARAM, "BiocParallelParam")))
        stop(wmsg("'BPPARAM' must be NULL or a BiocParallelParam object"))

    offset <- 0L
    while (TRUE) {
        if (chunksize >= 1L)
            cat("Reading lines (", chunksize, " max) ... ", sep="")
        json_lines <- readLines(con, n=chunksize)
        nline <- length(json_lines)
        if (nline == 0L) {
            if (chunksize >= 1L)
                cat("no more lines to read!\n")
            break
        }
        if (chunksize >= 1L) {
            from <- offset + 1L
            to <- offset + nline
            if (is.null(BPPARAM)) {
                workers <- ""
            } else {
                workers <- paste0(" (using ", bpnworkers(BPPARAM), " workers)")
            }
            cat("ok; processing lines ", from, "-", to,
                workers, " ... ", sep="")
        }
        st <- system.time(
            summarized_snvs <-.json_apply(json_lines, .summarize_snv,
                                          BPPARAM=BPPARAM)
        )
        if (chunksize >= 1L)
            cat("ok (in ", st[["elapsed"]], " sec.)\n", sep="")
        .dump_summarized_snvs(summarized_snvs, dump_dir)
        if (chunksize >= 1L)
            cat("\n")
        offset <- offset + nline
        if (chunksize >= 1L && nline < chunksize)
            break
    }
    cat("DONE.\n")
    invisible(offset)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OBSOLETE - extract_raw_snps_from_RefSNP_json()
###

### Returns NULL if not a preferred top level placement (PTLP).
.extract_alleles_if_placement_is_ptlp <- function(placement)
{
    .check_placement(placement)
    if (!placement$is_ptlp)
        return(NULL)
    alleles <- .extract_placement_alleles(placement)
    spdi_records <- .extract_spdi_records_from_alleles(alleles)
    ptlp_alleles <- .collapse_spdi_records(spdi_records)
    stopifnot(identical(ptlp_alleles$seq_id, placement$seq_id))
    ptlp_alleles
}

### Assumes that 'placements' contains exactly **one** preferred top level
### placement (PTLP) and returns the alleles corresponding to that placement.
.extract_ptlp_alleles <- function(placements)
{
    stopifnot(is.list(placements))
    alleles <- lapply(placements, .extract_alleles_if_placement_is_ptlp)
    alleles <- S4Vectors:::delete_NULLs(alleles)
    stopifnot(length(alleles) == 1L)
    alleles[[1L]]
}

### Returns NULL if 'movements' is empty.
.extract_alleles_from_movements <- function(movements)
{
    stopifnot(is.list(movements))
    if (length(movements) == 0L)
        return(NULL)
    spdi_records <- lapply(movements,
        function(mov) {
            stopifnot(is.list(mov))
            expected_mov_fields <- c(
                "component_ids", "observation", "allele_in_cur_release",
                "other_rsids_in_cur_release", "previous_release",
                "last_added_to_this_rs"
            )
            if (length(mov) == 5L)
                expected_mov_fields <- expected_mov_fields[-5L]
            stopifnot(identical(names(mov), expected_mov_fields))
            .make_spdi_record(mov$allele_in_cur_release)
        })
    .collapse_spdi_records(spdi_records)
}

### Returns a 1x6 data frame.
.extract_raw_snp <- function(snp, chrominfo, paranoid=FALSE)
{
    if (!isTRUEorFALSE(paranoid))
        stop(wmsg("'paranoid' must be TRUE or FALSE"))
    stopifnot(is.data.frame(chrominfo))

    data <- .extract_primary_snapshot_data(snp)
    variant_type <- data$variant_type
    stopifnot(isSingleString(variant_type))

    ## Extract alleles from the preferred top level placement (PTLP).
    ptlp_alleles <- .extract_ptlp_alleles(data$placements_with_allele)

    if (paranoid) {
        ## Extract alleles from "present_obs_movements".
        alleles2 <- .extract_alleles_from_movements(snp$present_obs_movements)
        if (!is.null(alleles2)) {
            stopifnot(identical(ptlp_alleles$seq_id, alleles2$seq_id))
            stopifnot(identical(ptlp_alleles$position, alleles2$position))
            stopifnot(identical(ptlp_alleles$deleted_sequence,
                                ptlp_alleles$deleted_sequence))
            stopifnot(all(alleles2$inserted_sequences %in%
                          ptlp_alleles$inserted_sequences))
        }
    }

    ## Replace 'seq_id' with official sequence name from NCBI assembly.
    seq_id <- ptlp_alleles$seq_id
    m <- match(seq_id, chrominfo[ , "RefSeqAccn"])
    if (is.na(m)) {
        seqname <- seq_id
    } else {
        seqname <- chrominfo[m , "SequenceName"]
    }

    inserted_sequences <- paste(ptlp_alleles$inserted_sequences, collapse=",")
    data.frame(refsnp_id=snp$refsnp_id,
               variant_type=variant_type,
               seqname=seqname,
               position=ptlp_alleles$position,
               deleted_sequence=ptlp_alleles$deleted_sequence,
               inserted_sequences=inserted_sequences)
}

### Returns a 6-col data frame with 1 row per JSON line.
.extract_raw_snps_from_json_lines <- function(json_lines, chrominfo,
                                              paranoid=FALSE)
{
    stopifnot(is.character(json_lines))
    raw_snps <- .json_apply(json_lines,
        function(snp) {
            raw_snp <- try(.extract_raw_snp(snp, chrominfo, paranoid=paranoid))
            if (inherits(raw_snp, "try-error"))
                stop(wmsg("Failed to extract raw snp from line ", snp))
            raw_snp
        })
    do.call(rbind, raw_snps)
}

### Using 'paranoid=TRUE' performs additional sanity checks on the extracted
### alleles but at the cost of a 50-60% slowdown!
### Typical use:
###   json_file <- "refsnp-chrMT.json"
###   out_file <- "chrMT_raw_snps.tab"
###   extract_raw_snps_from_RefSNP_json(json_file, out_file, chunksize=1000)
extract_raw_snps_from_RefSNP_json <- function(con, out="", chunksize=50000,
                                              assembly="GRCh38.p13",
                                              paranoid=FALSE)
{
    if (!isSingleNumber(chunksize))
        stop(wmsg("'chunksize' must be a single integer"))
    if (!is.integer(chunksize))
        chunksize <- as.integer(chunksize)
    if (is.character(con)) {
        if (!isSingleString(con))
            stop(wmsg("'con' must be a single string or a connection"))
        con <- .open_local_file(con)
        on.exit(close(con))
    }
    if (is.character(out)) {
        if (!isSingleString(out))
            stop(wmsg("'out' must be a single string or a connection"))
        if (out != "") {
            out <- file(out, "w")
            on.exit(close(out), add=TRUE)
        }
    }
    chrominfo <- getChromInfoFromNCBI(assembly)
    offset <- 0L
    while (TRUE) {
        if (chunksize >= 1L)
            cat("Reading lines (", chunksize, " max) ... ", sep="")
        json_lines <- readLines(con, n=chunksize)
        nline <- length(json_lines)
        if (nline == 0L) {
            if (chunksize >= 1L)
                cat("no more lines to read!\n")
            break
        }
        if (chunksize >= 1L) {
            from <- offset + 1L
            to <- offset + nline
            cat("ok; processing lines ", from, "-", to, " ... ", sep="")
        }
        raw_snps <- .extract_raw_snps_from_json_lines(json_lines, chrominfo,
                                                      paranoid=paranoid)
        if (chunksize >= 1L)
            cat("ok; writing ", nrow(raw_snps), " raw snps ... ", sep="")
        write.table(raw_snps, file=out, append=TRUE,
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        if (chunksize >= 1L)
            cat("ok\n")
        offset <- offset + nline
        if (chunksize >= 1L && nline < chunksize)
            break
    }
    cat("DONE.\n")
    invisible(offset)
}

