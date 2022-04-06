### RefSNP JSON (partly) documented at:
###   https://github.com/USF-HII/snptk/wiki/dbSNPJson

.EXPECTED_SPDI_FIELDS <- c(
    "seq_id",
    "position",
    "deleted_sequence",
    "inserted_sequence"
)

### 'spdi_records' must be a list of 1x4 data frames.
### Returns a 4-component **list** (NOT data.frame).
.collapse_spdi_records <- function(spdi_records)
{
    stopifnot(is.list(spdi_records))
    df <- do.call(rbind, spdi_records)
    stopifnot(identical(colnames(df), .EXPECTED_SPDI_FIELDS))
    seq_id <- unique(df[ , "seq_id"])
    stopifnot(length(seq_id) == 1L)
    position <- unique(df[ , "position"])
    stopifnot(length(position) == 1L)
    deleted_sequence <- unique(df[ , "deleted_sequence"])
    stopifnot(length(deleted_sequence) == 1L)
    inserted_sequences <- unique(df[ , "inserted_sequence"])
    list(seq_id=seq_id, position=position,
         deleted_sequence=deleted_sequence,
         inserted_sequences=inserted_sequences)
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
            spdi <- mov$allele_in_cur_release
            stopifnot(is.list(spdi),
                      identical(names(spdi), .EXPECTED_SPDI_FIELDS))
            as.data.frame(spdi)
        })
    .collapse_spdi_records(spdi_records)
}

### Returns NULL if not a preferred top level placement (PTLP).
.extract_alleles_from_placement <- function(placement)
{
    expected_fields <- c("seq_id", "is_ptlp", "placement_annot", "alleles")
    stopifnot(is.list(placement), identical(names(placement), expected_fields))
    if (!placement$is_ptlp)
        return(NULL)
    stopifnot(is.list(placement$alleles))
    spdi_records <- lapply(placement$alleles,
        function(allele) {
            stopifnot(is.list(allele),
                      identical(names(allele), c("allele", "hgvs")))
            allele <- allele$allele
            stopifnot(is.list(allele), identical(names(allele), "spdi"))
            spdi <- allele$spdi
            stopifnot(is.list(spdi),
                      identical(names(spdi), .EXPECTED_SPDI_FIELDS))
            as.data.frame(spdi)
        })
    alleles <- .collapse_spdi_records(spdi_records)
    stopifnot(identical(alleles$seq_id, placement$seq_id))
    alleles
}

### Assumes that 'placements' contains exactly **one** preferred top level
### placement (PTLP) and returns the alleles corresponding to that placement.
.extract_alleles_from_placements <- function(placements)
{
    stopifnot(is.list(placements))
    alleles <- lapply(placements, .extract_alleles_from_placement)
    alleles <- S4Vectors:::delete_NULLs(alleles)
    stopifnot(length(alleles) == 1L)
    alleles[[1L]]
}

.extract_raw_snp <- function(snp, chrominfo, paranoid=FALSE)
{
    if (!isTRUEorFALSE(paranoid))
        stop(wmsg("'paranoid' must be TRUE or FALSE"))
    stopifnot(is.data.frame(chrominfo))
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
    variant_type <- data$variant_type
    stopifnot(isSingleString(variant_type))

    ## Extract alleles from "placements_with_allele".
    alleles <- .extract_alleles_from_placements(data$placements_with_allele)

    if (paranoid) {
        ## Extract alleles from "present_obs_movements".
        alleles2 <- .extract_alleles_from_movements(snp$present_obs_movements)
        if (!is.null(alleles2)) {
            stopifnot(identical(alleles$seq_id, alleles2$seq_id))
            stopifnot(identical(alleles$position, alleles2$position))
            stopifnot(identical(alleles$deleted_sequence,
                                alleles2$deleted_sequence))
            stopifnot(all(alleles2$inserted_sequences %in%
                          alleles$inserted_sequences))
        }
    }

    ## Replace 'seq_id' with official sequence name from NCBI assembly.
    seq_id <- alleles$seq_id
    m <- match(seq_id, chrominfo[ , "RefSeqAccn"])
    if (is.na(m)) {
        seqname <- seq_id
    } else {
        seqname <- chrominfo[m , "SequenceName"]
    }

    inserted_sequences <- paste(alleles$inserted_sequences, collapse="/")
    c(snp$refsnp_id, variant_type, seqname, alleles$position,
      alleles$deleted_sequence, inserted_sequences)
}

### Based on rjson::fromJSON().
### rjson::fromJSON() and jsonlite::parse_json() are both fast. However,
### there's a caveat with the latter: using 'simplifyVector=TRUE' makes it
### very slow and do weird things!
### Very slow: it makes jsonlite::parse_json() about 20x slower
### on 'refsnp-chrMT.json'! This is because the simplification is
### implemented in pure R (via jsonlite:::simplify()).
### Weird things: when parsing 'refsnp-chrMT.json' the "present_obs_movements"
### field gets transformed in a weird way that makes it hard to work with.
### Using 'paranoid=TRUE' performs additional sanity checks on the extracted
### alleles but at the cost of a 50-60% slowdown!
### Typical use:
###   json_file <- "refsnp-chrMT.json"
###   out_file <- "chrMT_raw_snps.tab"
###   extract_raw_snps_from_json(json_file, out_file, n=1000)
extract_raw_snps_from_json <- function(con, out="", n=50000,
                                       assembly="GRCh38.p13", paranoid=FALSE)
{
    if (!isSingleNumber(n))
        stop(wmsg("'n' must be a single integer"))
    if (!is.integer(n))
        n <- as.integer(n)
    if (is.character(con)) {
        con <- file(con, "rb")
        on.exit(close(con))
    }
    if (is.character(out) && out != "") {
        out <- file(out, "w")
        on.exit(close(out), add=TRUE)
    }
    chrominfo <- getChromInfoFromNCBI(assembly)
    lineno <- 1L
    while (length(json_lines <- readLines(con, n=n)) != 0L) {
        if (n >= 1L) {
            to <- lineno + length(json_lines) - 1L
            cat("processing lines ", lineno, "-", to, " ... ", sep="")
        }
        for (json_line in json_lines) {
            snp <- rjson::fromJSON(json_line)
            raw_snp <- try(.extract_raw_snp(snp, chrominfo, paranoid=paranoid))
            if (inherits(raw_snp, "try-error"))
                stop(wmsg("Failed to extract raw snp from line ", json_line))
            cat(paste(raw_snp, collapse="\t"), "\n", sep="", file=out)
            lineno <- lineno + 1L
        }
        if (n >= 1L)
            cat("ok\n")
    }
    cat("DONE.\n")
    invisible(lineno - 1L)
}

