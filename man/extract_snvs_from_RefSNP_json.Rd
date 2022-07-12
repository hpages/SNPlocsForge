\name{extract_snvs_from_RefSNP_json}

\alias{extract_snvs_from_RefSNP_json}

\title{Extract SNPs of type "snv" from a RefSNP JSON file}

\description{
  Extract SNPs of type "snv" from a RefSNP JSON file.
}

\usage{
extract_snvs_from_RefSNP_json(con, dump_dir,
                              chunksize=10000, BPPARAM=NULL)
}

\arguments{
  \item{con}{
    File path or connection to a RefSNP JSON file (compressed files
    are supported).
  }
  \item{dump_dir}{
    Path to the directory where to dump the snvs.
  }
  \item{chunksize}{
    How many JSON lines to load at once in memory. Set to -1 to load the
    entire file in memory (strongly discouraged on a big JSON file!)
  }
  \item{BPPARAM}{
    \code{NULL} or a \link[BiocParallel]{BiocParallelParam} instance (from
    the \pkg{BiocParallel} package).
    This controls how the individual chunks are going to be processed
    once loaded in memory. Note that chunks are always loaded in memory
    sequentially. After being loaded in memory, the JSON lines in the
    current chunk are either processed sequentially (if \code{BPPARAM}
    is \code{NULL}) or in parallel (if \code{BPPARAM} is a
    \link[BiocParallel]{BiocParallelParam} instance).
  }
}

\details{
  RefSNP JSON files are made available by dbSNP for each release (a.k.a.
  build). For example, the RefSNP JSON files for dbSNP build 155 are
  available at \url{https://ftp.ncbi.nih.gov/snp/archive/b155/JSON/}.
  These files are compressed and have one RefSNP id per line.

  \code{extract_snvs_from_RefSNP_json()} will only consider RefSNP ids
  of variant type \code{"snv"}. Furthermore, for each RefSNP id of variant
  type \code{"snv"}, it will only consider its placements on sequences of
  type \code{"refseq_chromosome"}. The function will extract these placements
  plus their alleles and write them to output files in \code{dump_dir}.
  One output file will get created per sequence id.

  All the output files are tab-delimited files with one row per snv and
  the following columns:
  \enumerate{
    \item RefSNP id
    \item is preferred top level placement (PTLP)
    \item alleles position (zero-based)
    \item deleted sequence
    \item inserted sequences
  }
}

\value{
  The number of RefSNP ids processed (as an invisible integer).
  This is equal to the number of lines in the RefSNP JSON file.
}

\note{
  \code{extract_snvs_from_RefSNP_json()} is **very** slow! Depending on
  the chromosome and your machine, it will only process between 40 and 80
  RefSNP ids per second. At this speed it would take about 12 days just to
  process refsnp-chr1.json.bz2 (83578784 RefSNP ids). Using 9 workers (e.g.
  by setting \code{BPPARAM} to \code{MulticoreParam(9)}) makes this only
  about 3 times faster (i.e.  4 days instead of 12), which is also
  disappointing.
}

\examples{
json_file <- system.file("extdata", "refsnp-chrMT.json",
                         package="SNPlocsForge")
dump_dir <- file.path(tempdir(), "chrMT")
dir.create(dump_dir)

## Should take about 15 sec.:
extract_snvs_from_RefSNP_json(json_file, dump_dir, BPPARAM=MulticoreParam(6))

## Let's take a look at the output:
old_wd <- setwd(dump_dir)
dir()
cat(head(readLines("NC_012920.1.tab")), sep="\n")
read.delim("NC_012920.1.tab", header=FALSE, nrows=15)
setwd(old_wd)
}
\keyword{utilities}
