\name{extract_raw_snps_from_json}

\alias{extract_raw_snps_from_json}

\title{Extract "raw" SNPs from the JSON files provided by dbSNP}

\description{
  Extract "raw" SNPs from the RefSNP JSON files provided by dbSNP
  at \url{https://ftp.ncbi.nih.gov/snp/latest_release/JSON/}.
}

\usage{
extract_raw_snps_from_json(con, out="", n=50000,
                           assembly="GRCh38.p13", paranoid=FALSE)
}

\arguments{
  \item{con}{
    File path or connection to a RefSNP JSON file (compressed files
    are supported).
  }
  \item{out}{
    Output file (will be written in tab-delimited format).
  }
  \item{n}{
    How many JSON lines to load at once in memory. Set to -1 to load the
    entire file in memory (strongly discouraged on a big JSON file!)
  }
  \item{assembly}{
    The NCBI assembly that the alleles to extract from the file are
    expected to be placed on.
  }
  \item{paranoid}{
    Set to \code{TRUE} to perform additional sanity checks on the
    extracted alleles. Note that this is at the cost of a 50-60\% slowdown!
  }
}

\details{
  RefSNP JSON files have one RefSNP id per line.
  \code{extract_raw_snps_from_json()} will write one line to the output file
  (tab-delimited) for each line in the input file (RefSNP JSON).

  Columns written to the tab-delimited file:
  \enumerate{
    \item RefSNP id
    \item variant type
    \item sequence name
    \item alleles position
    \item deleted sequence
    \item inserted sequences
  }
}

\value{
  The number of RefSNP ids processed (as an invisible integer).
  This is equal to the number of lines in the input and output files (both
  files have one line per RefSNP id).

}

\note{
  \code{extract_raw_snps_from_json()} is very slow!

  Here are some timings:
  \itemize{
    \item It takes about 12 seconds on refsnp-chrMT.json (7333 RefSNP ids).
    \item It takes about 1 h 12 min on refsnp-chrY.json (2480654 RefSNP ids).
    \item It takes about 7 h 42 min on refsnp-chr22.json (14305736 RefSNP ids).
    \item It takes about 21 h on refsnp-chr10.json.bz2 (??? RefSNP ids).
    \item It takes about 35 h on refsnp-chr1.json.bz2 (??? RefSNP ids).
    \item So all in all, that means it takes about 16-17 days to process
          the RefSNP JSON files for all Human chromosomes (25 files,
          ??? RefSNP ids).
  }
}

\examples{
json_file <- system.file("extdata", "refsnp-chrMT.json",
                         package="SNPlocsForge")
out_file <- tempfile(fileext=".tab")

## Should take about 15 sec.:
extract_raw_snps_from_json(json_file, out_file, n=1000)

cat(head(readLines(out_file)), sep="\n")

read.delim(out_file, header=FALSE, nrows=15)
}
\keyword{utilities}
