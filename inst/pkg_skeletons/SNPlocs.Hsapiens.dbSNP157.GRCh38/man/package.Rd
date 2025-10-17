\name{SNPlocs.Hsapiens.dbSNP157.GRCh38}
\docType{package}

\alias{SNPlocs.Hsapiens.dbSNP157.GRCh38-package}
\alias{SNPlocs.Hsapiens.dbSNP157.GRCh38}
\alias{.loadLoc}
\alias{.loadAlleles}
\alias{COMPATIBLE_BSGENOMES}


\title{The SNPlocs.Hsapiens.dbSNP157.GRCh38 package}

\description{
  Human SNP locations and alleles extracted from dbSNP Build 157 and
  placed on the GRCh38/hg38 assembly.

  IMPORTANT NOTE: You need a machine with at least 16G of RAM in
  order to use this package!
}

\details{
  The @TOTAL_SNPS@ SNPs in this package were extracted from the
  RefSNP JSON files for chromosomes 1-22, X, Y, and MT, located at
  \url{https://ftp.ncbi.nih.gov/snp/archive/b157/JSON/} (these files
  were created by NCBI in Nov 2022).

  These SNPs are compatible with packages \pkg{BSgenome.Hsapiens.NCBI.GRCh38}
  and \pkg{BSgenome.Hsapiens.UCSC.hg38}, that is, they can be "injected" in
  the \link[BSgenome]{BSgenome} objects contained in these packages.

  SNP positions and alleles are reported with respect to the \emph{plus}
  strand.

  Only SNPs of type \emph{snv} (\emph{single-nucleotide} variant a.k.a.
  \emph{single-base substitution}) were kept. Other variant types supported
  by dbSNP are: \emph{delins} (indel), \emph{ins} (insertion), \emph{del}
  (deletion), and \emph{mnv} (multiple nucleotide variation). These other
  variants are NOT included in \pkg{SNPlocs.Hsapiens.dbSNP157.GRCh38} but
  are available in the \pkg{XtraSNPlocs.Hsapiens.dbSNP157.GRCh38} package.
}

\note{
  The SNPs in this package can be "injected" in
  \code{BSgenome.Hsapiens.NCBI.GRCh38} or \code{BSgenome.Hsapiens.UCSC.hg38},
  and will land at the correct positions.

  See \code{?\link[BSgenome]{injectSNPs}} in the \pkg{BSgenome} software
  package for more information about the SNP injection mechanism.
}

\references{
  SNP Home at NCBI:
  \url{https://www.ncbi.nlm.nih.gov/snp}

  dbSNP Human Build 157 Release announcement (March 12, 2025):
  \url{https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2025q1/000239.html}

  The GRCh38.p14 assembly:
  \url{https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/}

  The hg38 genome at UCSC (based on GRCh38.p14, as of April 2022, but the
  UCSC folks could change this in the future and base hg38 on a more recent
  patch release of GRCh38):
  \url{http://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38}
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \pkg{XtraSNPlocs.Hsapiens.dbSNP157.GRCh38} package for SNPs of
          type other than \emph{snv}. Note that this package might not be
          available yet.

    \item \code{\link[BSgenome]{snpcount}} in the \pkg{BSgenome} software
          package for how to access the data stored in this package.

    \item \code{\link[Biostrings]{IUPAC_CODE_MAP}} in the \pkg{Biostrings}
          package.

    \item The \link[GenomicRanges]{GPos} class in the \pkg{GenomicRanges}
          package.

    \item \code{\link[BSgenome]{injectSNPs}} in the \pkg{BSgenome} software
          package for SNP injection.

    \item The \pkg{VariantAnnotation} software package to annotate variants
          with respect to location and amino acid coding.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BASIC USAGE
## ---------------------------------------------------------------------

snps <- SNPlocs.Hsapiens.dbSNP157.GRCh38
snpcount(snps)
seqinfo(snps)

## Get the positions and alleles of all SNPs on chromosomes 22 and MT:
snpsBySeqname(snps, seqnames=c("22", "MT"))

## Get the positions and alleles of all SNPs within some regions:
snpsByOverlaps(snps, GRanges(c("Y:230001-232000", "19:88501-89000")))

## ---------------------------------------------------------------------
## B. EXTRACT SNP INFORMATION FOR A SET OF RS IDS
## ---------------------------------------------------------------------

my_rsids <- c("rs2639606",  "rs75264089", "rs73396229", "rs55871206",
              "rs10932221", "rs56219727", "rs73709730", "rs55838886",
              "rs3734153",  "rs79381275", "rs1516535",  "rs74342513")

## Note that the first call to snpsById() takes a long time but
## subsequent calls are faster.
my_snps <- snpsById(snps, my_rsids)
my_snps

## Translate the IUPAC ambiguity codes used to represent the alleles
## into nucleotides:
IUPAC_CODE_MAP[mcols(my_snps)$alleles_as_ambig]

## ---------------------------------------------------------------------
## C. INJECTION IN THE REFERENCE GENOME
## ---------------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg38)
ref_genome <- BSgenome.Hsapiens.UCSC.hg38
ref_genome

alt_genome <- injectSNPs(ref_genome, "SNPlocs.Hsapiens.dbSNP157.GRCh38")
alt_genome  # note the additional line "with SNPs injected from..."

alphabetFrequency(ref_genome$chr22)
alphabetFrequency(alt_genome$chr22)

## Get the number of nucleotides that were modified by this injection
## (i.e. number of SNP loci on chr22):
neditAt(ref_genome$chr22, alt_genome$chr22)  # @NB_OF_SNP_LOCI_ON_CHR22@
}

\keyword{package}
