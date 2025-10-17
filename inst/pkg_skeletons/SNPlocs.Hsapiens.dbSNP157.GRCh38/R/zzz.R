###
###

.genome <- "GRCh38.p14"

.onLoad <- function(libname, pkgname)
{
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make Seqinfo object.
    seqinfo_path <- file.path(extdata_dirpath, "seqinfo.txt")
    seqinfo <- BSgenome:::read_seqinfo_table(seqinfo_path, .genome)

    ## Make GenomeDescription object. See
    ##   https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/
    ## for more information about GRCh38.p14.
    reference_genome <- GenomeDescription(
        organism="Homo sapiens",
        common_name="Human",
        provider="NCBI",
        provider_version=.genome,
        release_date="2022/02/03",
        release_name="GRCh38 Patch Release 14",
        seqinfo=seqinfo
    )

    ## Make list of "sequence name translation tables".
    seqlevels <- seqlevels(seqinfo)
    compatible_genomes <- list(
        BSgenome.Hsapiens.NCBI.GRCh38=
            structure(seqlevels,
                      names=c(1:22, "X", "Y", "MT")),
        BSgenome.Hsapiens.UCSC.hg38=
            structure(seqlevels,
                      names=paste0("chr", c(1:22, "X", "Y", "M")))
    )

    ## Make and export SNPlocs object.
    snps <- newSNPlocs(
        provider="dbSNP",
        provider_version="dbSNP Human Build 157",
        release_date="March 2025",
        release_name="dbSNP Human Build 157",
        source_data_url="https://ftp.ncbi.nih.gov/snp/archive/b157/JSON/",
        download_date="Oct 6-7, 2025",
        reference_genome=reference_genome,
        compatible_genomes=compatible_genomes,
        data_pkgname=pkgname,
        data_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, snps, envir=ns)
    namespaceExport(ns, objname)

    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## SNP injection interface (see BSgenome package for how this interface
    ## is used).
    ## TODO: The SNP injection mechanism needs to be modernized and the
    ## modernized version should not use this ugly interface anymore.
    ##

    ## Define and export .loadLoc() and .loadAlleles().
    .loadLoc <- function(seqname)
        BSgenome:::.load_raw_snplocs(snps, seqname, TRUE)$loc
    objname <- ".loadLoc"
    assign(objname, .loadLoc, envir=ns)
    namespaceExport(ns, objname)

    .loadAlleles <- function(seqname)
        BSgenome:::.load_raw_snplocs(snps, seqname, TRUE)$alleles
    objname <- ".loadAlleles"
    assign(objname, .loadAlleles, envir=ns)
    namespaceExport(ns, objname)

    ## Export 'compatible_genomes' as 'COMPATIBLE_BSGENOMES'.
    objname <- "COMPATIBLE_BSGENOMES"
    assign(objname, compatible_genomes, envir=ns)
    namespaceExport(ns, objname)
}

