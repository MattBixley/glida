# create_ld_datasets.R
#
# Functions to create various population-specific LD patterns from VCF files.
# Requires PLINK2 to be installed and in the PATH.
#
# Nick Burns
# April, 2016

#' Recursively extracts sample IDs from the 1000 Genomes panel file.
#' Sample IDs written to 'sample_list.txt'.
#'
#' @param populations A vector of target populations.
#' @param n Integer (DEFAULT = 1). Decrementor for recursion.
#' @param allPops Boolean (DEFAULT = FALSE). If TRUE, no filtering of samples is performed.
#' @param first Boolean (DEFAULT = TRUE). Controls whether to overwrite or append to existing samples file.
#' @param popfile Location of the 1000 Genomes panel file.
#' @param outputFile Location of the output file, 'sample_list.txt'.
#' @return recursively iterates until all populations extracted.
#' @export
ldPopulation <- function (populations,
                          n = 1,
                          allPops = FALSE,
                          first = TRUE,
                          popfile = system.file("extdata", "1KGPopulations.panel", package="glida"),
                          outputFile = system.file("extdata", "sample_list.txt", package="glida")) {

    # recursive base-case
    if (n == 0) {
        return ()
    }

    if (allPops == TRUE) {
        command <- sprintf("cut -f1 %s > %s", popfile, outputFile)
        decr <- 0

    } else {
        # extract by population
        pipeout <- ifelse (first == TRUE, ">", ">>")
        command <- sprintf("%s %s %s %s %s",
                           system.file("bash", "getSampleIDs.sh", package = "glida"),
                           popfile,
                           populations[n],
                           pipeout,
                           outputFile)
        decr <- n - 1
    }

    system(command)
    return (ldPopulation(populations,
                         n = decr,
                         first = FALSE,
                         popfile = popfile,
                         outputFile = outputFile))
}

#' Downloads a defined region from 1000 Genomes.
#'
#' @param chromosome integer.
#' @param start CHARACTER. Please note, not integer.
#' @param end CHARACTER. Please note, not integer.
#' @param outputVCF. Filename.
#' @export
ldDownload <- function (chromosome = NULL,
                         start = NULL,
                         end = NULL,
                         outputVCF = NULL) {

    lclVCF <- if (missing(outputVCF)) sprintf("./Genotype_%s_%s-%s", chromosome, start, end)
              else outputVCF

    print(lclVCF)
    ldCmd <- sprintf("%s %s %s %s %s",
                     system.file("bash", "get1KGRegion.sh", package="glida"),
                     chromosome,
                     start,
                     end,
                     lclVCF)
    print(ldCmd)
    system(ldCmd)
}


#' Downloads the 1000 Genomes data for the region of interest
#' and calculates LD between all pairs of SNPs using PLINK2.
#' Writes results out to ./Data/Genotype_<chromosome>_<start>_<end>.ld
#' Requires:
#'     - bash script: ./exec/getLD.sh
#'     - sample_list.txt must be present
#'
#' @param chromosome Integer Chromosome number
#' @param start Integer The beginning of the region of interest (in base pairs).
#' @param end Integer The end of the region of interest (in base pairs).
#' @export
ldByRegion <- function (chromosome = NULL,
                        start = NULL,
                        end = NULL,
                        download = FALSE,
                        vcfFile = NULL,
                        sampleFile = NULL) {

    lclSamples <- if (missing(sampleFile)) system.file("extdata", "sample_list.txt", package="glida")
                  else samplesFile
    lclVCF <- if (missing(vcfFile)) sprintf("Genotype_%s_%s-%s.vcf", chromosome, start, end)
              else vcfFile

    if (download == TRUE) {
        ldDownload(chromosome = chromosome,
                    start = start,
                    end = end,
                    outputVCF = lclVCF)
    }

    ldCmd <- sprintf("%s %s %s",
                     system.file("bash", "getLD.sh", package = "glida"),
                     lclVCF,
                     lclSamples)
    system(ldCmd)
}

#' Given a lead SNP, calculates LD of all other SNP (relative to the lead SNP)
#' within a defined window (distance).
#' Requires:
#'   - bash script: ../exec/getProxy.sh
#'   - a VCF file (1000 Genomes data) for given region.
#'     Will download if not present
#'
#' @param leadSNP Character. The lead SNP of interest.
#' @param vcf filename (DEFAULT = ""). Full (or relative) path of the VCF file
#' @param downloadVCF Boolean (DEFAULT = FALSE). If TRUE, download the 1000 Genomes data for given region.
#' @param chromosome Integer (DEFAULT = NULL. If donwloadVCF == TRUE, this must be set. See ldByRegion().
#' @param start Integer (DEFAULT = NULL). If downloadVCF == TRUE, this must be set. See ldByRegion().
#' @param end Integer (DEFAULT = NULL). If downloadVCF == TRUE, this must be set. See ldByRegion().
#' @export
ldProxy <- function(leadSNP,
                    vcfFile = NULL,
                    ldOutput = NULL,
                    download = FALSE,
                    chromosome = NULL,
                    start = NULL,
                    end = NULL,
                    sampleFile = NULL) {

    lclVCF <- if (missing(vcfFile))
                    sprintf("Genotype_%s_%s-%s.vcf", chromosome, start, end)
              else vcfFile
    lclSamples <- if (missing(sampleFile))
                    system.file("extdata", "sample_list.txt", package = "glida")
                  else sampleFile
    lclOutput <- if (missing(ldOutput))
                    unlist(strsplit(lclVCF, "\\."))[1]
                 else ldOutput

    # If download is TRUE, or vcf file does not exists, download.
    if (download == TRUE | !file.exists(lclVCF)) {
        ldDownload(chromosome, start, end, outputVCF = lclVCF)
    }

    # calculate LD using PLINK2 (see ../exec/getProxy.sh)
    ldCmd <- sprintf("%s %s %s %s %s",
                     system.file("bash", "getProxy.sh", package="glida"),
                     leadSNP,
                     lclVCF,
                     lclSamples,
                     lclOutput)
    system(ldCmd)
}
