#' Establishes a connection to UCSC's public MySQL server
openConn <- function (genomeBuild = "hg19") {

#     # turn off warnings that the DB connector throws
#     # (annoying non-informative messages)
#     preState <- options("warn")
#     options(warn = -1)

    # connect to the database
    #drv <- DBI::dbDriver("MySQL")
    conn <- RMySQL::dbConnect(RMySQL::MySQL(),
                              user="genome",
                              host="genome-mysql.cse.ucsc.edu",
                              dbname=genomeBuild,
                              password="")

#     # reset the warnings
#     options(warn = preState$warn)
    return (conn)
}

#' Close an existing connection to UCSC's public MySQL server
closeConn <- function (conn) {
    RMySQL::dbDisconnect(conn)
}

#' Given an existing database connection, retrieve query results from UCSC
#'
#' @param query character. ANSI-SQL compliant query string.
#' @return results data frame.
#' @export
queryUCSC <- function (query) {

    try({
        # turn off warnings that the DB connector throws
        # (annoying non-informative messages)
        preState <- options("warn")
        options(warn = -1)
        
        conn <- openConn()
        results <- DBI::dbGetQuery(conn, query)
        closeConn(conn)
        
        # reset the warnings
        options(warn = preState$warn)

        return (results)
    })

    return ("There was an error querying from UCSC")
}

#' Find all genes within a defined chromosomal region.
#' Queries UCSC's knownGene table.
#' THIS ISN'T THE ONE I WOULD RECOMMEND USING. Check out fromUCSCEnsemblGenes.
#'
#' @param chromosome int. Chromosome number [1..22]
#' @param start int. Start of the region, base position
#' @param end int. End of region, base position
#' @return query string. To be passed as argument to queryUCSC().
fromUCSCKnownGene <- function (chromosome = NULL, start = NULL, end = NULL) {

    if (any(missing(chromosome), missing(start), missing(end))) {
        return ("You must specify 'chromosome', 'start' and 'end'")
    }

    # Queries the knownGene table from the UCSC Genome Browser
    query <- paste0("  SELECT ", chromosome, " as CHR,
                              kg.name as GeneName,
                              kg.txStart as GeneStart,
                              kg.txEnd as GeneEnd
                      FROM knownGene kg
                      WHERE kg.chrom = 'chr", chromosome, "'
                        AND ( (kg.txStart BETWEEN ", start, " AND ", end, ") OR
                              (kg.txEnd BETWEEN ", start, " AND ", end, ") OR
                              (kg.txStart < ", start, " AND txEnd > ", start, "))
                      GROUP BY kg.txStart, kg.txEnd
                      ORDER BY kg.txStart ASC;
                   ")
    return (query)
}

#' Find all genes within a defined chromosomal region.
#' Queries UCSC's refGene table.
#'
#' @param chromosome int. Chromosome number [1..22]
#' @param start int. Start of the region, base position
#' @param end int. End of region, base position
#' @return query string. To be passed as argument to queryUCSC().
fromUCSCRefGene <- function (chromosome = NULL, start = NULL, end = NULL) {

    if (any(missing(chromosome), missing(start), missing(end))) {
        return ("You must specify 'chromosome', 'start' and 'end'")
    }

    query <- paste0("  SELECT ", chromosome, " as CHR, ", "
                              name2 as GeneName,
                              txStart as GeneStart,
                              txEnd as GeneEnd,
                              strand,
                              exonCount
                      FROM refGene as rg
                      WHERE chrom = 'chr", chromosome, "'
                        AND ( (txStart BETWEEN ", start, " AND ", end, ") OR
                              (txEnd BETWEEN ", start, " AND ", end, ") OR
                              (txStart < ", start, " AND txEnd > ", start, "))
                      GROUP BY name2
                      ORDER BY txStart ASC;
                   ")
    return (query)
}

#' Find all genes within a defined chromosomal region.
#' Queries UCSC's Ensembl tables and includes additional information such as,
#' geneType, geneStatus and transcriptClass.
#' This is the best of these queries I think.
#'
#' @param chromosome int. Chromosome number [1..22]
#' @param start int. Start of the region, base position
#' @param end int. End of region, base position
#' @return query string. To be passed as argument to queryUCSC().
#' @export
fromUCSCEnsemblGenes <- function (chromosome = NULL, start = NULL, end = NULL) {

    if (any(missing(chromosome), missing(start), missing(end))) {
        return ("You must specify 'chromosome', 'start' and 'end'")
    }

    query <- paste0("  SELECT ", chromosome, " as CHR,
                                 wg1.name2 as GeneName,
                                 wg1.txStart as GeneStart,
                                 wg1.txEnd as GeneEnd,
                                 wg1.exonCount,
                                 wg2.geneType,
                                 CASE
                                     WHEN (wg2.geneType = 'protein_coding') THEN 'protein_coding'
                                     WHEN (wg2.geneType LIKE '%RNA' OR wg2.geneType LIKE '%rna') THEN 'rna_gene'
                                     WHEN (wg2.geneType LIKE '%pseudogene') THEN 'psuedogene'
                                     ELSE 'other_gene'
                                 END as geneCategory,
                                 wg2.geneStatus
                      FROM wgEncodeGencodeBasicV19 as wg1
                          INNER JOIN wgEncodeGencodeAttrsV19 as wg2 ON wg1.name = wg2.transcriptId
                      WHERE wg1.chrom = 'chr", chromosome, "'
                        AND ( (wg1.txStart BETWEEN ", start, " AND ", end, ") OR
                              (wg1.txEnd BETWEEN ", start, " AND ", end, ") OR
                              (wg1.txStart < ", start, " AND wg1.txEnd > ", start, "))
                      GROUP BY wg1.name2
                      ORDER BY wg1.txStart ASC;
                   ")
    return (query)
}

#' Extracts genes within region of interest from biomaRt.
#' NOTE: that biomaRt is a suggested dependecy only. So the package will run
#' perfectly fine if you do not have biomaRt installed. However, this
#' function will not execute.
#'
#' @param chromosome int. Chromosome number [1..22]
#' @param start int. Start of the region, base position
#' @param end int. End of region, base position
#' @param genomeBuild string. (DEFAULT = 'grch37'). Specify the build to use. Note the build can make a dramatic difference to the returned gene set.
#' @return gene set.
#' @export
fromBiomaRt <- function (chromosome = NULL, start = NULL, end = NULL, genomeBuild = "grch37") {

    if (any(missing(chromosome), missing(start), missing(end))) {
        return ("You must specify 'chromosome', 'start' and 'end'")
    }

    lclHost <- if (genomeBuild == "grch38") "ensembl.org"
               else sprintf("%s.ensembl.org", genomeBuild)

    ensemblMart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                    host=lclHost,
                                    path="/biomart/martservice" ,
                                    dataset="hsapiens_gene_ensembl")
    genes <- biomaRt::getBM(attributes=c("chromosome_name",
                                          "start_position",
                                          "end_position",
                                          "external_gene_name",
                                          "transcript_count",
                                          "gene_biotype"),
                             filters = c("chromosome_name",
                                         "start",
                                         "end"),
                             values = list(chromosome=chromosome,
                                           start=start,
                                           end=end),
                             mart = ensemblMart)
    return (genes)
}
