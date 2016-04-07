#' Produces a heatmap based on an LD dissimilarity matrix.
#'
#' @param ldMatrix matrix. Dissimilarity matrix of LD (R2) between SNPs.
#' @return heatmap + dendrogram
#' @export
ldHeatmap <- function (ldMatrix, clustMethod = "ward.D", distMethod = "euclidean") {

    colourPalette <- RColorBrewer::brewer.pal(5, "Greys")[5:1]
    getClust <- function (x) hclust(x, method = clustMethod)
    getDist <- function (x) dist(x, method = distMethod)

    ldMap <- heatmap(ldMatrix,
                     col = colourPalette,
                     Colv = 1:ncol(ldMatrix),
                     hclustfun = getClust,
                     distfun = getDist)

    return (ldMap)
}

#' Clusters SNPs based on an LD dissimilarity matrix. Returns a pvclust cluster object.
#' The pvclust object can be used to find:
#'     - similar (and dissimilar) SNPs (based on LD)
#'     - can be passed to pvclust::pvpick to find clusters with high/low pvalues
#'     - can be plotted with ldDendrogram()
#'
#' @param ldMatrix matrix. LD dissimilarity matrix.
#' @param clustMethod character. Method for clustering ("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid")
#' @param distMethod character. Method for calculating distances between SNPs. see pvclust help.
#' @return pvclust cluster object.
ldCluster <- function (ldMatrix, clustMethod = "ward.D", distMethod = "euclidean") {

    lclCluster <- pvclust::pvclust(ldMatrix,
                                   method.hclust = clustMethod,
                                   method.dist = distMethod,
                                   quiet=TRUE,
                                   nboot=100)

    return (lclCluster)
}

#' Performs heirarchical clustering of LD dissimilarity matrix and returns a dendrogram.
#'
#' @param ldClusters pvclust object. Return values from ldClsuter().
#' @param plotTitle character (DEFAULT = LD-based Clustering). Title to be displayed with the dendrogram.
#' @return dendrogram
#' @export
ldDendrogram <- function (ldClusters, plotTitle = "LD-based Clustering") {
    ldPlot <- plot(ldClusters,
                   print.pv = TRUE,
                   main = plotTitle,
                   cex = 0.75)

    return (ldPlot)
}

#' Produces a locus zoom-like plot, where the y-axis is LD (R2) rather than p-values.
#' LD data should be generated using ldProxy(). See relevant documentation.
#' Plots can be annotated with nearby genes and eQTLs (to be implemented in the future).
#'
#' @param ldData data frame. Standard PLINK LD output. See ldRead().
#' @param ldThreshold float (DEFAULT = 0.9). SNPs names are displayed for SNPs in LD above this threshold.
#' @return zoom: a ggplot2 plot object.
#' @export
ldZoom <- function (ldData, ldThreshold=0.9) {

    ldData$POS <- ldData$BP_B / 1000000

    # base plot, (POS, R2) with some custom colouring based on R2.
    zoom <- ggplot2::ggplot(ldData, ggplot2::aes(x = POS, y = R2)) +
        ggplot2::geom_point(ggplot2::aes(size = R2, colour = -R2),
                            alpha = 0.8) +
        ggplot2::scale_colour_gradientn(colours=rainbow(3)) +
        ggplot2::xlab("Position (Mb)") +
        ggplot2::guides(size=FALSE, colour=FALSE) +
        ggplot2::theme_bw()


    # include SNP names, where R2 > ldThreshold
    ldWindow <- ldData[ldData$R2 > ldThreshold, ]
    zoom <- zoom +
        ggrepel::geom_text_repel(data = ldWindow,
                                 ggplot2::aes(x = POS, y = R2, label = SNP_B),
                                 colour="grey10",
                                 size=3)

    # annotate plot with 2 vertical lines that mark the window of association
    # the window of association is the region where all SNPs R2 > ldThreshold
    zoom <- zoom +
        ggplot2::geom_vline(xintercept = min(ldWindow$POS),
                            color="darkgray",
                            linetype = "dashed") +
        ggplot2::geom_vline(xintercept = max(ldWindow$POS),
                            color="darkgray",
                            linetype = "dashed")

    return (zoom)
}

#' Adds gene annotations to an existing locus zoom, or ldZoom plot.
#'
#' @param zoom GGPLOT2 plot object. See ldZoom() as example.
#' @param genes data frame. Genes within the target region. See fromUCSCEnsemblGenes().
#' @param geneType character (DEFAULT = 'All'). Filter the displayed genes by gene type (not implemented yet).
#' @return zoom. A GGPLOT2 locus-zoom like plot, with gene annotations.
#' @export
geneAnnotation <- function (zoom, genes, geneType = "All") {

    # wrangle the genes into a plotable format
    lclGenes <- if (geneType == "All") genes
                else genes[genes$geneType == genetype, ]
    lclGenes$GeneStart <- lclGenes$GeneStart / 1000000
    lclGenes$GeneEnd <- lclGenes$GeneEnd / 1000000

    # set positions on y-axis.
    # This is a dirty fudge in order to get a semi-readable layout
    # STRATEGY:
    #     1. create K rows of genes, maximum of 5 rows.
    #     2. genes are ordered by geneStart, so layout succesive genes
    #        row-by-row (thus nearby genes should be on separate rows)
    #     3. Set a bound on the row position to be between 0 and 1, and then
    #        further limit this to 0:-0.5.
    N <- nrow(lclGenes)
    K <- ceiling(N / 5)
    bounds <- N / K
    lclGenes$Yvalues <- - (0.5 / bounds) * rep(1:bounds, length.out = N)

    # finally, colour the gene labels by geneType
    colourByType <- factor(lclGenes$geneCategory)
    colourByType <- as.integer(relevel(colourByType,
                                       ref = "protein_coding"))

    zoom <- zoom +
        ggplot2::geom_segment(data = lclGenes,
                              ggplot2::aes(x = GeneStart, xend = GeneEnd,
                                           y = Yvalues, yend = Yvalues),
                              size = 1, colour = "grey30", alpha = 0.5) +
        ggrepel::geom_text_repel(data = lclGenes,
                                 ggplot2::aes(x = (GeneStart + GeneEnd) / 2,
                                              y = Yvalues,
                                              label = GeneName),
                                 colour = colourByType,
                                 size = 3)

    return (zoom)
}
