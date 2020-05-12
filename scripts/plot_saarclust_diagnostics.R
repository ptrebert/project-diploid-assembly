#!/usr/bin/env Rscript

suppressMessages(library(SaaRclust))

args = commandArgs(trailingOnly=TRUE)

bed.file = args[1]
ref.genome = args[2]
output.folder = args[3]
plot.title = args[4]
haploid.assembly = args[5]

if (is.na(haploid.assembly)) {
    haploid.assembly = FALSE
} else {
    haploid.assembly = as.logical(haploid.assembly)
}

stopifnot(ref.genome == 'hg38')

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

plot.clustering <- NULL
plot.ordering <- NULL
plot.orienting <- NULL

if (!haploid.assembly) {

    plot.clustering <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'clustering',
        title = paste('Clustering', plot.title, sep=': ')
    )

    plot.orienting <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'orienting',
        title = paste('Orientation', plot.title, sep=': ')
    )
} else {

    plot.clustering <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'clustering',
        info.delim = '_',
        info.fields = c('cluster.SRC', 'contig.ID', 'order', 'cluster.ID'),
        col.by = 'cluster.ID',
        title = paste('Clustering', plot.title, sep=': ')
    )

    plot.ordering <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'ordering',
        info.delim = '_',
        info.fields = c('cluster.SRC', 'contig.ID', 'order', 'cluster.ID'),
        title = paste('Ordering', plot.title, sep=': ')
    )

    plot.orienting <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'orienting',
        info.delim = '_',
        info.fields = c('cluster.SRC', 'contig.ID', 'order', 'cluster.ID'),
        title = paste('Orientation', plot.title, sep=': ')
    )
}

if (!is.null(plot.clustering)) {
    ggsave(filename = paste(output.folder, 'clustering.pdf', sep='.'), plot = plot.clustering, width = 16, height = 8)
}

if (!is.null(plot.ordering)) {
    ggsave(filename = paste(output.folder, 'ordering.pdf', sep='.'), plot = plot.ordering, width = 16, height = 8)
}

if (!is.null(plot.orienting)) {
    ggsave(filename = paste(output.folder, 'orienting.pdf', sep='.'), plot = plot.orienting, width = 16, height = 8)
}

warnings()

quit(save='no')