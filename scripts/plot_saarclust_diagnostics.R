#!/usr/bin/env Rscript

suppressMessages(library(SaaRclust))

args <- commandArgs(trailingOnly=TRUE)

bed.file <- args[1]
ref.genome <- args[2]
output.folder <- args[3]
plot.title <- args[4]
sample.sex <- args[5]
haploid.assembly <- args[6]

if (is.na(haploid.assembly)) {
    haploid.assembly <- FALSE
} else {
    haploid.assembly <- as.logical(haploid.assembly)
}

if (is.na(sample.sex)) {
    sample.sex <- 'unknown'
}

if (sample.sex == 'unknown') {
    plot.chromosomes = paste0('chr', c(1:22))
} else if (sample.sex == 'male') {
    plot.chromosomes = paste0('chr', c(1:22, 'X', 'Y'))
} else if (sample.sex == 'female') {
    plot.chromosomes = paste0('chr', c(1:22, 'X'))
} else {
    stop(paste0('Unknown sample sex ', sample.sex))
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
        info.delim = '_',
        info.fields = c('cluster.ID', 'contig.ID', 'sample.sex', 'scl.info'),
        col.by = 'cluster.ID',
        title = paste('Clustering', plot.title, sep=': '),
        chromosomes = plot.chromosomes
    )

    plot.orienting <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'orienting',
        info.delim = '_',
        info.fields = c('cluster.ID', 'contig.ID', 'sample.sex', 'scl.info'),
        col.by = 'cluster.ID',
        title = paste('Orientation', plot.title, sep=': '),
        chromosomes = plot.chromosomes
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
        title = paste('Clustering', plot.title, sep=': '),
        chromosomes = plot.chromosomes
    )

    plot.ordering <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'ordering',
        info.delim = '_',
        info.fields = c('cluster.SRC', 'contig.ID', 'order', 'cluster.ID'),
        title = paste('Ordering', plot.title, sep=': '),
        chromosomes = plot.chromosomes
    )

    plot.orienting <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'orienting',
        info.delim = '_',
        info.fields = c('cluster.SRC', 'contig.ID', 'order', 'cluster.ID'),
        title = paste('Orientation', plot.title, sep=': '),
        chromosomes = plot.chromosomes
    )
}

if (!is.null(plot.clustering)) {
    ggsave(filename = paste(output.folder, 'clustering.pdf', sep='.'), plot = plot.clustering$plot, width = 16, height = 8)
}

if (!is.null(plot.ordering)) {
    ggsave(filename = paste(output.folder, 'ordering.pdf', sep='.'), plot = plot.ordering$plot, width = 16, height = 8)
}

if (!is.null(plot.orienting)) {
    ggsave(filename = paste(output.folder, 'orienting.pdf', sep='.'), plot = plot.orienting$plot, width = 16, height = 8)
}

warnings()

quit(save='no')