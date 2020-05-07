#!/usr/bin/env Rscript

suppressMessages(library(SaaRclust))

args = commandArgs(trailingOnly=TRUE)

bed.file = args[1]
ref.genome = args[2]
output.folder = args[3]
contig.ordering = args[4]
plot.title = args[5]

if (is.na(contig.ordering)) {
    contig.ordering = FALSE
}

stopifnot(ref.genome == 'hg38')

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))


plot.clustering <- plotClusteredContigs(
    bedfile = bed.file,
    min.mapq = 10,
    bsgenome = BSgenome.Hsapiens.UCSC.hg38,
    report = 'clustering',
    title = paste('Clustering', plot.title, sep=': ')
)

if (as.logical(contig.ordering)) {

    plot.ordering <- plotClusteredContigs(
        bedfile = bed.file,
        min.mapq = 10,
        bsgenome = BSgenome.Hsapiens.UCSC.hg38,
        report = 'ordering',
        info.delim = '_',
        info.fields = c('assm', 'ctg', 'order', 'ID'),
        title = paste('Ordering', plot.title, sep=': ')

    )
}

plot.orienting <- plotClusteredContigs(
    bedfile = bed.file,
    min.mapq = 10,
    bsgenome = BSgenome.Hsapiens.UCSC.hg38,
    report = 'orienting',
    title = paste('Orientation', plot.title, sep=': ')
)

ggsave(filename = paste(output.folder, 'clustering.pdf', sep='.'), plot = plot.clustering, width = 16, height = 8)

if (as.logical(contig.ordering)) {
    ggsave(filename = paste(output.folder, 'ordering.pdf', sep='.'), plot = plot.ordering, width = 16, height = 8)
}

ggsave(filename = paste(output.folder, 'orienting.pdf', sep='.'), plot = plot.orienting, width = 16, height = 8)

warnings()

quit(save='no')