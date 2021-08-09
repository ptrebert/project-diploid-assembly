#!/usr/bin/env Rscript

host.name = Sys.info()['nodename']

if (any(grepl('login', host.name, fixed=TRUE), grepl('submit', host.name, fixed=TRUE))) {
    stop('Running on cluster login/submit node - exiting')
}

suppressMessages(library(rjson))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Rsamtools))
suppressMessages(library(doParallel))
suppressMessages(library(matrixStats))
suppressMessages(library(assertthat))
suppressMessages(library(biovizBase))

args = commandArgs(trailingOnly=TRUE)

config = fromJSON(file = args[1])

r.source.folder = config$r_source_folder

hc.source.files = c(
    "calcProbs.R",
    "countDirectionalReads.R",
    "dumped_functions.R",
    "EMclust.R",
    "export.R",
    "findClusterPartners.R",
    "hardClust.R",
    "helperFuctions.R",
    "import.R",
    "importReads.R",
    "SaaRclust_evaluation_plots.R",
    "SaaRclust.R",
    "timedMessage.R",
    "utils.R",
    "wrapper_parallel.R",
    "wrapper.R"
)

hc.source.files = sapply(hc.source.files, function(x) {paste(r.source.folder, x, sep='/')})
hc.source.missing = sapply(hc.source.files, function(x) {if (!file.exists(x)) stop(paste('haploclust R source file does not exist:', x, sep=' '))})

hc.source.files = sapply(hc.source.files, function(x) suppressMessages(source(x)))

clust <- runSaaRclust(
    inputfolder=config$input_folder,
    outputfolder=config$output_folder,
    input_type=config$input_type,
    input.alignment.files=config$input_alignment_files,
    num.clusters=as.numeric(config$num_clusters),
    EM.iter=as.numeric(config$em_iter),
    numAlignments=as.numeric(config$num_alignments),
    hardclust.file=config$hard_clusters,
    softclust.file=config$soft_clusters,
    MLclust.file=config$MLclust,
    ss.clust.file=config$ss_clusters,
    clust.pairs.file=config$cluster_partners,
    wc.cells.file=config$wc_cells_clusters,
    numCPU=as.numeric(config$num_cpu)
)

warnings()

quit(save='no')