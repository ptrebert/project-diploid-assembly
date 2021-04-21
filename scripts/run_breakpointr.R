#!/usr/bin/env Rscript

host.name = Sys.info()['nodename']

if (any(grepl('login', host.name, fixed=TRUE), grepl('submit', host.name, fixed=TRUE))) {
    stop('Running on cluster login/submit node - exiting')
}

suppressMessages(library(breakpointR))

args = commandArgs(trailingOnly=TRUE)

commit.hash.expect = args[6]
if (!is.na(commit.hash.expect)) {
    pkg.desc = library(help=breakpointR)
    pkg.desc.version = grep("RemoteRef", pkg.desc$info[[1]], ignore.case=TRUE, value=TRUE)
    pkg.desc.version = strsplit(pkg.desc.version, split=" +", fixed=FALSE)[[1]][2]
    install.match = grepl(commit.hash.expect, pkg.desc.version)
    if (!install.match) {
        stop(paste('Library breakpointR version mismatch: installed >>>', pkg.desc.version, '<<< versus expected', commit.hash.expect, sep=' '))
    }
}

bam.folder = args[1]
config.file = args[2]
output.folder = args[3]
num.cpu = args[4]
output.file = args[5]

breakpointr(
    inputfolder=bam.folder,
    outputfolder=output.folder,
    configfile=config.file,
    numCPU=num.cpu
)

exportRegions(
    datapath=file.path(output.folder, "data"),
    file=output.file,
    collapseInversions=TRUE,
    collapseRegionSize=5000000,
    minRegionSize=5000000,
    state="wc"
)

warnings()

quit(save='no')