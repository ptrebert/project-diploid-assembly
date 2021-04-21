#!/usr/bin/env Rscript

host.name = Sys.info()['nodename']

if (any(grepl('login', host.name, fixed=TRUE), grepl('submit', host.name, fixed=TRUE))) {
    stop('Running on cluster login/submit node - exiting')
}

suppressMessages(library(StrandPhaseR))

args = commandArgs(trailingOnly=TRUE)

commit.hash.expect = args[7]
if (!is.na(commit.hash.expect)) {
    pkg.desc = library(help=StrandPhaseR)
    pkg.desc.version = grep("RemoteRef", pkg.desc$info[[1]], ignore.case=TRUE, value=TRUE)
    pkg.desc.version = strsplit(pkg.desc.version, split=" +", fixed=FALSE)[[1]][2]
    install.match = grepl(commit.hash.expect, pkg.desc.version)
    if (!install.match) {
        stop(paste('Library StrandPhaseR version mismatch: installed >>>', pkg.desc.version, '<<< versus expected', commit.hash.expect, sep=' '))
    }
}

bam.folder = args[1]
config.file = args[2]
variant.calls = args[3]
wc.regions = args[4]
output.folder = args[5]
sample.individual = args[6]

strandPhaseR(
    inputfolder=bam.folder,
    configfile=config.file,
    outputfolder=output.folder,
    positions=variant.calls,
    fillMissAllele=variant.calls,
    WCregions=wc.regions,
    exportVCF=sample.individual
)

warnings()

quit(save='no')