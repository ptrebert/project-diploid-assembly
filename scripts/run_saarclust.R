#!/usr/bin/env Rscript

suppressMessages(library(SaaRclust))

args = commandArgs(trailingOnly=TRUE)

commit.hash.expect = args[4]
if (!is.na(commit.hash.expect)) {
    pkg.desc = library(help=SaaRclust)
    pkg.desc.version = grep("RemoteRef", pkg.desc$info[[1]], ignore.case=TRUE, value=TRUE)
    pkg.desc.version = strsplit(pkg.desc.version, split=" +", fixed=FALSE)[[1]][2]
    install.match = grepl(commit.hash.expect, pkg.desc.version)
    if (!install.match) {
        stop(paste('Library SaaRclust version mismatch: installed >>>', pkg.desc.version, '<<< versus expected', commit.hash.expect, sep=' '))
    }
}

scaffoldDenovoAssembly(
    configfile = args[1],
    bamfolder = args[2],
    outputfolder = args[3]
)

warnings()

quit(save='no')