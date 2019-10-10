#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

devtools::install_git("git://github.com/daewoooo/StrandPhaseR.git", ref = git.commit)

quit(save="no")
