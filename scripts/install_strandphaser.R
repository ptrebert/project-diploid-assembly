#!/usr/bin/env Rscript

if (is.element('StrandPhaseR', installed.packages()[,1])) {
    remove.packages('StrandPhaseR')
}

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

devtools::install_git("git://github.com/daewoooo/StrandPhaseR.git", ref = git.commit)

quit(save="no")
