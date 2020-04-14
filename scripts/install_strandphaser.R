#!/usr/bin/env Rscript

if (is.element('StrandPhaseR', installed.packages()[,1])) {
    print('Removing previously installed version of StrandPhaseR')
    remove.packages('StrandPhaseR')
}

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

devtools::install_git(
    "git://github.com/daewoooo/StrandPhaseR.git",
    ref = git.commit,
    dependencies=FALSE,
    upgrade=FALSE
)

quit(save="no")
