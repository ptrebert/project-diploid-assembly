#!/usr/bin/env Rscript

if (is.element('SaaRclust', installed.packages()[,1])) {
    print('Removing previously installed version of SaaRclust')
    remove.packages('SaaRclust')
}

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

devtools::install_git(
    "git://github.com/daewoooo/SaaRclust.git",
    ref = git.commit,
    dependencies=FALSE,
    upgrade=FALSE
)

quit(save="no")
