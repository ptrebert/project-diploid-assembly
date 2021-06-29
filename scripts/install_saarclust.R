#!/usr/bin/env Rscript

if (is.element('SaaRclust', installed.packages()[,1])) {
    print('Removing previously installed version of SaaRclust')
    remove.packages('SaaRclust')
}

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]
git.repo = args[2]

if (is.na(git.repo)) {
    git.repo = "git://github.com/daewoooo/SaaRclust.git"
}

devtools::install_git(
    git.repo,
    ref = git.commit,
    dependencies=FALSE,
    upgrade=FALSE
)

quit(save="no")
