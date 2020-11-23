#!/usr/bin/env Rscript

suppressMessages(library(glue))
suppressMessages(library(GenomicRanges))
suppressMessages(library(LOLA))

path.regiondb.laptop = '/home/local/work/data/hgsvc/lola/regiondb/hg38'
path.userset = '/home/local/work/data/hgsvc/lola'
path.output = '/home/local/work/data/hgsvc/lola'

regiondb = loadRegionDB(path.regiondb.laptop)

# the following for loop produces the enrichments that are part of Supp. Table
# "assembly breaks / annotation enrichments"
# other variations below were exploring potential differences between HiFi/CLR
# for various ways of defining breaks (extending, merging etc.)
for (min.mapq in c(0, 10, 20)) {

    user.set = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_all.10kb.3c.bed", sep="/"))
    user.output = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_all.10kb.3c.lola.tsv", sep="/"))
    universe.set = glue("breaks_inv-{min.mapq}-60_all_all.10kb.3c.bed")

    col_idx = which(
                grepl(
                    universe.set,
                    regiondb$regionAnno$filename,
                    fixed=TRUE
                ),
                arr.ind=TRUE
            )
    

    universe = unlist(regiondb$regionGRL[col_idx])
    universe = disjoin(universe)

    user_roi = readBed(user.set)

    results = runLOLA(user_roi, universe, regiondb, cores=4, redefineUserSets=TRUE)

    writeCombinedEnrichment(results, outFolder=path.output, includeSplits=FALSE)

    generic.output = paste(path.output, "allEnrichments.tsv", sep="/")
    file.rename(generic.output, user.output)

}

# for (min.mapq in c(0)) {

#     for (tech in c('HiFi', 'CLR')) {

#         if (tech == 'HiFi') {
#             other.tech = 'CLR'
#         } else {
#             other.tech = 'HiFi'
#         }

#         user.set = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_{tech}.nogap.500bp.bed", sep="/"))
#         user.output = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_{tech}.nogap.500bp.lola.tsv", sep="/"))
#         universe.set = glue("breaks_inv-{min.mapq}-60_any_{other.tech}.nogap.500bp.bed")

#         col_idx = which(
#                     grepl(
#                         universe.set,
#                         regiondb$regionAnno$filename,
#                         fixed=TRUE
#                     ),
#                     arr.ind=TRUE
#                 )

#         universe = unlist(regiondb$regionGRL[col_idx])
#         universe = disjoin(universe)

#         user_roi = readBed(user.set)

#         results = runLOLA(user_roi, universe, regiondb, cores=4, redefineUserSets=TRUE)

#         writeCombinedEnrichment(results, outFolder=path.output, includeSplits=FALSE)

#         generic.output = paste(path.output, "allEnrichments.tsv", sep="/")
#         file.rename(generic.output, user.output)

#     }
# }

# #breaks_inv-10-60_any_CLR.ext500.mrg.bed
# for (min.mapq in c(10)) {

#     for (tech in c('HiFi', 'CLR')) {

#         if (tech == 'HiFi') {
#             other.tech = 'CLR'
#         } else {
#             other.tech = 'HiFi'
#         }

#         user.set = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_{tech}.ext500.mrg.bed", sep="/"))
#         user.output = glue(paste(path.userset, "breaks_inv-{min.mapq}-60_any_{tech}.ext500.mrg.lola.tsv", sep="/"))
#         universe.set = glue("breaks_inv-{min.mapq}-60_any_{other.tech}.ext500.mrg.bed")

#         col_idx = which(
#                     grepl(
#                         universe.set,
#                         regiondb$regionAnno$filename,
#                         fixed=TRUE
#                     ),
#                     arr.ind=TRUE
#                 )

#         universe = unlist(regiondb$regionGRL[col_idx])
#         universe = disjoin(universe)

#         user_roi = readBed(user.set)

#         results = runLOLA(user_roi, universe, regiondb, cores=4, redefineUserSets=TRUE)

#         writeCombinedEnrichment(results, outFolder=path.output, includeSplits=FALSE)

#         generic.output = paste(path.output, "allEnrichments.tsv", sep="/")
#         file.rename(generic.output, user.output)

#     }
# }
