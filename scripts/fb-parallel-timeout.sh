#!/bin/bash

if [[ $# -lt 3 ]];
then
    echo "usage: $0 [regions file] [ncpus] [timeout] [logfile] [freebayes arguments]"
    echo
    echo "Run freebayes in parallel over regions listed in regions file, using ncpus processors."
    echo "Will merge and sort output, producing a uniform VCF stream on stdout.  Flags to freebayes"
    echo "which would write to e.g. a particular file will obviously cause problms, so caution is"
    echo "encouraged when using this script."
    echo
    echo "This script: adapted by Tobias Marschall"
    echo
    echo "For original version, see this github repo:"
    echo
    echo "https://github.com/ekg/freebayes/blob/master/scripts/freebayes-parallel"
    echo
    exit
fi

regionsfile=$1
shift
ncpus=$1
shift
timeout=$1
shift
logfile=$1
shift

command=("freebayes" "$@")

( cat "$regionsfile" | parallel -k --joblog "$logfile" -j "$ncpus" "timeout ${timeout} ${command[@]}" --region {} ) | vcffirstheader | vcfstreamsort -w 10000 | vcfuniq