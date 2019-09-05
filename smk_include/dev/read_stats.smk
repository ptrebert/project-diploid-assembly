
WILDCARDS = glob_wildcards("/MMCI/TM/scratch/marschal/HG00514-assembly/pacbio/fastq/{filename}.bam")

rule master:
    input:
        expand('{filename}.pck',
                filename=WILDCARDS.filename)

rule compute_read_stats:
    input:
        '/MMCI/TM/scratch/marschal/HG00514-assembly/pacbio/fastq/{readset}.bam'
    output:
        '{readset}.pck'
    log: '{readset}.log'
    benchmark: '{readset}.rsrc'
    threads: 4
    run:
        exec = '/home/pebert/work/code/github/project-diploid-assembly/scripts/collect_read_stats.py'
        exec += ' --debug --chunk-size 25000 --num-cpu 4'
        exec += ' --output {output} --input-files {input}'
        exec += ' &> {log}'
        shell(exec)
