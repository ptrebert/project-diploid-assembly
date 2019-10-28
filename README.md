# Project repository: reference-free diploid assembly of a human genome

## Citation

DOI of this repository: TBA

## Scope of this repository

This repository contains only the Snakemake pipeline code plus some auxiliary scripts to go from raw
input data to polished haploid assemblies. Any self-contained, general purpose software tool used in
this pipeline is either available via conda/bioconda, or via github (in which case, it will be
installed as part of the pipeline workflow).

## Setup and naming requirements

1. Checkout this repository (git clone, or download zip)
2. Create the conda environment stored in `environment/conda/conda_dipassm.yml`
3. Adapt Snakemake configuration files and profiles according to your computing environment:
    - various profiles are stored under `environment/snakemake`, either for individual servers, or
      for compute clusters with job schedulers SGE, SLURM, or PBS/TORQUE
    - adapt the values for the maximal available number of CPU cores available in your environment
      in the Snakemake config file (prefix `smk_config/run_env/smk_cfg_`)
    - all resource-intensive jobs should be annotated with resource requirements. You likely need
      this info to write a proper cluster config: `mem_total_mb` and `mem_per_cpu_mb`,
      depending on whether you need to specify total amount of memory or memory per CPU core.
      **Important note**: for some jobs, the total amount of memory is not scaled with
      the number of CPUs used, because memory consumption is rather a function of the input data
      (say, number of reads and average read length). If jobs fail due to insufficient resources,
      then the rules need to be explicitly added to the cluster config.
4. If you want to use your own data, make sure to create the appropriate symlinks as needed in the
   file system folder hierarchy, and make Snakemake aware of these files by creating a rule in
   `smk_include/link_data_sources.smk`
5. If you use your own data, you probably need to adapt the build targets for Snakemake listed in
   `smk_include/results_child.smk` and `smk_include/results_parents.smk`. Another option is to change
   the `master` rule in the main `Snakefile`.
6. Run the pipeline using your configuration and profile.

#### File naming requirements

Currently, the pipeline has been designed for file names that follow this pattern:

**individual**\_**project**\_**platform**\_**suffix** + file extensions, usually separated by "."

Examples:

  - HG00733_sra_pbsq1-clr_1000.fastq.gz
    - individual *HG00733*
    - project *sra* (data downloaded from NCBI/SRA)
    - platform *pbsq1-clr*: PacBio Sequel-1 platform, CLR reads
    - suffix *1000*: since this is a read data sets, this indicates
      that the data is not downsampled (read as: 100.0%).
      Downsampling to 36.7% of all reads would be denoted with suffix *0367*
  - HG00733_sra_pbsq1-clr_sqa-wtdbg.fasta
    - *(all as above)*
    - suffix *sqa-wtdbg*: a squashed assembly file created with wtdbg
      created from the complete readset (1000 omitted)
  - HG00733_1kg_il25k-150di_P0IIL081_ERR1295766.filt.sam.bam
    - individual *HG00733*
    - project *1kg*: 1000 Genomes
    - platform *il25k-150di*: Illumina HiSeq 2500, 150 bp reads
    - suffix *P0IIL081_ERR1295766*: library and run ID; the run ID is
      stripped off in a later step
    - file extension *filt.sam.bam*: the pipeline distinguishes between
      regular "SAM/BAM" files, and the Pacbio-native variant of (unaligned) BAM
      files containing additional quality information (extension "pbn.bam")

For the evaluation of the results against a known reference, the pipeline distinguishes
two variants of the human genome reference assembly v38:

1. GRCh38_GCA_p13: untouched GRCh38 reference from Genbank, patch level 13; this version
   is used as evaluation reference
2. hg38_GCA_p13: modified version of the above reduced to the chromosomes as specified
   in the Snakemake configuration file (by default: 1-22, X, Y). This version is used
   for alignment steps (that are not part of the evaluation!)

## Common points of failure

1. Snakemake
  - The pipeline heavily uses Snakemake's `checkpoint` feature for data-dependent execution of
    subsequent rules. Unfortunately, this more or less recent feature is not yet as mature*
    as the older parts of Snakemake's code. A start-to-finish run of the pipeline may
    fail in between because of a faulty checkpoint evaluation. In this case, waiting for all
    remaining jobs to finish and restarting the pipeline often leads to a successful run.
      - *see, e.g., here:
        - [github issue #16](https://github.com/snakemake/snakemake/issues/16)
        - [github issue #55](https://github.com/snakemake/snakemake/issues/55)
2. URLs to data sources
  - The pipeline is designed as a self-contained piece of software that attempts downloading
    all necessary reference and data files from public sources/databases. It has happened
    in the past that unstable URLs lead to pipeline errors. This is obviously uncontrollable
    and requires manual correction of the affected URLs.
