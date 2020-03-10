# Pipeline tutorial: how to run a diploid genome assembly

The following step-by-step instructions describe how to configure and run the diploid genome assembly pipeline
either on new/custom data, or to reproduce the HGSVC assemblies. This tutorial covers all generic configuration steps,
and should thus also be read before running the pipeline [demo](demo.md) dataset. Pipeline execution is **only**
supported on **Linux systems**.

## Required input data
The pipeline has been tested with PacBio CLR, HiFi, and Oxford Nanopore ultra-long reads. The expected input formats
are as follows:
- PacBio CLR: ["pacbio-native"](https://pacbiofileformats.readthedocs.io/en/5.1/BAM.html) BAM
- PacBio HiFi/CCS: FASTQ
- Oxford Nanopore: FASTQ
- Strand-seq: FASTQ

## Get the pipeline code

For the rest of this guide, we assume `work_dir` to be our top-level directory:
```bash
/work_dir$
```

Clone the pipeline git repository and **(TODO)** switch to the development branch:

```bash
/work_dir$ git clone https://github.com/ptrebert/project-diploid-assembly.git
/work_dir$ cd project-diploid-assembly
/work_dir/project-diploid-assembly$ git checkout development
```

## Snakemake environment configuration

Running the pipeline requires [`Conda`](https://docs.conda.io/en/latest/miniconda.html) and a working
[`Snakemake`](https://snakemake.readthedocs.io/en/stable/) installation. For convenience, a suitable
Conda environment is shipped with the pipeline code and can be created as follows:

```bash
/work_dir$ conda env create \
    --file project-diploid-assembly/environment/conda/conda_dipassm.yml \
    --prefix ./smk_env
```

After successful setup, the Conda environment can be activated as follows:

```bash
/work_dir$ conda activate ./smk_env
```

Snakemake uses the concept of a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
to configure its own behavior depending on the compute environment.

#### Single server
**Only recommended for testing purposes or to run the demo data**: please refer to the Snakemake documentation
for all possible configurations. You can find an example of a single-server profile here:

```bash
/work_dir$ less project-diploid-assembly/environment/snakemake/server/d3compute/config.yaml
```

The following options are required to be set in the profile:

```yaml
cores: <NUM_CPUS>
use-conda: True
resources:
  mem_total_mb=<TOTAL_RAM_IN_MEGABYTE>
```

The option `mem_total_mb` is necessary to avoid running high-memory jobs, e.g., a whole-genome `flye`
assembly of CLR data, if only insufficient resources are available. 

#### Compute cluster (recommended)
Please refer to the Snakemake documentation for all options allowed in a Snakemake profile.
You can find an example of a compute cluster profile here:

```bash
/work_dir$ less project-diploid-assembly/environment/snakemake/cluster/hhu_pbs/config.yaml
```

The following options are required to be set in the profile:

```yaml
cluster: <SUBMIT_COMMAND>
cluster-config: <CLUSTER_CONFIG_JSON>
local-cores: <NUM_LOCAL_CPU>
jobs: <NUM_SUBMIT_JOBS>
use-conda: True
default-resources:
  - mem_per_cpu_mb=<RAM_PER_CPU_MB>
  - mem_total_mb=<TOTAL_RAM_PER_JOB_MB>
  - runtime_hrs=<RUNTIME_HOURS>
  - runtime_min=<RUNTIME_MINUTES>
```

The default resources should be chosen such that small jobs can be scheduled to the fastest queue available
on the cluster (typically with a walltime limit of only two or four hours). Currently, jobs that are
configured for cluster submission using only default resources require only one CPU core and less than 2048 MB
of memory.

Additionally, Snakemake supports checking for the status of submitted jobs using a so-called
[*cluster status* script](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status).
You can find an example for a PBS Pro-compatible script here:

```bash
/work_dir$ less project-diploid-assembly/scripts/cluster_status/hhu_hilbert.py
```

## Pipeline environment configuration
Since the efficiency of the pipeline strongly depends on how many independent jobs can be run in parallel,
and this is determined by the the type of machines available in a compute cluster, the pipeline requires an
additional Snakemake YAML configuration file with the following entries:

```yaml
num_cpu_max: <NUM>
num_cpu_high: <NUM>
num_cpu_medium: <NUM>
num_cpu_low: <NUM>
```

- max cpu: jobs like whole-genome assemblies (usually also requiring a lot of memory)
- high cpu: jobs like whole-genome read alignment with minimap2 or pbmm2
- medium cpu: jobs like a QUAST analysis of an assembly
- low cpu: jobs like Strand-seq alignments

Additionally, if the compute cluster supports [environment `modules`](http://modules.sourceforge.net), the module
name for the [`Singularity` container runtime](https://sylabs.io/singularity/) has to be specified here
(mandatory for using the Peregrine assembler and the DeepVariant variant caller). An example for a complete
environment configuration file can be found here:

```bash
/work_dir$ less project-diploid-assembly/smk_config/run_env/smk_cfg_env-hhu.yml
```

## Pipeline parameter and reference data configuration

**TODO** For human data (and either GRCh37 or GRCh38 for evaluation), the following files can be used as-is:

```bash
/work_dir$ less project-diploid-assembly/smk_config/params/smk_cfg_params_RV8.yml
/work_dir$ less project-diploid-assembly/smk_config/ref_data/reference_data_sources.yml
```

## Pipeline sample configuration

#### Sample and read set specification
Configuring sample data for the pipeline can be logically split into two parts, i.e., (i) specifying the sample
and what output to produce for it, and (ii) specifying where the input data are supposed to come from (the sources
for the long reads and the Strand-seq data).

Note that, internally, Snakemake merges all configuration files that are supplied via `--configfiles`, so the
following can all be put in a single YAML file (to be correct, also all of the above configuration except for the
Snakemake profile could be copied into the same file. This of course interferes with reusing some parts of
the configuration for other pipeline runs).

A minimal sample configuration YAML is structured as follows (using an HGSVC sample as example here):

```yaml
sample_description_HG02011:
  individual: HG02011
  super_population: AFR
  population: ACB
  family: BB13
  data_sources:
    - long_reads:
        readset: HG02011_hgsvc_pbsq2-clr
        technology: pacbio
        data_type: pacbio_native
        load_type: parts
    - strandseq:
        readset: HG02011_hgsvc_ilnxs-80pe_sseq
        library_fractions: one
```

The entry `individual` must match the sample name at the end of `sample_description_`.
The actual values for `super_population`, `population`, and `family` are irrelevant, and are only
used to sort the output for this sample.  
The `data_sources` entry must contain at least one of each `long_reads` and `strandseq`, but is otherwise
unlimited, i.e., many different read sets for the same individual can be configured at once.

For **long-read** read sets, the name has to be specified in the form `individual_project_platform`. In the above example,
that is:
- individual: HG02011
- project: hgsvc
- platform: pbsq2-clr (PacBio Sequel2, CLR data)

Currently supported platforms are:
- pb (for Pacbio) with any suffix (sq2, sq1 etc.) for the read types `-clr` and `-ccs`
- ont (for Oxford Nanopore) with any suffix for the read types `-ul` or `-any`.

The remaining information for long-read read sets basically determines the initial processing steps:

- technology: `pacbio` or `nanopore`; PacBio will always be aligned with `pbmm2` instead of `minimap2`
- data_type: `pacbio_native` or `fastq`
- load_type: `parts` or `complete`; input data coming in parts, e.g., one BAM file per SMRT cell, will be merged
to simplify downstream processing

For **Strand-seq** read sets, the same rules as above apply, except for the required suffix `sseq`. The platform
information for Strand-seq reads is (currently) not used in the pipeline, and is thus unrestricted. For the above
read set `HG02011_hgsvc_ilnxs-80pe_sseq`, the platform specification indicates 80bp paired-end reads generated on
an Illumina NextSeq machine.  
Since older Strand-seq data may have been prepared as two libraries that need to be merged into one, the pipeline
supports both current mono-fraction samples (`library_fractions: one`) and older libraries
(`library_fractions: two`).

#### Target specification
The pipeline is configured to produce a number of `targets` (a term borrowed from Snakemake), where each target
simply represents an output file, e.g., the phased assemblies or the variant calls. The current set of target
specifications can be examined in:

```bash
/work_dir$ less project-diploid-assembly/smk_include/targets.smk
```

Since Snakemake is a filename-driven workflow engine, it is mandatory to specify the full path of the requested
output file to trigger its production. Given the large number of useful output files produced by a single
pipeline run, this would be cumbersome. Hence, it is sufficient to specify the respective
[wildcard](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) values needed to
automatically build all output file paths. This is achieved as follows:

```yaml
sample_targets_HG02011:
  - defaults:
      hap_reads: HG02011_hgsvc_pbsq2-clr_1000
      vc_reads: HG02011_hgsvc_pbsq2-clr_1000
      sts_reads: HG02011_hgsvc_ilnxs-80pe_sseq
      pol_reads: HG02011_hgsvc_pbsq2-clr_1000
      pol_pass: arrow-p1
      hap_assm_mode: split
      hap:
        - h1-un
        - h2-un
        - h1
        - h2
  - target:
      nhr_assembler: flye
      hap_assembler: flye
      var_caller: longshot
  - target:
      nhr_assembler: shasta
      hap_assembler: shasta
      var_caller: freebayes
```

The sample target section of the config is composed of two entry types: `defaults` and `target`. Default values
are used repeatedly for each following target section (to avoid repetition). The above configures the pipeline
to produce two "target sets" for the same input reads (listed under defaults), one target set uses the `flye`
assembler and the `longshot` variant caller, and the other uses the `shasta` assembler and the `freebayes`
variant caller (whether or not that is a reasonable choice is beyond this tutorial).

> Note on "HG02011_hgsvc_pbsq2-clr" vs "HG02011_hgsvc_pbsq2-clr_**1000**": initially, the pipeline was
> designed to potentially also produce downsampled versions of the input long-read data sets. The suffix *1000*
> indicates that 100.0% (=1000) of input reads are to be used. The downsampling feature has not been
> implemented (up to now), but the suffix *1000* is still required in the target specification. 

There is no limit on the number of `defaults` and `target` sections per sample target specification.
The abbreviations refer to the following:

- hap_reads: (input) read set to be used for haplotyping/phasing
- vc_reads: read set used for variant calling
- sts_reads: Strand-seq read set
- pol_reads: reads used for assembly polishing
- pol_pass: number of passes for polisher, currently supported: `racon-p1`, `racon-p2`, and `arrow-p1`
- hap_assm_mode: haploid assembly mode, either per cluster/haplotype (`split`) or only per haplotype (`joint`)
- hap: haplotypes to assemble, four possible combinations as listed above
- nhr_assembler: tool for non-haplotype resolved (collapsed) assembly
- hap_assembler: tool for haploid assembly
- var_caller: tool for variant calling

#### Data source configuration