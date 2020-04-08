# Pipeline tutorial: how to run a diploid genome assembly

The following step-by-step instructions describe how to configure and run the diploid genome assembly pipeline.
This tutorial covers three different use cases:
1. running the pipeline on the pre-configured HGSVC samples to reproduce published results
2. running the pipeline on a small [demo](demo.md) dataset ("proof of operability")
3. running the pipeline on locally available data (custom run), either creating the necessary
configuration files manually or via using the `autoconf.py` script.

Since some parts of the pipeline setup are Snakemake-specific, certain sections of this tutorial
describe configurations that are mandatory for all three use cases (sections will be marked as such).
In general, pipeline execution is **only** supported on **Linux systems**.

## Required input data
*(This is given for the pre-configured HGSVC and the demo data. Required for all use cases)*

The pipeline has been tested with PacBio CLR, HiFi, and Oxford Nanopore ultra-long reads. The expected input formats
are as follows:
- PacBio CLR: ["pacbio-native"](https://pacbiofileformats.readthedocs.io/en/5.1/BAM.html) BAM
- PacBio HiFi/CCS: FASTQ (gzipped)
- Oxford Nanopore: FASTQ (gzipped)
- Strand-seq: FASTQ (gzipped)

## Get the pipeline code
*(Always required)*

For the rest of this guide, we assume `work_dir` to be our top-level directory:
```bash
/work_dir$
```

Clone the pipeline git repository and **(TODO: merge with master)** switch to the development branch:

```bash
/work_dir$ git clone https://github.com/ptrebert/project-diploid-assembly.git
/work_dir$ cd project-diploid-assembly
/work_dir/project-diploid-assembly$ git checkout development
```

## Snakemake execution environment
*(Always required)*

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
to configure its own behavior depending on the compute environment. A Snakemake *profile* is simply
a configuration file containing Snakemake command line parameters to avoid retyping them at every
Snakemake pipeline invocation. Please refer to the Snakemake help linked above for a complete list
of possible parameters for a Snakemake *profile*.

#### Execution environment: single server
**Only recommended for testing purposes or for running the demo data**

You can examine an example for a Snakemake *profile* suitable for single-server execution here:

```bash
/work_dir$ less project-diploid-assembly/environment/snakemake/server/d3compute/config.yaml
```

The following options must be set in the profile:

```yaml
cores: <NUM_CPUS>
use-conda: True
resources:
  mem_total_mb=<TOTAL_RAM_IN_MEGABYTE>
```

The option `mem_total_mb` is necessary to avoid running high-memory jobs, e.g., a whole-genome `flye`
assembly of CLR data, if only insufficient resources are available. 

#### Execution environment: compute cluster (recommended)

You can examine an example for a Snakemake *profile* suitable for single-server execution here:

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
on the cluster (typically with a time limit of only a few hours). Currently, pipeline jobs that are
configured for cluster submission using default resources require only one CPU core and less than 2048 MB
of memory.

Additionally, Snakemake supports checking for the status of submitted jobs using a so-called
[*cluster status* script](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status).
You can find an example for a PBS Pro-compatible script here:

```bash
/work_dir$ less project-diploid-assembly/scripts/cluster_status/hhu_hilbert.py
```

**Important reminder: the above configuration file is a Snakemake *profile* (supplied at the Snakemake
command line via `--profile`). All files described in the remainder of this tutorial are "regular"
configuration files (supplied at the Snakemake command line via `--configfiles`)**

## How to proceed...

If you want to **run the demo dataset** with just the default environment configuration,
you can immediately [switch to the respective instructions](demo.md).

If you want to reproduce published results using HGSVC samples, or run the pipeline on your
own data, please keep reading this tutorial. If you want to use the [`autoconf.py`](autoconf.md)
script to generate a suitable configuration file, you should still read the next section ("*Pipeline
run environment configuration*") to better understand resource settings. Afterwards, you
can switch to the [`autoconf.py`](autoconf.md) documentation.

## Pipeline run environment configuration
*(Not required for the demo or if the [`autoconf.py`](autoconf.md) script is used)*

Since the efficiency of the pipeline strongly depends on how many independent jobs can run in parallel,
the pipeline requires an additional Snakemake configuration file containing more info about CPU resources:

```yaml
num_cpu_max: <NUM>
num_cpu_high: <NUM>
num_cpu_medium: <NUM>
num_cpu_low: <NUM>
```

- max cpu: jobs like whole-genome assemblies (usually also requiring a lot of memory) will use this many
CPU cores. This should be the maximal number of CPU cores available in a single machine.
- high cpu: jobs like whole-genome read alignment with minimap2 or pbmm2 will use this many CPU cores.
Ideally, this number should be larger than 20.
- medium cpu: jobs like a QUAST analysis of an assembly will use this many CPU cores. Ideally, this
number should be between 10 and 20.
- low cpu: jobs like Strand-seq alignments will use this many CPU cores. Ideally, this number
should be between 4 and 10.

Additionally, if the compute cluster supports [environment `modules`](http://modules.sourceforge.net), the module
name for the [`Singularity` container runtime](https://sylabs.io/singularity/) has to be specified here
(mandatory for using the Peregrine assembler and the DeepVariant variant caller). Alternatively, you can set this
option to `False`:

```yaml
num_cpu_max: <NUM>
num_cpu_high: <NUM>
num_cpu_medium: <NUM>
num_cpu_low: <NUM>
env_module_singularity: False  # or name of Singularity module
```

You can examine an example for a complete run environment configuration file here:

```bash
/work_dir$ less project-diploid-assembly/smk_config/run_env/smk_cfg_env-hhu.yml
```

## Pipeline parameter and reference data configuration

The current development release of the pipeline supports only human data
(either GRCh37 or GRCh38 genome reference used for evaluation) out of the box.
The following configuration files can thus be used as-is:

```bash
/work_dir$ less project-diploid-assembly/smk_config/params/smk_cfg_params_RV9.yml
/work_dir$ less project-diploid-assembly/smk_config/ref_data/reference_data_sources.yml
```

TODO: continue from here

## Pipeline sample configuration
Configuring sample data for the pipeline can be logically split into three parts: (i) specifying the sample
and the types of associated input read sets; (ii) specifying which read sets to use for which step of the pipeline;
(iii) specifying the data source for each read set (e.g., local or FTP).

Note that, internally, Snakemake merges all configuration files that are supplied via `--configfiles`. The three
sections mentioned above are described separately, but can be put into the same configuration file if it seems
reasonable to do so. (Note: to be correct, also all of the above configuration except for the
Snakemake profile could be placed in the same file. This of course interferes with reusing some parts of
the configuration for other pipeline runs).

#### (i) Sample and read set specification

This section describes which read sets are associated with the same sample (there are a few places in the pipeline
where it is checked that only read sets from the same individual are processed together).

A minimal sample configuration is structured as follows (using an HGSVC sample as example here):

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

> You can find the original file here:
> ```bash
> /work_dir$ less project-diploid-assembly/smk_config/samples/AFR/ACB/hg02011.yml
> ```


The entry `individual` must match the sample name at the end of `sample_description_`.
The actual values for `super_population`, `population`, and `family` are irrelevant, and are only
used to sort the output for this sample together with other potential family members.  
The `data_sources` entry must contain at least one of each `long_reads` and `strandseq`, but is otherwise
unlimited, i.e., many different read sets for the same individual can be configured at once.

For **long-read** read sets, the name has to be of the form `individual_project_platform`.
In the above example, that is:
- individual: HG02011
- project: hgsvc
- platform: pbsq2-clr (PacBio Sequel2, CLR data)

Currently supported platforms are:
- pb (for Pacbio) with any suffix (sq2, sq1 etc.) for the read types `-clr` and `-ccs`
- ont (for Oxford Nanopore) with any suffix for the read types `-ul` or `-any`.

The remaining information for long-read read sets basically determines the (pre-) processing steps:

- technology: `pacbio` or `nanopore`; PacBio reads will always be aligned with `pbmm2` instead of `minimap2`
- data_type: `pacbio_native` or `fastq`; only PacBio-native input reads can be used for arrow polishing. Note that
"pacbio-native" BAM files are identified by the extension ".pbn.bam", whereas other BAM files have the extension
".sam.bam" (this is mandatory in the pipeline context).
- load_type: `parts` or `complete`; input data coming in parts, e.g., one BAM file per SMRT cell, will be merged
to simplify downstream processing

For **Strand-seq** read sets, the same rules as above apply, except for the required suffix `sseq`. The platform
information for Strand-seq reads is (currently) not used in the pipeline, and is thus unrestricted. For the above
read set `HG02011_hgsvc_ilnxs-80pe_sseq`, the platform specification indicates 80bp paired-end reads generated on
an Illumina NextSeq machine.  
Since older Strand-seq data may have been prepared as two libraries that need to be merged into one, the pipeline
supports both current mono-fraction samples (`library_fractions: one`) and older libraries
(`library_fractions: two`).

#### (ii) Sample target specification
The pipeline is configured to produce a number of `targets` (a term borrowed from Snakemake), where each target
simply represents an output file, e.g., the phased assemblies or the variant calls. The current set of all
preconfigured target specifications can be examined in:

```bash
/work_dir$ less project-diploid-assembly/smk_include/targets.smk
```

Since Snakemake is a filename-driven workflow engine, it is mandatory to specify the full path of the requested
output file to trigger its production. Given the large number of useful output files produced by a single
pipeline run, this would be cumbersome. Hence, it is sufficient to specify the
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

> You can find the original file here:
> ```bash
> /work_dir$ less project-diploid-assembly/smk_config/samples/AFR/ACB/hg02011.yml
> ```

The sample target section of the config is composed of two entry types: `defaults` and `target`. Default values
are used repeatedly for each following target section (to avoid repetition). The above configures the pipeline
to produce two "target sets" for the same input reads: one target set uses the `flye`
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
  - supported values: flye, wtdbg, pereg (= Peregrine), shasta
  - for expert users: any value matched by `[a-z0-9]+`
- hap_assembler: tool for haploid assembly
  - supported values: flye, wtdbg, pereg (= Peregrine), shasta, canu
- var_caller: tool for variant calling
  - supported values: longshot, freebayes, deepvar (= DeepVariant)

> **Important note**: if you happen to forget (or omit) a wildcard, the respective target cannot be produced
> by the pipeline. There won't be a message about this unless you set `show_warnings: True` in one of your
> config files (e.g., in the pipeline run environment configuration). Consider the
> case of omitting the wildcard "pol_pass: arrow-p1": no targets requiring a polished phased assembly can be
> build w/o that information, but all other targets will be created as usual. This is one possible way of
> stopping the pipeline early, e.g., to examine the draft haploid assembly before proceeding with the
> polishing step.

#### (iii) Data source configuration
The last piece of information pertains to the question where the pipeline should retrieve its input data from.
For all preconfigured HGSVC samples, this is currently either an FTP server or EBI/ENA.
The following descriptions thus focuses on locally available data (both long read and Strand-seq input files).

##### (iii-A) Long-read input data
We assume you have a local folder hierarchy where you collect the output data for all your sequencing experiments.
Because the output file names are meaningful for your internal file tracking, you cannot rename the original files.

For one single sample, it looks like this:

```bash
/seq_experiments/HG02011/HG02011_EDEVI_20200207_S64049_PL100149417-1_A01.subreads.bam
/seq_experiments/HG02011/HG02011_EDEVI_20200211_S64049_PL100149417-1_A01.subreads.bam
```

Because you cannot rename the files, you simply create symbolic links with appropriate names:

```bash
/linked_experiments/HG02011_hgsvc_pbsq2-clr/HG02011_hgsvc_pbsq2-clr.part1.pbn.bam
/linked_experiments/HG02011_hgsvc_pbsq2-clr/HG02011_hgsvc_pbsq2-clr.part2.pbn.bam
```

> Reminder: "pbn" is short for "pacbio-native" and is used within the pipeline context
> to separate "pacbio-native" from "non pacbio-native" BAM files (with file extension ".sam.bam")

Since the above symbolic links (names) are supposed to be used as-is by the pipeline, and the input data is coming
in parts (as indicated in the sample description, see part (i) above), it is necessary to
create the subfolder `HG02011_hgsvc_pbsq2-clr` where the individual `part` files will be placed.

Now you can configure the data source as follows:

```yaml
data_source_HG02011_local:
  comment: "OPTIONAL: annotate your data source"
  output: 'HG02011_local.json'
  server: 'localhost'
  data_source: '/linked_experiments'
  collect_files:
    - 'pbn.bam'
  sort_into:
    - 'bam'
  assume_correct_filenames: True
```

> **Important remark**: there is no need to have one data source configuration per sample. If there is a collection
> of samples locally available, they can all be configured in the same data source
> configuration file by placing a whole hierarchy of symbolic links under "/linked_experiments"

"collect_files" specifies which files to collect based on their file extension. "sort_into" tells the pipeline
where to put these files (here: under `input/bam`, plus the subfolder `HG02011_hgsvc_pbsq2-clr`). If the input data
were PacBio CCS/HiFi or Oxford Nanopore reads with the input format FASTQ, the above would change to `fastq.gz` and
`fastq`. Note that you can have BAM and FASTQ mixed in the same folder, there just has to be a one-to-one
correspondence between the entries in "collect_files" and "sort_into".
The import option here is `assume_correct_filenames: True`, which tells the pipeline to skip any attempt to
infer appropriate file names.

##### (iii-B) Strand-seq input data
In principle, the configuration of a local data source for Strand-seq data is quite similar to the above.
The only change required is that, because the individual Strand-seq FASTQ files are never merged into a single
FASTQ file for downstream processing, and Strand-seq data are assumed to be paired-end (short) reads, the file naming
scheme is slightly different:

```bash
/linked_experiments/HG02011_hgsvc_ilnxs-80pe_sseq/HG02011_hgsvc_ilnxs-80pe_RUN-ID_1.fastq.gz
/linked_experiments/HG02011_hgsvc_ilnxs-80pe_sseq/HG02011_hgsvc_ilnxs-80pe_RUN-ID_2.fastq.gz
(and so one)
```

> **Important note**: the above is only valid for mono-fraction Strand-seq libraries.

Since unique file names are required, there has to be a unique ID in each file name after the usual
`individual_project_platform` part, but before the mate indicators (`_1` and `_2`), i.e., what is indicated
above as "RUN-ID".  
Also note the suffix `sseq` at the end of the subfolder. It is noteworthy that both long-read and Strand-seq input
data can be collected (symlinked) in the same folder hierarchy provided that their file names can be used as-is

#### Q and A concerning data input

##### QA-1: how do I know which files were the original input for my pipeline?
Internally, the pipeline uses so-called *request* files to keep track of each original file source (remote or local).
These request files can be found in `input/FORMAT/READSET/requests`.

##### QA-2: my file names are quite well-behaved, do I have to create intermediate symbolic links?
You can try playing around with the data scraping script located here:
```bash
/work_dir$ less project-diploid-assembly/scripts/scan_remote_path.py
```
This script is for internal use only (mainly to rename the HGSVC data located on the 1000G FTP server).
To give an example, the long-read input files listed under (iii-A) could probably be renamed automatically
using the following call to the script:

```bash
scan_remote_path.py --debug \
    --server localhost \
    --data-source /seq_experiments \
    --collect-files "bam" \
    --sort-into "bam" \
    --assume-pacbio-native \
    --assume-clr-subreads \
    --file-infix "hgsvc-pbsq2-" \
    --local-path-suffix "{individual}_{file_infix}{tech}" \
    --output hg02011_local.json
```

If the output file names (in the JSON) adhere to the naming scheme as required by the pipeline, then the data
source can also be specified as follows (i.e., avoiding the intermediate step of symlinking with appropriate names):

```yaml
data_source_HG02011_local:
  comment: "OPTIONAL: annotate your data source"
  output: 'HG02011_local.json'
  server: 'localhost'
  data_source: '/seq_experiments'
  collect_files:
    - 'bam'
  sort_into:
    - 'bam'
  assume_pacbio_native: True
  assume_clr_subreads: True
  file_infix: "hgsvc_pbsq2-"
  local_path_suffix:  "{{individual}}_{{file_infix}}{{tech}}"
```

##### QA-3: are locally available input files copied?
Generally no. If symbolic links turn out to be problematic for some reason, placing a
`force_local_copy: True` in any of the Snakemake config files will trigger copying input files.

## Running the pipeline

#### Step 1 (optional)
Snakemake supports creating Conda environments on the fly to isolate software installations. This feature is used
in the pipeline, which leads to some overhead when the individual environments are created. Since software setup is
also a common point of failure, the whole Conda environment creation can be done before executing the pipeline.

This can be achieved as follows (always perform a dry run first):

```bash
/work_dir$ conda activate ./smk_env
(smk_env)/work_dir$ cd project-diploid-assembly
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile path_to_your_profile/ \
    --configfiles path_to_your_run_env/run_env.yml smk_config/params/smk_cfg_params_RV8.yml \
    --cluster-status path_to_your_cluster_status_script.py \
    setup_env
```

Note that the environment setup also checks the presence of the Singularity container runtime only needed for
the Peregrine assembler and the DeepVariant variant caller. If these two tools are not to be used for any
pipeline run, a failed check can be ignored. 

#### Step 2
After a successful environment setup, a pipeline run can be started as follows (again, do a dry run first):

```bash
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile path_to_your_profile/ \
    --configfiles path_to_your_run_env/run_env.yml \
                    smk_config/params/smk_cfg_params_RV8.yml \
                    smk_config/ref_data/reference_data_sources.yml \
                    path_to_your_sample_config.yml \
                    path_to_your_data_source_config.yml \
    --cluster-status path_to_your_cluster_status_script.py \
    master_custom
```

Note that the above assumes that you create two configuration files, one only for the sample description and
target specification (`path_to_your_sample_config.yml`), and a second one to define the data sources
(`path_to_your_data_source_config`). As already mentioned, all config files are merged by Snakemake, so it is
up to you how you organize your configuration files.

The given Snakemake target `master_custom` triggers a pipeline run for all possible targets (=output files) for
all samples that are found in the configuration (i.e., for each `sample_description_` entry with specified targets
and data sources).

#### Step 3
To simplify collecting the most important results from a successful pipeline run, the full (!) paths of all
produced targets are dumped to a single file located here:

```bash
/work_dir/run_folder/output/targets/<SUPER-POPULATION>_<POPULATION>_<FAMILY>/<INDIVIDUAL>.fofn
```

This simplifies copying the most important results to a permanent storage location by reading the file paths
from the above "file of filenames" (fofn).
