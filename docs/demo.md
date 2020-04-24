# DEMO

## Running the pipeline demo data

The following instructions assume that you have read the [tutorial](tutorial.md)
at least up to the point "*Snakemake execution environment*". Your working directory
should thus look as follows:

```bash
/work_dir$ ls -1
project-diploid-assembly/
smk_env/
```

Please download the [demo data (DOI: 10.5281/zenodo.3746293)](https://doi.org/10.5281/zenodo.3746293)
into your working directory (~ 6.2 GB), and extract the gzipped tar:

```bash
/work_dir$ tar xzvf pipeline_demo.tar.gz
```

After this operation, your working directory should look like this:

```bash
/work_dir$ ls -1
demo_data/
pipeline_demo.tar.gz
project-diploid-assembly/
smk_env/
```

The pipeline repository contains a Snakemake *profile* that specifies a compute environment
with **24 CPU cores** and **64 GB of main memory**. You can either use the Snakemake *profile* and the
pipeline run environment configuration that are shipped with the pipeline code in the repository,
or you can use your own based on the information given in the [tutorial](tutorial.md).
In both cases, please proceed to the instructions how to [execute the pipeline](execute.md).

## How to interpret the results of the demo

In all brevity, just don't. The demo data is a heavily downsampled version of a publicly
available PacBio Sequel-2 HiFi/CCS dataset retrieved from EBI/ENA (PRJNA540705),
and of the respective Strand-seq data (PRJEB14185). The objective was to create a dataset
that could be processed from start to finish with moderate resources and within a reasonable
amount of time (less than 24 hours). A successful run of the demo data is a "proof of function"
for the pipeline, but it does not generate "biologically interesting" results.