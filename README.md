# Project repository: reference-free diploid assembly of a human genome

## Citation

DOI of preprint: [10.1101/855049](https://doi.org/10.1101/855049)

Porubsky and Ebert et al.: "A fully phased accurate assembly of an individual human genome" (*in revision*)

## Scope of this repository

This repository contains the Snakemake pipeline code plus some auxiliary scripts to go from raw
input data to polished haploid assemblies. Any self-contained, general purpose software tool used in
the pipeline is either available via conda/bioconda, or via github. In any case, the pipeline
implementation covers the entire software setup required for a complete pipeline run. 

In particular, the code for the `SaaRclust` R package is
available in [David Porubsky's github](https://github.com/daewoooo/SaaRclust).

## Documentation

There are several step-by-step manuals available that describe all use cases currently supported
for this pipeline. First-time users should start by reading the [tutorial](docs/tutorial.md).
If you encounter any problems or "strange behaviour" during pipeline execution, please check
the [FAQ](docs/faq.md) for explanations and solutions. If this does not help, please open a
[github issue](https://guides.github.com/features/issues).


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
