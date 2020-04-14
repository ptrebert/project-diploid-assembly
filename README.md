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