# Project repository: Phased Genome Assembly using Strand-seq (PGAS)

## Citation

If you use this pipeline or extract and reuse original code/rules from this repository,
please cite the following two papers:

> Porubsky and Ebert et al.  
> "Fully Phased Human Genome Assembly without Parental Data Using Single-Cell Strand Sequencing and Long Reads."  
> Nature Biotechnology, December 2020  
> [DOI:10.1038/s41587-020-0719-5](https://doi.org/10.1038/s41587-020-0719-5)

> Ebert, Audano, Zhu and Rodriguez-Martin et al.  
> "De novo assembly of 64 haplotype-resolved human genomes of diverse ancestry and integrated analysis of structural variation"  
> bioRxiv preprint (*article in press*)  
> [DOI:10.1101/2020.12.16.423102](https://doi.org/10.1101/2020.12.16.423102)

#### Deprecated citation

Please do not reference the preprint ([10.1101/855049](https://doi.org/10.1101/855049)) of the *Nature Biotechnology* paper anymore.

## Scope of this repository

This repository contains the Snakemake pipeline code plus some auxiliary scripts to go from raw
input data to polished haploid assemblies. Any self-contained, general purpose software tool used in
the pipeline is either available via conda/bioconda, or via github. In any case, the pipeline
implementation covers the entire software setup required for a complete pipeline run. 

In particular, the code for the `SaaRclust`, `StrandPhaseR` and `breakpointR` R packages is
available in [David Porubsky's github](https://github.com/daewoooo/SaaRclust).

## Documentation

There are several step-by-step manuals available that describe all use cases currently supported
for this pipeline. First-time users should start by reading the [tutorial](docs/tutorial.md).
If you encounter any problems or "strange behaviour" during pipeline execution, please check
the [FAQ](docs/faq.md) for explanations and solutions. If this does not help, please open a
[github issue](https://guides.github.com/features/issues).