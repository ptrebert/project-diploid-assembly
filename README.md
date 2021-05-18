# Project repository: Phased Genome Assembly using Strand-seq (PGAS)

## Citation

If you use this pipeline or extract and reuse original code/rules from this repository,
please cite the following two papers:

> Porubsky and Ebert et al.  
> "Fully Phased Human Genome Assembly without Parental Data Using Single-Cell Strand Sequencing and Long Reads."  
> Nature Biotechnology, December 2020  
> [DOI: 10.1038/s41587-020-0719-5](https://doi.org/10.1038/s41587-020-0719-5)

> Ebert, Audano, Zhu and Rodriguez-Martin et al.  
> "Haplotype-resolved diverse human genomes and integrated analysis of structural variation"  
> Science, February 2021  
> [DOI: 10.1126/science.abf7117](https://doi.org/10.1126/science.abf7117)

#### Deprecated citations

Please do not reference the preprints ([10.1101/855049](https://doi.org/10.1101/855049) and [10.1101/2020.12.16.423102](https://doi.org/10.1101/2020.12.16.423102)) anymore.

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


## Changelog

### Current DEV

- Because of various serious problems related to Snakemake's "checkpoints", all checkpoints in the pipeline were removed and turned into regular rules.
This should solve a number of checkpoint-related bugs/problems (e.g., gh#55, gh#216, gh#727 and probably more) that are otherwise difficult to tackle
(upgrading Snakemake is not possible because of gh#883 in case this would fix some issues). However, this change removes dynamic updates from the workflow,
which interferes with dynamic inputs (= Strand-seq libraries that are excluded at runtime after library QC). It is likely that the pipeline requires
at least one restart to fix job failures that result from the discrepancy between the Strand-seq libraries loaded from the data source, and the Strand-seq
libraries that remain available after QC (provided that automatic QC is set for the particular sample). Effective after #ed628ae.
- automatic Strand-seq library QC using ASHLEYS can be activated by adding `library_qc: yes` to the Strand-seq readset config
- long-read support in the pipeline is officially limited to PacBio HiFi/CCS; using other long reads as input for assembly may or may not work
- support for scraping remote data sources has been dropped. At the moment, setting `use_legacy_data_scraping` can be set in a config to reactivate the behavior,
but only for code versions predating the removal of all checkpoints (in #ed628ae). Because of this incompatibility, all code related to remote data scraping / legacy data scraping
will be removed in one of the next commits. If access to the data sources used for the HGSVC2 paper is required, tagged version `v1.0.1` must be used.
- support for local data sources does not yet include Illumina short reads
- in addition to `show_warnings: true`, `show_debug_messages: true` can be added to a config to print info to stderr (mostly relevant for developing)