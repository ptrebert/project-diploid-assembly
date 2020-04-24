# FAQ

#### QA-1: how do I know which files were the original input for my pipeline?
Internally, the pipeline uses so-called *request* files to keep track of each original file source (remote or local).
These request files can be found in `input/FORMAT/READSET/requests`.

#### QA-2: my file names are quite well-behaved, do I have to create intermediate symbolic links?
You can try playing around with the data scraping script located here:
```bash
/work_dir$ less project-diploid-assembly/scripts/scan_remote_path.py
```
This script is for internal use only (mainly to rename the HGSVC data located on the 1000G FTP server).
To give an example, the long-read input files described in the [tutorial](tutorial.md) (*Pipeline sample
configuration*, point (III-A)) could probably be renamed automatically
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

#### QA-3: are locally available input files copied?
Generally no. If symbolic links turn out to be problematic for some reason, placing a
`force_local_copy: True` in any of the Snakemake config files will trigger copying input files.
**Important**: symbolic links do not work when using the Peregrine assembler, you have to copy
the long read input data in this case.

#### QA-4: where do I find information about failed jobs?
If the general Snakemake log file does not contain any info, you can find (for most of the jobs)
a log file under the identical path as the output file for that job in the top-level subfolder
`/work_dir/run_folder/log`. In a cluster setup, the default configuration assumes that the
`stdout` and `stderr` outputs of the jobs are directed to `/work_dir/run_folder/log/cluster_jobs/stdout`
and `*/stderr`, respectively.

#### QA-5: what to do if there is no obvious reason for a failed job?
Restart the pipeline. Given that a complete pipeline run encompasses several thousand jobs, there is
almost always a few that fail (most commonly due to I/O issues). If the same job keeps failing, set
`keep-incompete: True` in your Snakemake profile to prevent Snakemake from deleting any partial/incomplete
output. If there is no helpful evidence in the partial output, please contact the developer of the respective
tool. If you need to know which version of the tool is running in the pipeline, check the `conda`
environment files located here: `/work_dir/project-diploid-assembly/environment/conda/*.yml`

#### QA-6: what does the message "Caught Snakemake error #55: checkpoint evaluation failed" mean?
The pipeline heavily uses
[Snakemake's `checkpoint` feature](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution)
for data-dependent execution of subsequent rules. Unfortunately, this more or less
recent feature is not yet as mature as the older parts of Snakemake's code.
A start-to-finish run of the pipeline may encounter a failed checkpoint evaluation*.
The pipeline includes a workaround for this situation, but it will always display the above
message to indicate (mainly intended for the pipeline developer) that this Snakemake
error still exists. As soon as the Snakemake error has been fixed, several parts of the
pipeline, also causing a lot of overhead at the moment, can be simplified.
  - *see, e.g., here
    - [github issue #16](https://github.com/snakemake/snakemake/issues/16)
    - [github issue #55](https://github.com/snakemake/snakemake/issues/55)


#### QA-7: downloading EBI/ENA metadata (file reports) fails for unclear reasons
It sometimes happens that the URL for downloading metadata (file reports) from EBI/ENA
(temporarily) changes. You can find the URL that is used by the pipeline in this
Python script:

```bash
/work_dir$ less project-diploid-assembly/scripts/utilities/downloader.py
```

Look for the entry `ENA_FILE_REPORT_URL` at the beginning of the file. You can compare
this URL to the one displayed in your browser when you try to download the metadata
(file report) manually. In case it differs, please open a github issue.

#### QA-8: can I use precomputed results in the pipeline?
In general, yes, but this requires expert knowledge in using Snakemake and most often
some manual tweaking, e.g., to get around the various checkpoints in the pipeline
(see QA-6 above). One example that is comparatively straightforward can be discussed
in more detail here. Let's assume you want to use an existing genome assembly as your
"collapsed"/non-haplotype resolved (`nhr_assembly` in the pipeline context) assembly.
You can realize that under the following conditions:
1. you name the FASTA assembly file according to pipeline requirements
2. you specify the correct input and output paths for symlinking the FASTA file
3. you adapt your sample target configuration (see the [tutorial](tutorial.md), section
"Pipeline sample configuration", subsection "(ii) Sample target specification") such that
your external assembly is recognized

Point 1: default collapsed assemblies are named as follows (going by example):
`HG00732_hgsvc_pbsq2-clr_1000_nhr-flye.fasta`. Name your custom assembly in the same
way, but **do not use** assembler names as supported by the pipeline, e.g., `flye`,
`shasta`, or `pereg` (for Peregrine). Why not? If you use one of the default assembler
names, the pipeline will always look for the other output that should have been
produced (also indicating a successful assembly). Also make sure that the part of
the assembly file name indicating the readset (`HG00732_hgsvc_pbsq2-clr_1000`) matches
with data that you actually have available.

Point 2: create a new pipeline configuration file (YAML format) containing the following
info about how to link external files in the pipeline context:

```yaml
link_data_input:
  - /LOCAL/PATH/TO/YOUR/ASSEMBLY/HG00732_hgsvc_pbsq2-clr_1000_nhr-uw27b.fasta

link_data_output:
  - output/reference_assembly/non-hap-res/HG00732_hgsvc_pbsq2-clr_1000_nhr-uw27b.fasta
```

Note here that the non-default assembler name `uw27b` is used (also, assembler names
should be a single word, letters and numbers only). The above config tells the pipeline
which input file (your custom assembly) to link in the right place for the pipeline;
in this case, "the right place" is the folder where collapsed assemblies are expected
to be found. In principle, most (all?) of the pipeline results can be replaced by
externally linked data, but listing all corresponding pipeline paths is of course not
possible here.  
**Make sure that you also use this additional config file when you start the pipeline**.

Point 3: in your sample target specification, there is a Snakemake wildcard called
`nhr_assembler` that tells the pipeline which assembler to use for building the
collapsed assembly. Find that wildcard and replace the value (say, `nhr_assembler: flye`)
with your custom assembler name: `nhr_assembler: uw27b`  
When Snakemake is looking for the rule to produce your assembly
(`HG00732_hgsvc_pbsq2-clr_1000_nhr-uw27b.fasta`), it will find the rule that triggers
the symlinking of your assembly, which is sufficient for Snakemake to proceed with
executing the pipeline using your custom assembly.  
If you want to check that the symlinking worked as expected before starting the complete
pipeline run, start your pipeline (see the [execution](execute.md) part of the tutorial),
and state the following rule name as Snakemake target: `master_link_data_sources`. This
should only trigger linking your precomputed results to the right place in the pipeline
folder hierarchy.
