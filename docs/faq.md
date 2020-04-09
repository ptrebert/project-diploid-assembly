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
`keep-incompete: True` in your Snakemake profile to prevent that Snakemake deletes any partial/incomplete
output. If there is no helpful evidence in the partial output, please contact the developer of the respective
tool. If you need to know which version of the tool is running in the pipeline, check the `conda`
environment files located here: `/work_dir/project-diploid-assembly/environment/conda/*.yml`

