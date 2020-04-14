# EXECUTE

## Prepare the software environment

Strictly speaking, this step is not mandatory as Snakemake will [create the necessary `Conda` environments
on-the-fly](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management)
(**caveat**: `use-conda: True` must be set in the Snakemake *profile*, or at least at the command
line when starting the pipeline).

However, it is highly recommended to create all `Conda` environments and to perform all software installation steps
before starting the actual pipeline run. This can be achieved as follows:

```bash
/work_dir$ conda activate ./smk_env
(smk_env)/work_dir$ cd project-diploid-assembly
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile path_to_your_profile/ \
    --configfiles path_to_your_run_env/run_env.yml smk_config/params/smk_cfg_params_RV9.yml \
    setup_env
```

**Caveat**: the above command performs only a `--dry-run`, and does not trigger any software setup. After
a successful dry run, remove the `--dry-run` part from the command to start the process.

Please note that if your compute infrastructure is a cluster, you can add your cluster status script to
the above call because certain software installations are performed as independent jobs that can be executed
on a compute node in a cluster.

#### Demo

If you are running the **demo** dataset with the default environment configuration, the software setup
command looks like this:

```bash
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile environment/snakemake/demo/ \
    --configfiles smk_config/demo/run_env.yml smk_config/demo/params.yml \
    setup_env
```

Please note the command `setup_env` at the end. If this were omitted, a regular pipeline run would be triggered
which would ungracefully fail because the configuration is incomplete for a regular run.
The software setup can take quite some time depending on the speed of your storage back end.

#### Acceptable failures

The software setup part of the pipeline routinely checks for the availability of the `Singularity`
container runtime. If your infrastructure does not support `Singularity`, a failure of the respective
job can be ignored; this implies, however, that you cannot use tools that are executed in a `Singularity`
container (e.g., the Peregrine assembler, or the DeepVariant variant caller).


## Start the pipeline run

The command to start the actual diploid assembly process follows the same structure as the software setup above.
The following is a summary of all possible configuration files, and specific calls for all three use cases
(see the [tutorial](tutorial.md): reproduce, demo, or custom data) are given afterwards.

Configuration checklist:
1. Snakemake *profile*
2. one pipeline run environment config
3. one pipeline parameter config
4. one or many pipeline sample config(s)
5. pipeline sample target config (as many as sample configs)
6. one or many pipeline data source config(s)
7. one or many pipeline reference data source config(s)

**Caveat**: the pipeline config files 2. to 7. are all supplied at the Snakemake command line via `--configfiles`,
and are read one after the other by Snakemake and merged into a single configuration inside the pipeline.
During this process, configuration parameters will be overwritten **without warning** if they are specified
in multiple configuration files (only the parameter read last from the configuration files will be present).

**Simplification**: since Snakemake merges all configuration files into a single object anyway, you can put all
information into a single file if you wish (which is done by the [`autoconf.py`](autoconf.md) script). As an example,
the sample and the sample target configuration can be reasonably put in the same file (and potentially also
the respective data source). Note that, if your data is easily accessible (public download or similar), you can
"share" your results by simply sharing your pipeline configuration files with your collaboration partners,
provided they have enough computational resources to run the pipeline with your configuration.

Start the pipeline as follows (assumed infrastructure is a compute cluster):

```bash
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile PATH_TO_YOUR_SNAKEMAKE_PROFILE/ \
    --configfiles ALL_PIPELINE_CONFIG_FILES \
    --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT
```

As above, remove the `--dry-run` option to start the pipeline run.

#### Demo

For running the demo dataset with all default values for the environment configuration, run the following
(assumed to be executed on a single server):

```bash
(smk_env)/work_dir/project-diploid-assembly$ snakemake \
    --dry-run \
    --directory ../run_folder \
    --profile environment/snakemake/demo/ \
    --configfiles smk_config/demo/run_env.yml smk_config/demo/params.yml \
                  smk_config/demo/na12878.yml smk_config/ref_data/reference_data_sources.yml
```

Please note that the config file `na12878.yml` contains the parameter settings of the pipeline
configurations 4., 5. and 6. from the above checklist.

#### Autoconf

If you generated your pipeline configuration using the [`autoconf.py`](autoconf.md) script, a pipeline
run can be started as follows (assumed infrastructure is a compute cluster):

 ```bash
 (smk_env)/work_dir/project-diploid-assembly$ snakemake \
     --dry-run \
     --directory ../run_folder \
     --profile PATH_TO_YOUR_SNAKEMAKE_PROFILE \
     --configfiles ../autoconf_config/run_assembly.yml \
     --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT
 ```

#### Reproduce HGSVC results

All PacBio HiFi/CCS or CLR data produced in the HGSVC-consortium context have been preconfigured in the
pipeline to enable straightforward replication of published pipeline results such as phased and
polished assemblies. Note that using the HGSVC configuration presets shipped with the pipeline repository
always triggers the download of the raw data (at the time of writing: from the IGSR/HGSVC FTP).

Running all HGSVC data (3 family trios for PacBio HiFi/CCS, 26 samples for PacBio CLR) can only be
accomplished on a compute cluster (at least a small one for PacBio HiFi data), which is the assumed
infrastructure in the following examples.

##### Example: run all HGSVC samples with PacBio/HiFi or PacBio CLR data available

Note that reproducing the pipeline results for all samples with available PacBio HiFi data requires using
the Peregrine assembler and the DeepVariant variant caller. Both these tools are executed as `Singularity`
containers in the pipeline.

```bash
 (smk_env)/work_dir/project-diploid-assembly$ snakemake \
     --dry-run \
     --directory ../run_folder \
     --profile PATH_TO_YOUR_SNAKEMAKE_PROFILE \
     --configfiles PATH_TO_YOUR_PIPELINE_RUN_ENVIRONMENT_CONFIG \
                    smk_config/params/smk_cfg_params_RV9.yml \
                    smk_config/data_sources/hgsvc_ftp_sources.yml \
                    smk_config/ref_data/reference_data_sources.yml \
                    `ls smk_config/samples/hgsvc/*/*/*.yml` \  # note the backticks here
                    smk_config/selectors/hgsvc_ccs_run.yml \
     --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT
 ```

The line ``` `ls smk_config/samples/hgsvc/*/*/*.yml` ``` (note the backticks) is expanded by the shell
before the Snakemake call is executed; it simply collects all HGSVC sample configuration files. The
configuration file `smk_config/selectors/hgsvc_ccs_run.yml` contains a couple of keywords that select
those sample targets (= Snakemake wildcard values) that match PacBio/HiFi runs, i.e., that define using
the Peregrine assembler and the DeepVariant variant caller (this is advanced configuration of the
pipeline and not important in the context of this tutorial).

In order to switch to selecting only samples with PacBio CLR data available, replace this
configuration file

```bash
smk_config/selectors/hgsvc_ccs_run.yml
```

with this one

```bash
smk_config/selectors/hgsvc_clr_run.yml
```

##### Example: run one specific HGSVC sample with PacBio/HiFi or PacBio CLR data available

Replace the following line from the above example

```bash
`ls smk_config/samples/hgsvc/*/*/*.yml`
```
with the path to the sample configuration file you want to run, e.g.,

```bash
smk_config/samples/hgsvc/AMR/PUR/hg00733.yml
```
Switching between PacBio/HiFi and PacBio/CLR runs works in the same way as above.

##### Example (beta status): run specific HGSVC samples by family

If you know what family your sample belongs to, but you don't want to look up the sample
ID number (e.g., NA19240), you can select individual samples or families as follows:

```bash
 (smk_env)/work_dir/project-diploid-assembly$ snakemake \
     --dry-run \
     --directory ../run_folder \
     --profile PATH_TO_YOUR_SNAKEMAKE_PROFILE \
     --configfiles PATH_TO_YOUR_PIPELINE_RUN_ENVIRONMENT_CONFIG \
                    smk_config/params/smk_cfg_params_RV9.yml \
                    smk_config/data_sources/hgsvc_ftp_sources.yml \
                    smk_config/ref_data/reference_data_sources.yml \
                    `ls smk_config/samples/hgsvc/*/*/*.yml` \ # note: collect all configs
                    smk_config/selectors/hgsvc_ccs_run.yml \ # switch between CCS and CLR as before
     --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT \
     run_yri_child
 ```

The pipeline will automatically collect the correct sample. Since Snakemake supports specifying
multiple targets at once, you can list several family members or whole trios:

```bash
    [omitted]
                    `ls smk_config/samples/hgsvc/*/*/*.yml` \  # note: collect all configs
                    smk_config/selectors/hgsvc_ccs_run.yml \  # switch between CCS and CLR as before
     --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT \
     run_yri_child \
     run_chs_trio \
     run_pur_mother
 ```

Naturally, this approach is of limited use for individuals who are not part of a family and can thus
not be "identified" via their relationships as above.

```bash
    [omitted]
                    `ls smk_config/samples/hgsvc/*/*/*.yml` \  # note: collect all configs
                    smk_config/selectors/hgsvc_clr_run.yml \  # switch between CCS and CLR as before
     --cluster-status PATH_TO_YOUR_CLUSTER_STATUS_SCRIPT \
     run_hg01596_individual
 ```

The only benefit here is that you don't have to specify the one correct sample configuration file,
but rather let the pipeline sort that out for you.

## Collect results

To simplify collecting the most important results from a successful pipeline run, the full (!) paths of all
produced targets are dumped to a single file located here:

```bash
/work_dir/run_folder/output/targets/<SUPER-POPULATION>_<POPULATION>_<FAMILY>/<INDIVIDUAL>.fofn
```

This enables copying the most important results to a permanent storage location by reading the file paths
from the above "file of filenames" (fofn).
