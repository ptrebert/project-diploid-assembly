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

TODO

## Collect results

To simplify collecting the most important results from a successful pipeline run, the full (!) paths of all
produced targets are dumped to a single file located here:

```bash
/work_dir/run_folder/output/targets/<SUPER-POPULATION>_<POPULATION>_<FAMILY>/<INDIVIDUAL>.fofn
```

This enables copying the most important results to a permanent storage location by reading the file paths
from the above "file of filenames" (fofn).
