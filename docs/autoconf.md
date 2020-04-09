# Autoconf

## Using the autoconf.py script

The probably easiest way to start a pipeline run for your own data is to use the `autoconf.py` script
to generate the necessary configuration file. Since all pipeline configuration is realized with simple
textual [YAML](https://yaml.org/) files, you can edit the auto-generated configuration files if you
need more flexibility. Please note that the `autoconf.py` script only supports generating config files
for one sample at a time.

Run the `autoconf.py` script as follows:

```bash
/work_dir$ conda activate ./smk_env
(smk_env)/work_dir$ cd project-diploid-assembly
(smk_env)/work_dir/project-diploid-assembly$ ./autoconf.py
```

The script is interactively guiding you through the configuration process by asking a series of
basic questions about your data, e.g., the local storage path or the type of long reads. You can
always accept the default value (if one is provided!) by hitting `<enter>`. You can reduce
the number of questions to the bare minimum by accepting all defaults:

```bash
(smk_env)/work_dir/project-diploid-assembly$ ./autoconf.py --accept-defaults
```

After you successfully completed the autoconf process, you find two additional folders in your working
directory:

```bash
/work_dir$ ls -1
autoconf_config/
autoconf_linked_data/
project-diploid-assembly/
smk_env/
```

The `autoconf_config` folder contains the generated configuration file, and the `autoconf_linked_data`
folder contains symbolic links to the input data. The symbolic links are named following the pattern
required by the pipeline to process your data correctly.

**Caveat**: if the `autoconf.py` script fails at deriving well-behaved names for your input files, please
open a github issue showing a handful of examples of the file names that cannot be processed. However, since
file names are a matter of personal preference, or sometimes of project requirements, the worst case would be
that you have to create appropriately named symbolic links to your input files yourself.

Next, please proceed to the documentation on how to [execute the pipeline](execute.md).