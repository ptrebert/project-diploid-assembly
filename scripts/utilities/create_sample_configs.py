#!/usr/bin/env python3

import sys
import pathlib as pl
import argparse as argp

import pandas as pd
import yaml


CONFIG_TEMPLATE = """
sample_description_{sample}:
    individual: {sample}
    sex: {sex}
    super_population: {super_population}
    population: {population}
    family: {family_id}
    member: {relationship}
    data_sources:
        - long_reads:
            readset: {sample}_LR_READSET
            technology: pacbio
            data_type: fastq
            load_type: parts
            data_source_folder: EMPTY_PATH
        - strandseq:
            readset: {sample}_SSEQ_READSET
            library_fractions: one
            library_qc: yes
            data_source_folder: EMPTY_PATH
        - short_reads:
            readset: {sample}_SHORT_READSET
            load_type: complete
            data_source_folder: EMPTY_PATH
        - ont_reads:
            readset: {sample}_ONT_READSET
            load_type: parts
            data_source_folder: EMPTY_PATH

sample_targets_{sample}:
  - aliases:
      1: &long_reads {sample}_LR_READSET
      2: &strandseq_reads {sample}_SSEQ_READSET
      3: &ont_reads {sample}_ONT_READSET
      4: &short_reads {sample}_SHORT_READSET
  - defaults:
      hap_reads: *long_reads
      vc_reads: *long_reads
      sseq_reads: *strandseq_reads
      hap_assm_mode: split
      hap:
        - h1-un
        - h2-un
  - target:
      nhr_assembler: hifiasm
      hap_assembler: hifiasm
      var_caller: deepvar
"""


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--sample-list',
        '-sl',
        '-s',
        '-l',
        dest='sample_list',
        type=lambda p: pl.Path(p).resolve().absolute(),
        required=True,
        help='Path to list of sample IDs (TXT).'
    )
    parser.add_argument(
        '--sample-metadata',
        '-sm',
        '-t',
        '-m',
        dest='sample_metadata',
        type=lambda p: pl.Path(p).absolute(),
        default=pl.Path(pl.Path(__file__).resolve().absolute().parents[2], 'annotation', 'samples', 'samples.tsv'),
        help='Path to list of sample metadata table (TSV).',
    )
    parser.add_argument(
        '--out-folder',
        '-of',
        '-o',
        dest='out_folder',
        type=lambda p: pl.Path(p).resolve().absolute(),
        required=True,
        help='Path to output folder (will be created if it does not exist).'
    )
    parser.add_argument(
        '--ignore-unknown-samples',
        '-i',
        action='store_true',
        default=False,
        help='Ignore existing sample config YAML files even if the sample is unknown.'
    )
    args = parser.parse_args()
    if not (args.sample_list.exists() and args.sample_metadata.exists()):
        raise ValueError(
            f'At least one of the file paths for the sample listing "{args.sample_list}" '
            f'and the sample metadata "{args.sample_metadata}" is invalid.'
        )
    return args


def check_unknown_samples(samples, metadata):

    missing = set(samples) - set(metadata['sample'].values)
    if missing:
        raise ValueError(f'The following samples are not in the sample metadata table: {sorted(missing)}')
    return


def check_sample_description(sample_in_filename, sample_desc, metadata):

    err_msg = ''

    sample_in_desc = sample_desc['individual']
    if sample_in_filename != sample_desc['individual']:
        err_msg += f'Sample name mismatch (filename): {sample_in_filename} vs {sample_in_desc}\n'
    
    md_info = metadata.loc[metadata['sample'] == sample_in_filename, :]
    md_idx = md_info.index[0]
    config_keys = ['sex', 'population', 'super_population', 'family', 'member']
    md_keys = ['sex', 'population', 'super_population', 'family_id', 'relationship']
    for cfg_key, md_key in zip(config_keys, md_keys):
        # why str?
        # numeric family IDs will be parsed as int by pyyaml
        cfg_value = str(sample_desc[cfg_key])
        md_value = str(md_info.at[md_idx, md_key])
        if cfg_value.lower() == md_value.lower():
            continue
        # one historic exceptions
        if cfg_key == 'member':
            if cfg_value == 'parent' and md_value in ['mother', 'father']:
                continue
            if cfg_value in ['unrelated', 'unspecified'] and md_value in ['unrelated', 'unspecified']:
                continue
        err_msg += f'Sample {sample_in_filename} metadata mismatch for keys {cfg_key} / {md_key}: '
        err_msg += f'SAMPLE CFG [found] {cfg_value} / METADATA [expect] {md_value}\n'

    if err_msg:
        raise ValueError(err_msg)
    return


def check_sample_config_error(out_folder, metadata, ignore_extra_samples):

    config_files = sorted(out_folder.glob('**/*.yml'))
    config_files.extend(sorted(out_folder.glob('**/*.yaml')))
    existing_configs = []
    for cfg in config_files:
        sample_in_filename = cfg.with_suffix('').name.upper()
        existing_configs.append(sample_in_filename)
        if sample_in_filename not in metadata['sample'].values:
            if ignore_extra_samples:
                sys.stderr.write(f'\nWarning: ignoring unknown sample config at path {cfg}\n')
                continue
            raise ValueError(f'Unknown sample config at path {cfg}')

        with cfg.open(mode='rb') as cfg_file:
            sample_config = yaml.load(cfg_file, Loader=yaml.SafeLoader)
        
        # sample description must exist
        sample_desc = sample_config[f'sample_description_{sample_in_filename}']
        check_sample_description(sample_in_filename, sample_desc, metadata)

        # "sample targets" section may not exist, but if it does
        # it must end with the sample name
        sample_target_keys = sorted([k for k in sample_config.keys() if k.startswith('sample_targets')])
        if len(sample_target_keys) > 1:
            raise ValueError(f'More than one "sample targets" section: {sample_target_keys}')
        elif len(sample_target_keys) == 1:
            if not sample_target_keys[0].endswith(sample_in_filename):
                raise ValueError(f'Sample name in "sample targets" section does not match: '
                                 f'{sample_in_filename} vs {sample_target_keys[0]}')
        else:
            # 0: no sample targets section - that is ok
            continue
    return existing_configs


def build_sample_config(sample_name, metadata, out_folder):
    
    md_info = metadata.loc[metadata['sample'] == sample_name, :]
    md_info = dict(*md_info.to_dict('records'))

    sample_config = CONFIG_TEMPLATE.format(**md_info)

    super_pop = md_info['super_population']
    population = md_info['population']
    out_sub_folder = out_folder / super_pop / population
    out_sub_folder.mkdir(parents=True, exist_ok=True)

    out_file = out_sub_folder / sample_name.lower()
    out_file = out_file.with_suffix('.yml')
    with out_file.open(mode='w') as dump:
        _ = dump.write(sample_config)
    
    # make load test
    with out_file.open(mode='rb') as dump:
        _ = yaml.load(dump, Loader=yaml.SafeLoader)

    return
        

def create_sample_configs(samples, metadata, existing_configs, out_folder):

    for s in samples:
        if s in existing_configs:
            continue
        build_sample_config(s, metadata, out_folder)
    return


def main():
    args = parse_command_line()

    samples = args.sample_list.read_text(encoding='ascii', errors='strict').strip().split()
    metadata = pd.read_csv(
        args.sample_metadata,
        sep='\t',
        header=0,
        index_col=None
    )
    check_unknown_samples(samples, metadata)

    args.out_folder.mkdir(parents=True, exist_ok=True)
    existing_configs = check_sample_config_error(args.out_folder, metadata, False)

    create_sample_configs(samples, metadata, existing_configs, args.out_folder)




if __name__ == '__main__':
    main()
