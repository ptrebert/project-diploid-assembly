#!/usr/bin/env python

import os
import shutil


def get_haplotype(file_name):

    if 'h1-un' in file_name or '_h1_' in file_name:
        return 'h1-un'
    elif 'h2-un' in file_name or '_h2_' in file_name:
        return 'h2-un'
    else:
        raise ValueError('Unrecognized haplotype: {}'.format(file_name))


def get_read_info(sample, file_name):
    if sample == 'NA24385':
        return 'hpg_pbsq2-ccs_1000'
    elif sample == 'NA12878':
        return 'giab_pbsq2-ccs_1000'
    else:
        if 'ccs' in file_name.lower():
            return 'hgsvc_pbsq2-ccs_1000'
        elif 'clr' in file_name.lower():
            return 'hgsvc_pbsq2-clr_1000'
        else:
            raise ValueError('Unrecognized read type: {}'.format(file_name))


def get_assembler_info(read_info):
    if 'ccs' in read_info:
        return 'pereg', 'racon-p2'
    elif 'clr' in read_info:
        return 'flye', 'arrow-p1'
    else:
        raise ValueError('Cannot match assembler to reads: {}'.format(read_info))


def get_new_file_ext(file_name):

    if 'not_scaffolded' in file_name.lower():
        return 'bng-unsupported.fasta'
    elif file_name.endswith('.agp'):
        return 'bng-hybrid.agp'
    elif file_name.endswith('.fasta'):
        return 'bng-scaffolds.fasta'
    else:
        raise ValueError('Cannot handle file name: {}'.format(file_name))
    

def build_new_name(file_name):

    sample = file_name.split('_', 1)[0]
    if sample.startswith('GM'):
        if sample == 'GM00864':
            sample = sample.replace('GM', 'HG')
        else:
            sample = sample.replace('GM', 'NA')
    read_info = get_read_info(sample, file_name)
    hap = get_haplotype(file_name)
    assembler, polisher = get_assembler_info(read_info)
    new_file_ext = get_new_file_ext(file_name)

    new_name = '{}_{}-{}.{}.{}.{}'.format(sample, read_info, assembler, hap, polisher, new_file_ext)
    return new_name


target_path = '/gpfs/project/ebertp/projects/rfdga/production/EVAL/run_folder/output/evaluation/scaffolded_assemblies'

for root, dirs, files in os.walk(os.getcwd()):
    if not files:
        continue
    for f in files:
        if not (f.endswith('.fasta') or f.endswith('.agp')):
            continue
        old_path = os.path.join(root, f)
        new_path = os.path.join(target_path, build_new_name(f))
        if os.path.isfile(new_path):
            continue
        shutil.copy(old_path, new_path)

        
