
import os
import sys
import re

from snakemake.exceptions import WildcardError as WildcardError


TARGET_PATHS = {
    "BUILD_NHR_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "non-hap-res",
        "{hap_reads}_nhr-{nhr_assembler}.fasta"
    ),

    "BUILD_CLUSTERED_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "clustered",
        "{sts_reads}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}.fasta",
    ),

    "BUILD_DRAFT_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "draft", "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.fasta"
    ),

    "BUILD_POLISHED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "polishing",
        "{pol_reads}",
        "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.fasta"
    ),

    "REPORT_NHR_ASSEMBLY": os.path.join(
        "output", "evaluation", "quastlg_busco",
        "{eval_known_ref}-{eval_gene_model}",
        "reference_assembly", "non-hap-res",
        "{hap_reads}_nhr-{nhr_assembler}/report.pdf"
    ),

    "REPORT_DRAFT_HAPLOID_ASSEMBLY": os.path.join(
        "output", "evaluation", "quastlg_busco",
        "{eval_known_ref}-{eval_gene_model}",
        "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "draft", "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}/report.pdf"
    ),

    "REPORT_POLISHED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "evaluation", "quastlg_busco",
        "{eval_known_ref}-{eval_gene_model}",
        "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "polishing",
        "{pol_reads}",
        "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/report.pdf"
    ),

    "STATS_SAMPLE_SUMMARY": os.path.join(
        "output", "statistics", "stat_dumps",
        "{hap_reads}.{file_ext}.pck"
    ),

    "STATS_VARIANT_CALLING_INITIAL": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sts_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.vcf.stats"
    ),

    "STATS_VARIANT_CALLING": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sts_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.GQ{filter_vcf_gq}.vcf.stats"
    ),

    "STATS_STRANDPHASER": os.path.join(
        "output", "statistics", "phasing",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "{hap_reads}.spr-phased.stats.tsv"
    ),

    "STATS_INTEGRATIVE_PHASING": os.path.join(
        "output", "statistics", "phasing",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "{hap_reads}.wh-phased.stats.tsv"
    ),

    "STATS_READ_HAPLO_TAGGING": os.path.join(
        "output", "statistics", "tag_split",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "{hap_reads}.tags.{tag_source}.tsv"
    ),

    "PLOT_INPUT_SAMPLE_STATS": os.path.join(
        "output", "plotting", "statistics", "input_reads",
        "{hap_reads}.{file_ext}.stats.pdf"
    )
}


def extract_wildcard_config_parameters():
    """
    :return:
    """
    wildcard_re = '{\w+}'
    global_settings = {}
    config_keys = set(config.keys())

    for target_name, target_path in TARGET_PATHS.items():
        target_wildcards = [m.strip('{}') for m in re.findall(wildcard_re, target_path)]

        shared_keys = set(target_wildcards).intersection(config_keys)
        for sk in shared_keys:
            global_settings[sk] = config[sk]
    return global_settings


CONFIG_TARGETS_GLOBAL_SETTINGS = extract_wildcard_config_parameters()


def extract_selected_targets():
    """
    :return:
    """
    try:
        select_targets = config['select_targets']
    except KeyError:
        select_targets = dict()
    else:
        for key, values in select_targets.items():
            if not isinstance(values, list):
                sys.stderr.write('\nERROR: select_targets have to be specified as a name-to-list mapping.\n')
                sys.stderr.write('Found in config: {} to {} ({}) mapping.\n'.format(key, values, type(values)))
                raise ValueError('Invalid select targets specified (expected list).')
    return set(select_targets.keys()), select_targets


CONFIG_TARGETS_SELECTED_KEYS, CONFIG_TARGETS_SELECTED_VALUES = extract_selected_targets()


def annotate_readset_data_types(sample_desc):
    """
    :param sample_desc:
    :return:
    """
    try:
        data_sources = sample_desc['data_sources']
    except KeyError:
        sys.stderr.write('\nERROR: no data sources defined for sample: {}\n'.format(sample_desc))
        raise ValueError('No data sources defined')

    source_annotation = dict()
    for data_record in data_sources:
        if 'long_reads' not in data_record:
            continue
        this_record = data_record['long_reads']
        readset = this_record['readset']
        dt = this_record['data_type']
        if dt == 'pacbio_native':
            v = {'file_ext': 'pbn.bam',
                 'tag_source': 'pbn'}
            source_annotation[readset] = v
        elif dt == 'fastq':
            v = {'file_ext': 'fastq',
                 'tag_source': 'fq'}
            source_annotation[readset] = v
        else:
            raise ValueError('Unexpected data type: {} / {}'.format(dt, this_record))
    return source_annotation


def define_file_targets(wildcards):
    """
    :param wildcards:
    :return:
    """
    individual = wildcards.individual
    try:
        sample_desc = config['sample_description_' + individual]
        if individual != sample_desc['individual']:
            sys.stderr.write('\nIndividual mismatch: {} vs {}\n'.format(individual, sample_desc['individual']))
            raise ValueError('Sample description individual does not '
                             'match requested individual: {} vs {}'.format(individual, sample_desc['individual']))
    except KeyError as ke:
        sys.stderr.write('\nWARNING: no sample description for individual [sample_description_] {} in config\n'.format(individual))
        return []

    try:
        sample_targets = config['sample_targets_' + individual]
    except KeyError as ke:
        sys.stderr.write('\nNo targets specified for individual [target_specification_] {} in config\n'.format(individual))
        raise ke

    # make a copy here to have global settings for potential debugging output
    target_values = dict(CONFIG_TARGETS_GLOBAL_SETTINGS)
    readset_annotation = annotate_readset_data_types(sample_desc)

    file_targets = []

    keep_path = False
    selected_path = None
    if 'select_target_path' in config:
        keep_path = True
        selected_path = config['select_target_path']

    for target_specification in sample_targets:
        if 'aliases' in target_specification:
            continue
        elif 'defaults' in target_specification:
            target_spec = target_specification['defaults']
            target_values.update(target_spec)
        else:
            # Copy dict with default values for each target spec.
            # Note that target specs are processed in order,
            # so switching to other defaults for a second set of
            # target specs is possible
            target_spec = target_specification['target']
            tmp = dict(target_values)
            tmp.update(target_spec)
            keep_target = True
            for key in set(tmp.keys()).intersection(CONFIG_TARGETS_SELECTED_KEYS):
                current_value = tmp[key]
                if isinstance(current_value, str):
                    keep_target &= current_value in CONFIG_TARGETS_SELECTED_VALUES[key]
                elif isinstance(current_value, list):
                    # we check that only lists are in selected_targets
                    current_values = set(current_value)
                    keep_target &= len(current_values.intersection(CONFIG_TARGETS_SELECTED_VALUES[key])) > 0
                else:
                    raise ValueError('Cannot handle data type of target parameter: '
                                     '{} / {} / {}'.format(key, current_value, type(current_value)))
            if not keep_target:
                continue
            for target_name, target_path in TARGET_PATHS.items():
                if keep_path and selected_path != target_name:
                    continue
                if '{hap_reads}' in target_path:
                    hap_readset = tmp['hap_reads']
                    for readset, annotation in readset_annotation.items():
                        if hap_readset.startswith(readset):
                            # prefix-matching because the actual readset likely carries
                            # additional info such as _1000
                            tmp.update(annotation)
                            break
                try:
                    complete_targets = expand(target_path, **tmp)
                except (KeyError, WildcardError):
                        raise ValueError('Missing parameter values for target {}: '
                                         '(known: {})'.format(target_name, tmp))
                else:
                    for entry in complete_targets:
                        assert len(entry) > 1, 'Define file targets: looks like iterating over ' \
                                               'exploded string instead of list: {}'.format(complete_targets)
                        file_targets.append(entry)

    file_targets = sorted(file_targets)

    return file_targets


rule dump_build_targets:
    input:
         define_file_targets
    output:
        'output/targets/{super_population}_{population}_{family}/{individual}.fofn'
    run:
        import os

        root_dir = os.getcwd()

        # The internals of InputFiles are still unclear to me.
        # Looks like the formerly used attribute "(input.)targets" did only
        # exist if there were targets, i.e., the list had to be non-empty
        # for the attribute to exist. Weird behavior...
        if input and len(input) > 0:
            with open(output[0], 'w') as dump:
                for file_target in input:
                    full_path = os.path.join(root_dir, file_target)
                    assert os.path.isfile(full_path), \
                        'Non-existing file as build target: {} / {}'.format(full_path, output[0])
                    _ = dump.write(full_path + '\n')
