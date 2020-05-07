
import os
import sys
import re

from snakemake.exceptions import WildcardError as WildcardError


TARGET_PATHS = {
    "INIT_COMPLETE_DATA": os.path.join(
        "input",
        "{input_format}",
        "{hap_reads}.{file_ext}{ext_modifier}"
    ),
    "BUILD_NHR_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "non-hap-res",
        "{hap_reads}_nhr-{nhr_assembler}.fasta"
    ),

    "BUILD_CLUSTERED_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "clustered",
        "{sseq_reads}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}.fasta",
    ),

    "BUILD_DRAFT_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "draft", "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.fasta"
    ),

    "BUILD_POLISHED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.fasta"
    ),

    "BUILD_POLISHED_CLUSTERED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "clustering",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.scV{git_commit_version}.fasta"
    ),

    "HAPLOID_READ_COVERAGE_HAP1": os.path.join(
        "output", "cov_tracks", "hap_reads",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}_map-to_{eval_align_ref}.h1.bigWig"
    ),

    "HAPLOID_READ_COVERAGE_HAP2": os.path.join(
        "output", "cov_tracks", "hap_reads",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}_map-to_{eval_align_ref}.h2.bigWig"
    ),

    "HAPLOID_READ_COVERAGE_UN": os.path.join(
        "output", "cov_tracks", "hap_reads",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}_map-to_{eval_align_ref}.un.bigWig"
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
        "{sseq_reads}",
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
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/report.pdf"
    ),

    "REPORT_POLISHED_CLUSTERED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "evaluation", "quastlg_busco",
        "{eval_known_ref}-{eval_gene_model}",
        "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "clustering",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.scV{git_commit_version}/report.pdf"
    ),

    "STATS_SAMPLE_SUMMARY": os.path.join(
        "output", "statistics", "stat_dumps",
        "{hap_reads}.{file_ext}.pck"
    ),

    "STATS_VARIANT_CALLING_INITIAL": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sseq_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.vcf.stats"
    ),

    "STATS_VARIANT_CALLING": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sseq_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.GQ{filter_vcf_gq}.vcf.stats"
    ),

    "STATS_STRANDPHASER": os.path.join(
        "output", "statistics", "phasing",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}.spr-phased.stats.tsv"
    ),

    "STATS_INTEGRATIVE_PHASING": os.path.join(
        "output", "statistics", "phasing",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}.wh-phased.stats.tsv"
    ),

    "STATS_READ_HAPLO_TAGGING": os.path.join(
        "output", "statistics", "tag_split",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "{hap_reads}.tags.{tag_source}.tsv"
    ),

    "PLOT_INPUT_SAMPLE_STATS": os.path.join(
        "output", "plotting", "statistics", "input_reads",
        "{hap_reads}.{file_ext}.stats.pdf"
    ),

    "PLOT_SAARCLUST_DIAG_ASSEMBLY_CLUSTERING": os.path.join(
        "output", "plotting", "saarclust_diagnostics", "reference_assembly", "clustered",
        "{sseq_reads}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}_map-to_{eval_align_ref}.clustering.pdf",
    ),

    "PLOT_SAARCLUST_DIAG_ASSEMBLY_ORIENTING": os.path.join(
        "output", "plotting", "saarclust_diagnostics", "reference_assembly", "clustered",
        "{sseq_reads}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}_map-to_{eval_align_ref}.orienting.pdf",
    ),

    "PLOT_SAARCLUST_DIAG_HAPLOID_ASSEMBLY_CLUSTERING": os.path.join(
        "output", "plotting", "saarclust_diagnostics", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "clustering",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.scV{git_commit_version}_map-to_{eval_align_ref}.clustering.pdf",
    ),

    "PLOT_SAARCLUST_DIAG_HAPLOID_ASSEMBLY_ORDERING": os.path.join(
        "output", "plotting", "saarclust_diagnostics", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "clustering",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.scV{git_commit_version}_map-to_{eval_align_ref}.ordering.pdf",
    ),

    "PLOT_SAARCLUST_DIAG_HAPLOID_ASSEMBLY_ORIENTING": os.path.join(
        "output", "plotting", "saarclust_diagnostics", "diploid_assembly",
        "strandseq_{hap_assm_mode}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sseq_reads}",
        "polishing",
        "{pol_reads}",
        "clustering",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.scV{git_commit_version}_map-to_{eval_align_ref}.orienting.pdf",
    ),

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


def extract_selected_or_skipped_targets(handle_targets):
    """
    :param handle_targets:
    :return:
    """
    try:
        modify_targets = config[handle_targets]
    except KeyError:
        modify_targets = dict()
    else:
        for key, values in modify_targets.items():
            if not isinstance(values, list):
                sys.stderr.write('\nERROR: {} have to be specified as a name-to-list mapping.\n')
                sys.stderr.write('Found in config: {} to {} ({}) mapping.\n'.format(handle_targets, key, values, type(values)))
                raise ValueError('Invalid {} specified (expected list).'.format(handle_targets))
    return set(modify_targets.keys()), modify_targets


CONFIG_TARGETS_SELECTED_KEYS, CONFIG_TARGETS_SELECTED_VALUES = extract_selected_or_skipped_targets('select_targets')

CONFIG_TARGETS_SKIPPED_KEYS, CONFIG_TARGETS_SKIPPED_VALUES = extract_selected_or_skipped_targets('skip_targets')


def check_target_modifier_match(target_spec, check_keys, check_values, check_keep):
    """
    If targets are selected in the config / by the user, check if this target_spec
    is to be kept or skipped

    :param target_spec:
    :param check_keys:
    :param check_values:
    :param check_keep:
    :return:
    """
    modify_target = check_keep
    shared_keys = set(target_spec.keys()).intersection(check_keys)
    if len(shared_keys) < len(check_keys) and check_keep:
        # for keeping a target spec, every attribute has to be true;
        # if a target spec is missing a required key, it must not
        # be selected
        modify_target = False  # this is a misnomer
        shared_keys = set()  # to skip over loop
    for key in shared_keys:
        current_value = target_spec[key]
        if isinstance(current_value, str):
            if check_keep:
                modify_target &= current_value in check_values[key]
            else:
                modify_target |= current_value in check_values[key]
        elif isinstance(current_value, list):
            # we check that only lists are in selected_targets
            current_values = set(current_value)
            if check_keep:
                modify_target &= len(current_values.intersection(check_values[key])) > 0
            else:
                modify_target |= len(current_values.intersection(check_values[key])) > 0
        else:
            raise ValueError('Cannot handle data type of target parameter: '
                             '{} / {} / {}'.format(key, current_value, type(current_value)))
    return modify_target


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

    individual = sample_desc['individual']
    source_annotation = dict()
    for data_record in data_sources:
        readset_type = list(data_record.keys())[0]
        readset_spec = data_record[readset_type]
        readset_name = readset_spec['readset']
        if not readset_name.startswith(individual):
            raise ValueError('Readset {} does not match with individual: {}'.format(readset_name, individual))
        if readset_type != 'long_reads':
            continue
        dt = readset_spec['data_type']
        if dt == 'pacbio_native':
            v = {
                'input_format': 'bam',
                'file_ext': 'pbn.bam',
                'ext_modifier': '',
                'tag_source': 'pbn'
            }
            source_annotation[readset_name] = v
        elif dt == 'fastq':
            v = {
                'input_format': 'fastq',
                'file_ext': 'fastq',
                'ext_modifier': '.gz',
                'tag_source': 'fq'
            }
            source_annotation[readset_name] = v
        else:
            raise ValueError('Unexpected data type: {} / {}'.format(dt, readset_spec))
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
        if bool(config.get('show_warnings', False)):
            sys.stderr.write('\nWARNING: no sample description for individual [sample_description_] {} in config\n'.format(individual))
        return []

    try:
        sample_targets = config['sample_targets_' + individual]
    except KeyError as ke:
        sys.stderr.write('\nNo targets specified for individual [target_specification_] {} in config\n'.format(individual))
        raise ke

    # make a copy here to have global settings for potential debugging output
    global_settings = dict(CONFIG_TARGETS_GLOBAL_SETTINGS)
    target_settings = dict(global_settings)
    readset_annotation = annotate_readset_data_types(sample_desc)

    file_targets = []

    keep_path = False
    selected_path = None
    if 'select_target_path' in config:
        keep_path = True
        selected_path = config['select_target_path']
        if not selected_path in TARGET_PATHS:
            raise ValueError('Selected target path does not exist: '
                             '{} // {}'.format(selected_path, sorted(TARGET_PATHS.keys())))

    for target_specification in sample_targets:
        if 'aliases' in target_specification:
            continue
        elif 'defaults' in target_specification:
            target_spec = target_specification['defaults']
            target_settings = dict(global_settings)
            target_settings.update(target_spec)
        else:
            # Copy dict with default values for each target spec.
            # Note that target specs are processed in order,
            # so switching to other defaults for a second set of
            # target specs is possible
            target_spec = target_specification['target']
            tmp = dict(target_settings)
            tmp.update(target_spec)

            if not check_target_modifier_match(tmp,
                    CONFIG_TARGETS_SELECTED_KEYS,
                    CONFIG_TARGETS_SELECTED_VALUES,
                    True):
                if bool(config.get('show_warnings', False)):
                    sys.stderr.write('\nWARNING: discarding target spec: {}\n'.format(target_spec))
                continue

            if check_target_modifier_match(tmp,
                CONFIG_TARGETS_SKIPPED_KEYS,
                CONFIG_TARGETS_SKIPPED_VALUES,
                False):
                if bool(config.get('show_warnings', False)):
                    sys.stderr.write('\nWARNING: skipping over target spec: {}\n'.format(target_spec))
                continue

            for target_name, target_path in TARGET_PATHS.items():
                if keep_path and selected_path != target_name:
                    if target_name.startswith('INIT'):
                        # Init targets are always needed
                        pass
                    else:
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
                except (KeyError, WildcardError) as error:
                    if bool(config.get('show_warnings', False)):
                        sys.stderr.write('\nMissing parameter values for target {}: '
                                         '(known: {}) - {} [Skipping]\n'.format(target_name, tmp, str(error)))
                    continue
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
                    if file_target.startswith('input/'):
                        # skip over INIT targets
                        continue
                    full_path = os.path.join(root_dir, file_target)
                    assert os.path.isfile(full_path), \
                        'Non-existing file as build target: {} / {}'.format(full_path, output[0])
                    _ = dump.write(full_path + '\n')
