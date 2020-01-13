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

    "STATS_SAMPLE_SUMMARY": os.path.join(
        "output", "statistics", "stat_dumps",
        "{hap_reads}.{datatype}.pck"
    ),

    "STATS_VARIANT_CALLING": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sts_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.GQ{filter_vcf_gq}.vcf.stats"
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
        "{hap_reads}.tags.{tag_type}.tsv"
    ),
}


def extract_wildcard_config_parameters(sample_desc, wildcard_re):
    """
    :param sample_desc:
    :param wildcard_re:
    :return:
    """

    target_values = {}
    config_keys = set(config.keys())

    for target_name, target_path in TARGET_PATHS.items():
        target_wildcards = [m.strip('{}') for m in re.findall(wildcard_re, target_path)]
        if 'datatype' in target_wildcards or 'tag_type' in target_wildcards:
            if 'pbn' in sample_desc['datatype']:
                target_values['datatype'] = 'pbn.bam'
                target_values['tag_type'] = 'pbn'
            elif 'fastq' in sample_desc['datatype']:
                target_values['datatype'] = 'pbn.bam'
                target_values['tag_type'] = 'pbn'
            else:
                raise ValueError('Cannot handle data type of sample {}: {}'.format(sample_desc['individual'],
                                                                                   sample_desc['datatype']))
        shared_keys = set(target_wildcards).intersection(config_keys)
        for sk in shared_keys:
            target_values[sk] = config[sk]
    return target_values


def define_file_targets(wildcards):
    """
    :param wildcards:
    :return:
    """
    wildcard_re = '{\w+}'
    individual = wildcards.individual
    try:
        sample_desc = config['sample_description_' + individual]
        if individual != sample_desc['individual']:
            sys.stderr.write('\nIndividual mismatch: {} vs {}\n'.format(individual, sample_desc['individual']))
            raise ValueError('Sample description individual does not '
                             'match requested individual: {} vs {}'.format(individual, sample_desc['individual']))
    except KeyError as ke:
        sys.stderr.write('\nNo sample description for individual [sample_description_] {} in config\n'.format(individual))
        raise ke

    try:
        target_spec = config['sample_targets_' + individual]
    except KeyError as ke:
        sys.stderr.write('\nNo targets specified for individual [target_specification_] {} in config\n'.format(individual))
        raise ke

    specs = extract_wildcard_config_parameters(sample_desc, wildcard_re)
    file_targets = []

    for spec_type, spec_params in target_spec:
        if spec_type == 'aliases':
            continue
        elif spec_type == 'defaults':
            specs.update(spec_params)
        else:
            # Copy dict with default values for each target spec.
            # Note that target specs are processed in order,
            # so switching to other defaults for a second set of
            # target specs is possible
            tmp = dict(specs)
            tmp.update(spec_params)
            for target_name, target_path in TARGET_PATHS.items():
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
    return sorted(file_targets)


rule dump_build_targets:
    input:
         targets = define_file_targets
    output:
        'output/targets/{population}_{family}/{individual}.fofn'
    run:
        if not input.targets:
            try:
                os.unlink(output[0])
            except (OSError, IOError):
                pass
        else:
            with open(output[0], 'w') as dump:
                for file_target in input.targets:
                    _ = dump.write(file_target + '\n')
