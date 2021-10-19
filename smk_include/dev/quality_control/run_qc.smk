
include: 'qc_aux.smk'
include: 'qc_collect_input.smk'
include: 'qc_preprocess.smk'
include: 'qc_ont_ec.smk'
include: 'qc_kmers.smk'
include: 'qc_qv_estimate.smk'
include: 'qc_read_cov.smk'

localrules: run_read_cov, run_qv_estimate

wildcard_constraints:
    sample = '(' + '|'.join([s for s in ONTUL_SAMPLES]) + ')'


def build_read_cov_targets(wildcards):

    template = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.cov.cache.h5'

    long_read_types = ['ONTUL', 'HIFIEC', 'HIFIAF']
    long_read_sets = ['guppy-5.0.11-sup-prom', 'hifiasm-v0.15.4', 'pgas-v14-dev']

    targets = []
    for sample in ONTUL_SAMPLES:
        for read_type, readset in zip(long_read_types, long_read_sets):
            formatter = {
                'sample': sample,
                'read_type': read_type,
                'readset': readset,
                'reference': config['reference']
            }
            fmt_target = template.format(**formatter)
            targets.append(fmt_target)
    return sorted(targets)


rule run_read_cov:
    input:
        build_read_cov_targets


def build_seqence_qv_estimate_targets(wildcards):

    template_seq = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_SHORT_{short_reads}_nohpc.seq-qv.h5'
    template_global = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_SHORT_{short_reads}_nohpc.qv.tsv'

    long_read_types = ['ONTUL', 'HIFIEC', 'HIFIAF']
    long_read_sets = ['guppy-5.0.11-sup-prom', 'hifiasm-v0.15.4', 'pgas-v14-dev']

    targets = []
    for sample in ONTUL_SAMPLES:
        short_readset = SAMPLE_INFOS[sample]['SHORT_RS']
        for long_reads, read_set in zip(long_read_types, long_read_sets):
            for hpc in ['ishpc', 'nohpc']:
                formatter = {
                    'sample': sample,
                    'long_reads': long_reads,
                    'readset': read_set,
                    'short_reads': short_readset,
                    'hpc': hpc
                }
                fmt_target = template_seq.format(**formatter)
                targets.append(fmt_target)
                fmt_target = template_global.format(**formatter)
                targets.append(fmt_target)

    template_seq = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_HIFIEC_hifiasm-v0.15.4_nohpc.seq-qv.h5'
    template_global = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_HIFIEC_hifiasm-v0.15.4_nohpc.qv.tsv'

    for sample in ONTUL_SAMPLES:
        for long_reads, read_set in zip(long_read_types, long_read_sets):
            for hpc in ['ishpc', 'nohpc']:
                if long_reads == 'HIFIEC':
                    continue
                formatter = {
                    'sample': sample,
                    'long_reads': long_reads,
                    'readset': read_set,
                    'hpc': hpc
                }
                fmt_target = template_seq.format(**formatter)
                targets.append(fmt_target)
                fmt_target = template_global.format(**formatter)
                targets.append(fmt_target)

    return sorted(targets)


rule run_qv_estimate:
    input:
        build_seqence_qv_estimate_targets


rule run_ont_error_correction:
    input:
        expand(
            'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.foobar',
            sample=ONTUL_SAMPLES,
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIAF'],
            graph_readset=['pgas-v14-dev']
        ),
        expand(
            'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.foobar',
            sample=ONTUL_SAMPLES,
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIEC'],
            graph_readset=['hifiasm-v0.15.4']
        ),
