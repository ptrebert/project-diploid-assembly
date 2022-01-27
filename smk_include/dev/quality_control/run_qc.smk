
include: 'qc_aux.smk'
include: 'qc_collect_input.smk'
include: 'qc_preprocess.smk'
include: 'qc_ont_ec.smk'
include: 'qc_kmers.smk'
include: 'qc_qv_estimate.smk'
include: 'qc_read_cov.smk'

localrules: run_read_cov, run_qv_estimate, run_ont_error_correction

COMPLETE_SAMPLES = sorted(set(ONTUL_SAMPLES).intersection(HIFIEC_SAMPLES))

wildcard_constraints:
    sample = '(' + '|'.join(COMPLETE_SAMPLES) + ')'


def build_read_cov_targets(wildcards):

    template_cov_cache = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_T2TXYM.cov.cache.h5'
    template_summary = 'input/{read_type}/{sample}_{read_type}_{readset}.stats-summary.tsv'
    template_dump = 'input/{read_type}/{sample}_{read_type}_{readset}.stats-dump.pck'

    long_read_types = ['ONTUL', 'HIFIEC', 'HIFIAF', 'ONTEC']
    long_read_sets = [RS_ONTUL, RS_HIFIEC, RS_HIFIAF, RS_ONTEC]

    targets = []
    for sample in COMPLETE_SAMPLES:
        for read_type, readset in zip(long_read_types, long_read_sets):
            formatter = {
                'sample': sample,
                'read_type': read_type,
                'readset': readset,
            }
            fmt_target = template_cov_cache.format(**formatter)
            targets.append(fmt_target)

            fmt_target = template_dump.format(**formatter)
            targets.append(fmt_target)

            fmt_target = template_summary.format(**formatter)
            targets.append(fmt_target)
    return sorted(targets)


rule run_read_cov:
    input:
        build_read_cov_targets


def build_seqence_qv_estimate_targets(wildcards):

    template_seq = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_SHORT_{short_reads}_nohpc.seq-qv.h5'
    template_global = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_SHORT_{short_reads}_nohpc.qv.tsv'

    long_read_types = ['ONTUL', 'HIFIEC', 'HIFIAF', 'ONTEC']
    long_read_sets = [RS_ONTUL, RS_HIFIEC, RS_HIFIAF, RS_ONTEC]

    targets = []
    for sample in COMPLETE_SAMPLES:
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

    template_seq = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_HIFIEC_' + f'{RS_HIFIEC}_nohpc.seq-qv.h5'
    template_global = 'output/qv_estimate/{sample}_{long_reads}_{readset}_{hpc}_REF_HIFIEC_' + f'{RS_HIFIEC}_nohpc.qv.tsv'

    for sample in COMPLETE_SAMPLES:
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
            'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.stats.h5',
            sample=COMPLETE_SAMPLES,
            read_type=['ONTEC'],
            readset=[RS_ONTUL],
            graph_reads=['HIFIEC'],
            graph_readset=[RS_HIFIEC],
            kmer=config['mbg_init_kmer'],
            window=config['mbg_window_size'],
            resolve=config['mbg_resolve_kmer']
        ),
        expand(
            'input/ONTEC/{sample}_ONTEC_{readset}.fasta.gz',
            sample=COMPLETE_SAMPLES,
            readset=RS_ONTEC
        )
