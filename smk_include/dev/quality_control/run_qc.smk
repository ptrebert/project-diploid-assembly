
include: 'qc_aux.smk'
include: 'qc_collect_input.smk'
include: 'qc_preprocess.smk'
include: 'qc_ont_ec.smk'
include: 'qc_kmers.smk'
include: 'qc_qv_estimate.smk'
include: 'qc_read_cov.smk'

localrules: run_read_cov, run_seq_stats

wildcard_constraints:
    sample = '(' + '|'.join([s for s in ONTUL_SAMPLES]) + ')'

rule run_read_cov:
    input:
        read_cov_hifiec = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.paf.gz',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['HIFIEC'],
            readset=['hifiasm-v0.15.4'],
            reference=config['reference']
        ),
        read_cov_ontul = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.paf.gz',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTUL'],
            readset=['guppy-5.0.11-sup-prom'],
            reference=config['reference']
        ),
        # bed_cov_ontec_k501 = expand(
        #     'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.aln.bed',
        #     sample=[s for s in ONTUL_SAMPLES],
        #     read_type=['ONTEC'],
        #     readset=['guppy-5.0.11-sup-prom'],
        #     graph_reads=['HIFIEC'],
        #     graph_readset=['hifiasm-v0.15.4'],
        #     kmer=[501],
        #     window=[100]
        # ),
        # bed_cov_ontec_k1001 = expand(
        #     'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.aln.bed',
        #     sample=[s for s in ONTUL_SAMPLES],
        #     read_type=['ONTEC'],
        #     readset=['guppy-5.0.11-sup-prom'],
        #     graph_reads=['HIFIEC'],
        #     graph_readset=['hifiasm-v0.15.4'],
        #     kmer=[1001],
        #     window=[200]
        # )


rule run_seq_stats:
    input:
        stats_hifiec = expand(
            'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz',
            read_type=['HIFIEC'],
            sample=[s for s in ONTUL_SAMPLES],
            readset=['hifiasm-v0.15.4']
        ),
        stats_ontul = expand(
            'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz',
            read_type=['ONTUL'],
            sample=[s for s in ONTUL_SAMPLES],
            readset=['guppy-5.0.11-sup-prom']
        ),
        # stats_ontec_k501 = expand(
        #     'input/{read_type}/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.stats.tsv.gz',
        #     sample=[s for s in ONTUL_SAMPLES],
        #     read_type=['ONTEC'],
        #     readset=['guppy-5.0.11-sup-prom'],
        #     graph_reads=['HIFIEC'],
        #     graph_readset=['hifiasm-v0.15.4'],
        #     kmer=[501],
        #     window=[100]
        # ),
        # stats_ontec_k1001 = expand(
        #     'input/{read_type}/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.stats.tsv.gz',
        #     sample=[s for s in ONTUL_SAMPLES],
        #     read_type=['ONTEC'],
        #     readset=['guppy-5.0.11-sup-prom'],
        #     graph_reads=['HIFIEC'],
        #     graph_readset=['hifiasm-v0.15.4'],
        #     kmer=[1001],
        #     window=[200]
        # ),


def build_seqence_qv_estimate_targets(wildcards):

    template_seq = 'output/qv_estimate/{sample}_{long_reads}_{readset}_REF_SHORT_{short_reads}.seq-qv.h5'
    template_global = 'output/qv_estimate/{sample}_{long_reads}_{readset}_REF_SHORT_{short_reads}.qv.tsv'

    long_read_types = ['ONTUL', 'HIFIEC', 'HIFIAF']
    long_read_sets = ['guppy-5.0.11-sup-prom', 'hifiasm-v0.15.4', 'pgas-v14-dev']

    targets = []
    for sample in SHORT_SAMPLES:
        short_readset = SAMPLE_INFOS[sample]['SHORT_RS']
        for long_reads, read_set in zip(long_read_types, long_read_sets):
            formatter = {
                'sample': sample,
                'long_reads': long_reads,
                'readset': read_set,
                'short_reads': short_readset
            }
            fmt_target = template_seq.format(**formatter)
            targets.append(fmt_target)
            fmt_target = template_global.format(**formatter)
            targets.append(fmt_target)
    return sorted(targets)


rule run_qv_estimate:
    input:
        build_seqence_qv_estimate_targets