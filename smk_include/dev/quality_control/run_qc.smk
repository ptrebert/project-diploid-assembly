
include: 'qc_aux.smk'
include: 'qc_collect_input.smk'
include: 'qc_preprocess.smk'
include: 'qc_kmers.smk'
include: 'qc_read_cov.smk'

localrules: run_read_cov, run_seq_stats

wildcard_constraints:
    sample = '(' + '|'.join([s for s in ONTUL_SAMPLES]) + ')'

rule run_read_cov:
    input:
        bed_cov_hifiec = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.aln.bed',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['HIFIEC'],
            readset=['hifiasm-v0.15.4'],
            reference=config['reference']
        ),
        bed_cov_ontul = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.aln.bed',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTUL'],
            readset=['guppy-5.0.11-sup-prom'],
            reference=config['reference']
        ),
        bed_cov_ontec_k501 = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.aln.bed',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIEC'],
            graph_readset=['hifiasm-v0.15.4'],
            kmer=[501],
            window=[100]
        ),
        bed_cov_ontec_k1001 = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.aln.bed',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIEC'],
            graph_readset=['hifiasm-v0.15.4'],
            kmer=[1001],
            window=[200]
        )


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
        stats_ontec_k501 = expand(
            'input/{read_type}/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.stats.tsv.gz',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIEC'],
            graph_readset=['hifiasm-v0.15.4'],
            kmer=[501],
            window=[100]
        ),
        stats_ontec_k1001 = expand(
            'input/{read_type}/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.stats.tsv.gz',
            sample=[s for s in ONTUL_SAMPLES],
            read_type=['ONTEC'],
            readset=['guppy-5.0.11-sup-prom'],
            graph_reads=['HIFIEC'],
            graph_readset=['hifiasm-v0.15.4'],
            kmer=[1001],
            window=[200]
        ),
