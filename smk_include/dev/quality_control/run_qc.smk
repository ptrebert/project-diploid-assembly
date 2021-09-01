
include: 'qc_aux.smk'
include: 'qc_collect_input.smk'
include: 'qc_preprocess.smk'
include: 'qc_kmers.smk'
include: 'qc_read_cov.smk'

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
