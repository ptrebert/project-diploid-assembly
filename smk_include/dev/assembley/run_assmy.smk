include: 'assmy_collect_input.smk'
include: 'assmy_separate_tigs.smk'
include: 'assmy_connect_tigs.smk'

localrules: run_all

rule run_all:
    input:
        tig_aln = expand(
            'output/tig_aln/{sample_long}_MAP-TO_{reference}.psort.bam',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ASSEMBLED_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13']
        ),
        tigs_to_chrom = expand(
            'output/tig_aln/chrom_groups/{sample_long}_MAP-TO_{reference}.tigs{ext}',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ASSEMBLED_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            ext=[
                '.AM.pass.txt',
                '.AM.fail.txt',
                '.XY.pass.txt',
                '.XY.fail.txt',
                '.AMXY.fail.txt',
                '.UN.fail.txt'
            ]
        ),
        gono_ref = expand(
            'output/gonosomal_reference/{sample_long}.{reference}.AMXYUN.tigs.stats.tsv',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ASSEMBLED_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
        ),
        lin_ontul_gono_ref = expand(
            'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            kmer_size=[15],
            hpc=['nohpc'],
            ont_type=['ONTUL']
        ),
        lin_ontec_gono_ref = expand(
            'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTEC_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            kmer_size=[15],
            hpc=['nohpc'],
            ont_type=['ONTEC']
        ),
        graph_ontul_gono_ref = expand(
            'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.{ont_type}.ga.gaf',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL']
        ),