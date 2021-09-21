include: 'assmy_collect_input.smk'
include: 'assmy_preprocess.smk'
include: 'assmy_separate_tigs.smk'
include: 'assmy_connect_tigs.smk'
include: 'assmy_target_xypar.smk'
include: 'assmy_finnish_hybrid.smk'

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
            'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.cov.bed',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            kmer_size=[15],
            hpc=['nohpc'],
            ont_type=['ONTUL']
        ),
        lin_ontec_gono_ref = expand(
            'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.cov.bed',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTEC_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            kmer_size=[15],
            hpc=['nohpc'],
            ont_type=['ONTEC']
        ),
        # graph_ontul_gono_ref = expand(
        #     'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.{ont_type}.ga.gaf',
        #     sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
        #     reference=['T2Tv11_38p13Y_chm13'],
        #     ont_type=['ONTUL']
        # )


rule run_hybrid_assembly:
    input:
        hybrid_ont_align = expand(
            'output/hybrid/110_final_graph/{sample_long}.{ont_type}.{tigs}.final.gfa',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL'],
            tigs=['TIGRAW']
        )


# rule run_extract_confirmed_reads:
#     input:
#         xypar_reads = expand(
#             'output/references/{sample_long}.XYPAR.reads.fastq.gz',
#             sample_long=[
#                 'AFR-YRI-Y117-M_NA19239',
#                 'AMR-PUR-PR05-M_HG00731',
#                 'EAS-CHS-SH032-M_HG00512',
#                 'EUR-ASK-3140-M_NA24385'   
#             ]
#         ),
#         xypar_assm = 'output/target_assembly/xypar_reads/xypar.p_utg.gfa'


# rule run_mmap_to_aug_reference:
#     input:
#         xypar_reads = expand(
#             'output/read_subsets/xypar/{sample_long}.{reference}.augY.{ont_type}.{chrom}-reads.fasta.gz',
#             sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
#             reference=['T2Tv11_38p13Y_chm13'],
#             ont_type=['ONTUL'],
#             chrom=['chrX', 'chrY']
#         )


rule run_extract_aligned_chry_reads:
    input:
        fasta = expand(
            'output/read_subsets/chry/{sample_long}_{read_type}.chrY-reads.{mapq}.fasta.gz',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            read_type=['HIFIEC', 'ONTUL'],
            mapq=['mq00', 'mq60'],
        ),
        fastq = expand(
            'output/read_subsets/chry/{sample_long}_{read_type}.chrY-reads.{mapq}.fastq.gz',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            read_type=['HIFIAF'],
            mapq=['mq00', 'mq60'],
        ),
