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
        ),
        hybrid_contig_ref_align = expand(
            'output/hybrid/210_align_ref/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.h5',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL'],
            tigs=['TIGRAW'],
            reference=['T2Tv11_hg002Yv2_chm13']
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


def make_combinations(samples, sample_long_desc, read_types, mapq_thresholds):

    wc_name_samplelong = set(t[0] for t in sample_long_desc).pop()
    for wc_name_sample, sample in samples:
        for wc_name_readtype, read_type in read_types:
            for wc_name_mapqt, mapq_trhreshold in mapq_thresholds:
                formatter = {
                    wc_name_sample: sample,
                    wc_name_samplelong: SAMPLE_INFOS[sample]['long_id'],
                    wc_name_readtype: read_type,
                    wc_name_mapqt: mapq_trhreshold
                }
                yield formatter
    return


rule run_chry_targeted_t2t:
    input:
        target_chry = expand(
            'output/target_assembly/chry_reads/{sample}/{sample_long}_{read_type}.{mapq}.p_ctg.gfa',
            make_combinations,
            sample=[sample for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            read_type=['HIFIAF', 'HIFIEC'],
            mapq=['mq00', 'mq60'],
        ),
        aln_cache = expand(
            'output/hybrid/220_gfa_annotation/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL'],
            tigs=['ECMQ0YRAW', 'AFMQ0YRAW'],
            reference=['T2Tv11_hg002Yv2_chm13']
        )


rule run_chry_targeted_lja:
    input:
        target_chry = expand(
            'output/target_assembly/chry_reads/lja/{sample_long}_{read_type}.{mapq}/assembly.fasta',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ['HG02666', 'HC02666', 'NA18989', 'HC18989']],
            read_type=['HIFIAF', 'HIFIEC'],
            mapq=['mq00']
        )


rule run_afr_mix_tests:
    input:
        mbg_duo = expand(
            'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.gfa',
            zip,
            sample=['NA193N7'] * 3,
            sample_info=['AFR-LWK-DUO-M'] * 3,
            read_type=['HIFIEC'] * 3,
            mapq=['mq00'] * 3,
            kmer=[501, 1001, 2001],
            window=[100, 200, 400]
        ),
        mbg_trio = expand(
            'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.gfa',
            zip,
            sample=['NA193NN'] * 3,
            sample_info=['AFR-LWK-TRIO-M'] * 3,
            read_type=['HIFIEC'] * 3,
            mapq=['mq00'] * 3,
            kmer=[501, 1001, 2001],
            window=[100, 200, 400]
        ),
        mbg_quart = expand(
            'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.gfa',
            zip,
            sample=['AFR4MIX'] * 3,
            sample_info=['AFR-NNN-QUART-M'] * 3,
            read_type=['HIFIEC'] * 3,
            mapq=['mq00'] * 3,
            kmer=[501, 1001, 2001],
            window=[100, 200, 400]
        ),
        hybrid_assm = expand(
            'output/hybrid/210_align_ref/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.h5',
            sample_long=['AFR-LWK-DUO-M_NA193N7', 'AFR-LWK-TRIO-M_NA193NN', 'AFR-NNN-QUART-M_AFR4MIX'],
            ont_type=['ONTUL'],
            tigs=['ECMQ0YRAW', 'AFMQ0YRAW'],
            reference=['T2Tv11_hg002Yv2_chm13']
        )
