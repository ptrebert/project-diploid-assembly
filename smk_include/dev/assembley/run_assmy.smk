include: 'assmy_collect_input.smk'
include: 'assmy_preprocess.smk'
include: 'assmy_separate_tigs.smk'
include: 'assmy_connect_tigs.smk'
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
        hybrid_ont_align = expand(
            'output/hybrid/ont_to_graph/{sample_long}.{ont_type}.{tigs}.gaf',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL'],
            tigs=['TIGRAW']
        )
        # graph_ontul_gono_ref = expand(
        #     'output/read_aln/{sample_long}.{reference}.AMXYUN.tigs.{ont_type}.ga.gaf',
        #     sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
        #     reference=['T2Tv11_38p13Y_chm13'],
        #     ont_type=['ONTUL']
        # )


rule run_extract_confirmed_reads:
    input:
        xypar_reads = expand(
            'output/references/{sample_long}.XYPAR.reads.fastq.gz',
            sample_long=[
                'AFR-YRI-Y117-M_NA19239',
                'AMR-PUR-PR05-M_HG00731',
                'EAS-CHS-SH032-M_HG00512',
                'EUR-ASK-3140-M_NA24385'   
            ]
        )


rule run_mmap_to_aug_reference:
    input:
        mmap_align = expand(
            'output/read_aln/{sample_long}.{reference}.augY.{ont_type}.mmap.paf',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ONTUL_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            reference=['T2Tv11_38p13Y_chm13'],
            ont_type=['ONTUL']
        )