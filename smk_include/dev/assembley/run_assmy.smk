include: 'assmy_collect_input.smk'
include: 'assmy_separate_tigs.smk'

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
        )