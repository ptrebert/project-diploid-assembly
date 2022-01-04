include: 'assmy_aux.smk'
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


COMPLETE_SAMPLES = sorted(set(ONTUL_SAMPLES).intersection(set(HIFIEC_SAMPLES)))
if 'NA19317' in COMPLETE_SAMPLES and 'NA19347' in COMPLETE_SAMPLES:
    COMPLETE_SAMPLES.append('NA193N7')


rule run_extract_read_subsets:
    input:
        fasta = expand(
            'output/read_subsets/{chrom}/{sample_long}_{read_type}.{chrom}-reads.{mapq}.fasta.gz',
            chrom=['chrX', 'chrY', 'chrXY'],
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in COMPLETE_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            read_type=['HIFIEC', 'ONTUL', 'ONTEC', 'OHEC'],
            mapq=['mq00'],
        ),
        fastq = expand(
            'output/read_subsets/{chrom}/{sample_long}_{read_type}.{chrom}-reads.{mapq}.fastq.gz',
            chrom=['chrX', 'chrY', 'chrXY'],
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in COMPLETE_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            read_type=['HIFIAF'],
            mapq=['mq00'],
        )


HIFIASM_HYBRID_TIGS = [
    'HAS-UTGRAW-HECMQ0Y',
    'HAS-UTGPRI-HECMQ0Y',
    'HAS-UTGRAW-HECMQ0XY',
    'HAS-UTGPRI-HECMQ0XY',
    'HAS-UTGRAW-OHECMQ0Y',
    'HAS-UTGPRI-OHECMQ0Y',
    'HAS-UTGRAW-OHECMQ0XY',
    'HAS-UTGPRI-OHECMQ0XY',
]

rule run_targeted_hifiasm_hybrid:
    input:
        gfa_labels_regular = expand(
            'output/hybrid/220_gfa_annotation/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in COMPLETE_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M'],
            ont_type=['ONTUL', 'UNSET'],
            tigs=HIFIASM_HYBRID_TIGS,
            reference=['T2Tv11_hg002Yv2_chm13']
        ),
        gfa_labels_duo = expand(
            'output/hybrid/220_gfa_annotation/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv',
            sample_long=['AFR-LWK-DUO-M_NA193N7'],
            ont_type=['ONTUL', 'UNSET'],
            tigs=HIFIASM_HYBRID_TIGS,
            reference=['T2Tv11_hg002Yv2_chm13']
        )


def define_mbg_hybrid_targets(wildcards):

    target_files = []

    params_info = [(phash, *pvalues) for phash, pvalues in MBG_PARAMS.items()]
    params_info = ['assembler_params/MBG_{}_k{}-w{}-r{}.info'.format(*p) for p in params_info]
    target_files.extend(params_info)

    tigs_hifiec = ['HECMQ0Y', 'HECMQ0XY']
    tigs_ohec = ['OHECMQ0Y', 'OHECMQ0XY']
    tigs_spec = tigs_hifiec + tigs_ohec

    samples_long = [SAMPLE_INFOS[sample]['long_id'] for sample in COMPLETE_SAMPLES if SAMPLE_INFOS[sample]['sex'] == 'M']
    if 'NA193N7' in COMPLETE_SAMPLES:
        samples_long.append('AFR-LWK-DUO-M_NA193N7')  # AFR duo mix

    ont_types = ['ONTUL', 'UNSET']  # unset = input assemblies, non-hybrid
    reference = 'T2Tv11_hg002Yv2_chm13'

    mbg_spec = [(int(get_mbg_param(spec_hash, 'resolve')), spec_hash) for spec_hash in MBG_PARAMS.keys()]

    template = 'output/hybrid/220_gfa_annotation/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv'

    for sample in samples_long:
        for ont_type in ont_types:
            for tigs in tigs_spec:
                for rk, spec_key in mbg_spec:
                    if rk > 20000 and tigs.startswith('HEC'):
                        # large resolveK for HIFIEC-only read sets
                        # is pointless
                        continue
                    formatter = {
                        'sample_long': sample,
                        'ont_type': ont_type,
                        'tigs': f'MBG-{spec_key}-{tigs}',
                        'reference': reference
                    }
                    out_file = template.format(**formatter)
                    target_files.append(out_file)
    return target_files


rule run_targeted_mbg_hybrid:
    input:
        define_mbg_hybrid_targets


rule run_chry_targeted_lja_input:
    input:
        gfa_labels = expand(
            'output/hybrid/220_gfa_annotation/{sample_long}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv',
            sample_long=[SAMPLE_INFOS[sample]['long_id'] for sample in ['HG02666', 'HC02666', 'NA18989', 'HC18989']],
            ont_type=['UNSET'],
            tigs=[f'LJA-{lja_params}-ECMQ0Y' for lja_params in LJA_PARAMS.keys()],
            reference=['T2Tv11_hg002Yv2_chm13']
        )
