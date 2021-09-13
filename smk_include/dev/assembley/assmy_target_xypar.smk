

rule hifiasm_xypar_targeted_assembly:
    input:
        fastq = expand(
            'output/references/{sample_long}.XYPAR.reads.fastq.gz',
            sample_long=[
                'AFR-YRI-Y117-M_NA19239',
                'AMR-PUR-PR05-M_HG00731',
                'EAS-CHS-SH032-M_HG00512',
                'EUR-ASK-3140-M_NA24385'   
            ]
        )
    output:
        primary_unitigs = 'output/target_assembly/xypar_reads/xypar.p_utg.gfa',
        primary_contigs = 'output/target_assembly/xypar_reads/xypar.p_ctg.gfa',
        raw_unitigs = 'output/target_assembly/xypar_reads/xypar.r_utg.gfa',
        discard = multiext(
                'output/target_assembly/xypar_reads/xypar',
                '.a_ctg.gfa', '.a_ctg.noseq.gfa',
                '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
                '.r_utg.noseq.gfa'
            ),
        ec_reads = 'output/target_assembly/xypar_reads/xypar.ec.fa',
        reads_ava = 'output/target_assembly/xypar_reads/xypar.ovlp.paf'
    log:
        hifiasm = 'log/output/target_assembly/xypar_reads.hifiasm.log',
    benchmark:
        os.path.join('rsrc/output/target_assembly/xypar_reads.hifiasm' + '.t{}.rsrc'.format(config['num_cpu_medium']))
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-ec --write-paf --primary {input.fastq} &> {log.hifiasm}'


# rule mbg_xypar_targeted_assembly:
#     input:
#         fastq = expand(
#             'output/references/{sample_long}.XYPAR.reads.fastq.gz',
#             sample_long=[
#                 'AFR-YRI-Y117-M_NA19239',
#                 'AMR-PUR-PR05-M_HG00731',
#                 'EAS-CHS-SH032-M_HG00512',
#                 'EUR-ASK-3140-M_NA24385'   
#             ]
#         )
#     output:
#         'output/target_assembly/xypar_kmers/'
#     log:

#     benchmark:

#     conda:
#         '../../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_medium']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 180224 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 12 * attempt
#     shell:
#         'MBG '