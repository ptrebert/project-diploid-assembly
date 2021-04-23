localrules: master

ASSEMBLY_GRAPHS_FOLDER = '/gpfs/project/projects/medbioinf/projects/sig_chrY/run_folder/input/graphs'
REFERENCES_FOLDER = '/gpfs/project/projects/medbioinf/data/references'

rule dump_unitgs_to_fasta:
    input:
        gfa = os.path.join(ASSEMBLY_GRAPHS_FOLDER, '{assembly_graph}.{tigs}.gfa')
    output:
        fasta = 'output/sequences/{assembly_graph}.{tigs}.fasta',
        stats = 'output/stats/{assembly_graph}.{tigs}.contig.stats',
        mapping = 'output/stats/{assembly_graph}.{tigs}.read-contig.map'
    log:
        'log/output/sequences/{assembly_graph}.{tigs}.dump.log'
    benchmark:
        'rsrc/output/sequences/{assembly_graph}.{tigs}.dump.rsrc'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
        '--out-fasta {output.fasta} --out-map {output.mapping} '
        '--out-stats {output.stats} &> {log}'


rule align_tig_sequences_to_reference:
    input:
        reference = os.path.join(REFERENCES_FOLDER, '{reference}.fasta'),
        rep_kmer = os.path.join(REFERENCES_FOLDER, '{reference}.k{kmer}.rep-grt09998.txt'),
        tig_seq = 'output/sequences/{assembly_graph}.{tigs}.fasta',
    output:
        bam = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.wmap.bam'
    benchmark:
        'rsrc/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.wmap.rsrc'
    log:
        'log/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.wmap.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 36864 + 36864 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        sample = lambda wildcards: wildcards.assembly_graph.split('_')[0],
        sort_threads = 4,
        align_threads = int(config['num_cpu_high'] - 4)
    shell:
        'winnowmap -W {input.rep_kmer} -t {params.align_threads} -ax map-pb -R "@RG\\tID:1\\tSM:{params.sample}" {input.reference} {input.tig_seq} | '
        'samtools sort -m 4096M --threads {params.sort_threads} --output-fmt BAM -l 9 > {output.bam}'


rule dump_tig_alignments_to_bed:
    input:
        bam = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.wmap.bam'
    output:
        bed = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.bed.gz'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        compress_threads = int(config['num_cpu_low'] - 1)
    shell:
        'bedtools bamtobed -splitD -i {input.bam} | pigz -p {params.compress_threads} > {output.bed}'


OUTPUT_FILES = [
    'output/alignments/tig_to_ref/HG00512_hgsvc_pbsq2-ccs_1000.r_utg_MAP-TO_T2Tv1_38p13Y_chm13.k15.bed.gz',
    'output/alignments/tig_to_ref/HG00731_hgsvc_pbsq2-ccs_1000.r_utg_MAP-TO_T2Tv1_38p13Y_chm13.k15.bed.gz',
    'output/alignments/tig_to_ref/NA19239_hgsvc_pbsq2-ccs_1000.r_utg_MAP-TO_T2Tv1_38p13Y_chm13.k15.bed.gz',
    'output/alignments/tig_to_ref/NA24385_hpg_pbsq2-ccs_1000.r_utg_MAP-TO_T2Tv1_38p13Y_chm13.k15.bed.gz',
]

rule master:
    input:
        OUTPUT_FILES