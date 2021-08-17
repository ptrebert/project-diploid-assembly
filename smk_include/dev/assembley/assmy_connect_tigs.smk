
rule determine_frequent_kmers_gono_reference:
    input:
        fasta = 'output/gonosomal_reference/fasta/{sample_long}.{reference}.AMXYUN.tigs.fasta'
    output:
        db = directory('output/gonosomal_reference/kmer_db/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl'),
        repkmer = 'output/gonosomal_reference/kmer_db/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl.repkmer-grt09998.txt'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_gb = lambda wildcards, attempt: 3 * attempt,
    params:
        hpc = lambda wildcards: 'compress' if wildcards.hpc == 'ishpc' else '',
    shell:
        'meryl count k={wildcards.kmer_size} threads={threads} memory={resources.mem_total_gb} {params.hpc} output {output.db} {input.fasta}'
        ' && '
        'meryl print greater-than distinct=0.9998 {output.db} > {output.rep_kmer}'


rule wmap_align_ont_to_gono_reference:
    input:
        fasta = 'output/gonosomal_reference/fasta/{sample}_{sample_info}.{reference}.AMXYUN.tigs.fasta',
        rep_kmer = 'output/gonosomal_reference/kmer_db/{sample}_{sample_info}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl.repkmer-grt09998.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        bam = 'output/read_aln/{sample}_{sample_info}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam',
        bai = 'output/read_aln/{sample}_{sample_info}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam.bai',
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
    params:
        preset = lambda wildcards: 'map-pb' if wildcards.ont_type == 'ONTEC' else 'map-ont',
        sort_threads = 4,
        sort_mem = 2048
    shell:
        'winnowmap -k {wildcards.kmer_size} -W {input.rep_kmer} -x {params.preset} --MD -Y --eqx -L -a --secondary=no '
            '-R "@RG\\tID:{wildcards.sample}_{wildcards.sample_info}_{wildcards.ont_type}\\tSM:{wildcards.sample}" '
            '{input.fasta} {input.reads}'
        ' | '
        'samtools view -F 260 -u -q 20'
        ' | '
        'samtools sort -@ {params.sort_threads} -m {params.sort_mem}M -O BAM --no-PG -l 6 > {output.bam}'
        ' && '
        'samtools index {output.bam}'


rule ga_align_ont_to_gono_reference:
    input:
        container = ancient('graphaligner.sif'),
        graph = 'output/gonosomal_reference/graph/{sample}_{sample_info}.{reference}.AMXYUN.tigs.gfa',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        gaf = 'output/read_aln/{sample}_{sample_info}.{reference}.AMXYUN.tigs.{ont_type}.ga.gaf',
        hybrid_reads = 'output/read_aln/{sample}_{sample_info}.{reference}.AMXYUN.{ont_type}.ONTHY.fasta',
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
    params:
        preset = 'vg'
    shell:
        'module load Singularity && singularity exec {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 1 '
            '--corrected-out {output.hybrid_reads} '
            '-a {output.gaf} &> {log}'
