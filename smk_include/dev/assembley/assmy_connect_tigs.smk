
rule determine_frequent_kmers_gono_reference:
    input:
        fasta = 'output/gonosomal_reference/fasta/{sample_long}.{reference}.AMXYUN.tigs.fasta'
    output:
        db = directory('output/gonosomal_reference/kmer_db/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl'),
        rep_kmer = 'output/gonosomal_reference/kmer_db/{sample_long}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl.repkmer-grt09998.txt'
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
        fasta = 'output/gonosomal_reference/fasta/{sample_info}_{sample}.{reference}.AMXYUN.tigs.fasta',
        rep_kmer = 'output/gonosomal_reference/kmer_db/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.meryl.repkmer-grt09998.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        bam = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam',
        bai = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam.bai',
    log:
        'log/output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.log'
    benchmark:
        'rsrc/output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.rsrc'
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 * attempt,
        runtime_hrs = lambda wildcards, attempt: 167,
    params:
        preset = lambda wildcards: 'map-pb' if wildcards.ont_type == 'ONTEC' else 'map-ont',
        sort_threads = 4,
        sort_mem = 4096
    shell:
        'winnowmap -k {wildcards.kmer_size} -W {input.rep_kmer} -x {params.preset} --MD -Y --eqx -L -a --secondary=no '
            '-R "@RG\\tID:{wildcards.sample}_{wildcards.sample_info}_{wildcards.ont_type}\\tSM:{wildcards.sample}" '
            '{input.fasta} {input.reads} 2> {log}'
        ' | '
        'samtools view -F 260 -u -q 20'
        ' | '
        'samtools sort -@ {params.sort_threads} -m {params.sort_mem}M -O BAM --no-PG -l 6 > {output.bam}'
        ' && '
        'samtools index {output.bam}'


rule dump_ont_to_gono_reference:
    input:
        bam = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam',
        bai = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.psort.bam.bai',
    output:
        bed = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.k{kmer_size}.{hpc}.{ont_type}.wmap.cov.bed',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 1
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'bedtools bamtobed -i {input.bam} > {output.bed}'


rule ga_align_ont_to_gono_reference:
    input:
        container = ancient('graphaligner.sif'),
        graph = 'output/gonosomal_reference/graph/{sample_info}_{sample}.{reference}.AMXYUN.tigs.gfa',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        gaf = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.{ont_type}.ga.gaf',
        hybrid_reads = 'output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.{ont_type}.ONTHY.fasta',
    log:
        'log/output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.{ont_type}.ga.log'
    benchmark:
        'rsrc/output/read_aln/{sample_info}_{sample}.{reference}.AMXYUN.tigs.{ont_type}.ga.rsrc'
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES
#    conda: '../../../environment/conda/conda_biotools.yml'
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


rule mmap_align_ont_to_aug_reference:
    input:
        fasta = 'output/references/{reference}.augY.fasta',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        paf = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.paf',
    log:
        'log/output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.log',
    benchmark:
        'rsrc/output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.rsrc'
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 167,
    params:
        preset = lambda wildcards: 'map-pb' if wildcards.ont_type == 'ONTEC' else 'map-ont',
    shell:
        'minimap2 -x {params.preset} --secondary=no -o {output.paf} {input.fasta} {input.reads} &> {log}'
        