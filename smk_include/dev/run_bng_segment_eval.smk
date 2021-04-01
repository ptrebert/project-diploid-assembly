localrules: master

samples = [
    'HG00512',
    'HG00514',
    'HG00731',
    'HG00732',
    'HG00733',
    'NA19238',
    'NA19239',
    'NA19240'
]

output_files = []

for hap in ['H1', 'H2']:
    for s in samples:
        tmp = 'output/segment_align/{sample}_{hap}.wmap-k19.bam'.format(**{'sample': s, 'hap': hap})
        output_files.append(tmp)

rule extract_hap1_contigs:
    input:
        tig_names = '/beeond/data/hifiasm_v13_assemblies/tig_names/{sample}_H1_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta',
    output:
        assm_tigs = 'output/tigs/{sample}_H1_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'


rule extract_hap2_contigs:
    input:
        tig_names = '/beeond/data/hifiasm_v13_assemblies/tig_names/{sample}_H2_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta',
    output:
        assm_tigs = 'output/tigs/{sample}_H2_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'


rule count_hap1_assembly_kmers:
    input:
        fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
    output:
        kmer_db = directory('output/kmer/{sample}_H1.k19.db/'),
        rep_kmer = 'output/kmer/{sample}_H1.k19.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer/{sample}_H1.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 * attempt / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt
    shell:
        'meryl count k=19 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


rule count_hap2_assembly_kmers:
    input:
        fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'
    output:
        kmer_db = directory('output/kmer/{sample}_H2.k19.db/'),
        rep_kmer = 'output/kmer/{sample}_H2.k19.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer/{sample}_H2.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 * attempt / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt
    shell:
        'meryl count k=19 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


rule align_segments_to_assemblies:
    input:
        tigs = 'output/tigs/{sample}_{hap}_tigs.fasta',
        segments = '/beeond/data/hifiasm_v13_assemblies/chm13_H0_segments.fasta',
        kmers = 'output/kmer/{sample}_{hap}.k19.rep-grt09998.txt'
    output:
        'output/segment_align/{sample}_{hap}.wmap-k19.bam'
    benchmark:
        'rsrc/output/segment_align/{sample}_{hap}.CHM13.wmap.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 * attempt / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
    shell:
        'winnowmap -W {input.kmers} -t {threads} -ax asm20 {input.tigs} {input.segments} | samtools sort | samtools view -b > {output}'


rule master:
    input:
        output_files
