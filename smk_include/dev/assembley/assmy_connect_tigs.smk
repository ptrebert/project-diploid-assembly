
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
        sample = CONSTRAINT_REGULAR_SAMPLES
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
        sample = CONSTRAINT_REGULAR_SAMPLES
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
        sample = CONSTRAINT_REGULAR_SAMPLES
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 167,
    params:
        preset = lambda wildcards: 'map-pb' if wildcards.ont_type == 'ONTEC' else 'map-ont',
    shell:
        'minimap2 -x {params.preset} --secondary=no -o {output.paf} {input.fasta} {input.reads} &> {log}'


PAF_HEADER = [
    'query_name',
    'query_length',
    'query_start',
    'query_end',
    'orientation',
    'target_name',
    'target_length',
    'target_start',
    'target_end',
    'res_matches',
    'block_length',
    'mapq',
    'aln_type',
    'tag_cm',
    'tag_s1',
    'tag_s2',
    'divergence',
    'tag_rl'
]


rule collect_xy_read_statistics:
    input:
        paf = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.paf',
    output:
        stats = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.stats.tsv',
        chrx_reads = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.chrX-reads.txt',
        chry_reads = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.chrY-reads.txt',
    wildcard_constraints:
        sample = CONSTRAINT_REGULAR_SAMPLES
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    run:
        import pandas as pd
        import numpy as np
    
        autosomes = [f'chr{i}' for i in range(1,23)]

        df = pd.read_csv(
            input.paf,
            sep='\t',
            header=None,
            names=PAF_HEADER,
            usecols=['query_name', 'query_length', 'target_name', 'mapq', 'divergence']
        )
        df['divergence'] = df['divergence'].apply(lambda x: float(x.split(':')[-1]))
        df['x_or_y'] = 1
        df.loc[df['target_name'].isin(autosomes), 'x_or_y'] = -1
        df.loc[df['target_name'].isin(['chrX']), 'x_or_y'] = 0

        is_chry = set()
        chry_lengths = []
        chry_divergences = []
        is_chrx = set()
        chrx_lengths = []
        chrx_divergences = []
        ambig = 0

        for read, alignments in df.groupby('query_name'):
            if (alignments['x_or_y'] == -1).all():
                continue
            if (alignments['x_or_y'] == -1).any():
                ambig += 1
                continue
            if (alignments['x_or_y'] == 1).all():
                is_chry.add(read)
                chry_lengths.append(alignments['query_length'].values[0])
                chry_divergences.extend(alignments['divergence'].tolist())
            else:
                is_chrx.add(read)
                chrx_lengths.append(alignments['query_length'].values[0])
                chrx_divergences.extend(alignments['divergence'].tolist())

        chry_divergences = np.array(chry_divergences, dtype=np.float16)
        chrx_divergences = np.array(chrx_divergences, dtype=np.float16)

        with open(output.chrx_reads, 'w') as dump:
            _ = dump.write('\n'.join(sorted(is_chrx)) + '\n')
        with open(output.chry_reads, 'w') as dump:
            _ = dump.write('\n'.join(sorted(is_chry)) + '\n')
        
        with open(output.stats, 'w') as dump:
            _ = dump.write(f'total_alignments\t{df.shape[0]}\n')
            _ = dump.write(f'total_reads\t{df["query_name"].nunique()}\n')
            _ = dump.write(f'ambiguous_reads\t{ambig}\n')

            _ = dump.write(f'chrY_num_reads\t{len(is_chry)}\n')
            _ = dump.write(f'chrY_sum_length\t{sum(chry_lengths)}\n')
            _ = dump.write(f'chrY_mean_divergence\t{chry_divergences.mean()}\n')
            _ = dump.write(f'chrY_stddev_divergence\t{chry_divergences.std()}\n')

            _ = dump.write(f'chrX_num_reads\t{len(is_chrx)}\n')
            _ = dump.write(f'chrX_sum_length\t{sum(chrx_lengths)}\n')
            _ = dump.write(f'chrX_mean_divergence\t{chrx_divergences.mean()}\n')
            _ = dump.write(f'chrX_stddev_divergence\t{chrx_divergences.std()}\n')
    # END OF RUN BLOCK


rule extract_xy_reads:
    input:
        read_names = 'output/read_aln/{sample_info}_{sample}.{reference}.augY.{ont_type}.mmap.{chrom}-reads.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        'output/read_subsets/xypar/{sample_info}_{sample}.{reference}.augY.{ont_type}.{chrom}-reads.fasta.gz'
    message: 'DEPRECATED'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: 2
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt ** 3,
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'seqtk subseq {input.reads} {input.read_names} | pigz -p 2 --best > {output}'

