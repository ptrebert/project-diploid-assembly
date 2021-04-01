localrules: master, merge_segment_stats

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
        tmp = 'output/intersect/{sample}_{hap}_isect.tsv'.format(**{'sample': s, 'hap': hap})
        output_files.append(tmp)

output_files.append('output/score_stats/merged_stats.tsv')

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


rule dump_to_bed:
    input:
        'output/segment_align/{sample}_{hap}.wmap-k19.bam'
    output:
        'output/segment_align/{sample}_{hap}.wmap-k19.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools bamtobed -ed -i {input} > {output}'


rule intersect_alignment_with_annotation:
    input:
        aln = 'output/segment_align/{sample}_{hap}.wmap-k19.bed',
        annotation = '/beeond/data/hifiasm_v13_assemblies/coordinates/{sample}_{hap}_segments.bed'
    output:
        'output/intersect/{sample}_{hap}_isect.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -wao -a {input.annotation} -b {input.aln} > {output}'


rule score_bng_segments:
    input:
        'output/intersect/{sample}_{hap}_isect.tsv'
    output:
        'output/score_stats/{sample}_{hap}_stats.tsv'
    run:
        import collections as col
        import pandas as pd

        colors = ['pink', 'green', 'cyan', 'unknown', 'blue',
            'orange', 'purple', 'red']

        def extract_segment_color(seg_name):
            for c in colors:
                if c in seg_name:
                    return c
            return 'blank'

        def compute_match_score(row):
            score = 0
            if row['assm_color'] == row['seg_color']:
                score += 1
                if row['assm_strand'] == row['strand']:
                    score += 1
            return score

        header = ['assm_cluster', 'assm_start', 'assm_end', 'assm_name',
            'seg_cluster', 'seg_start', 'seg_end', 'seg_name',
            'edist', 'strand', 'overlap']

        df = pd.read_csv(input[0], sep='\t', names=header,
            usecols=['assm_name', 'seg_name', 'edist', 'strand', 'overlap'])

        df['assm_color'] = df['assm_name'].apply(lambda x: x.split('_')[-1])
        df['assm_strand'] = df['assm_name'].apply(lambda x: x.split('_')[-2])
        df['strand'].replace({'-': 'minus', '+': 'plus'}, inplace=True)
        df['seg_color'] = df['seg_name'].apply(extract_segment_color)

        df['score'] = df.apply(compute_match_score, axis=1)

        full_match = col.Counter()
        color_match = col.Counter()
        missed = col.Counter()
        unknown_pair = col.Counter()
        mismatch_pair = col.Counter()
        num_counts = col.Counter()

        select_columns = ['assm_color', 'seg_color', 'score']

        for assm_region, region_alns in df.groupby('assm_name'):
            num_counts['segments'] += 1

            row_indexer = region_alns['edist'] == region_alns['edist'].min()
            
            assm_color, seg_color, score = region_alns.loc[row_indexer, select_columns].values[0]
            if score == 2:
                full_match[assm_color] += 1
                num_counts['full_match'] += 1
                num_counts['any_match'] += 1
            elif score == 1:
                color_match[assm_color] += 1
                num_counts['color_match'] += 1
                num_counts['any_match'] += 1
            else:
                if seg_color == 'unknown':
                    unknown_pair[assm_color] += 1
                elif seg_color == 'blank':
                    missed[assm_color] += 1
                elif assm_color != seg_color:
                    mismatch_pair[assm_color] += 1
                else:
                    raise ValueError('unexpected ', assm_color, seg_color, score)

        with open(output[0], 'w') as dump:
            _ = dump.write('num_segments\t{}\n'.format(num_counts['segments']))
            _ = dump.write('num_any_match\t{}\n'.format(num_counts['any_match']))
            _ = dump.write('pct_any_match\t{}\n'.format(round(num_counts['any_match'] / num_counts['segments'] * 100, 1)))
            _ = dump.write('num_full_match\t{}\n'.format(num_counts['full_match']))
            _ = dump.write('pct_full_match\t{}\n'.format(round(num_counts['full_match'] / num_counts['segments'] * 100, 1)))
            _ = dump.write('num_color_match\t{}\n'.format(num_counts['color_match']))
            _ = dump.write('pct_color_match\t{}\n'.format(round(num_counts['color_match'] / num_counts['segments'] * 100, 1)))

            for counter, label in [(unknown_pair,'unknown'),(mismatch_pair,'mismatch'),(missed, 'missed')]:
                
                for k in sorted(counter.keys()):
                    num_k = counter[k]
                    row_label = label + '_' + k
                    _ = dump.write('{}\t{}\n'.format(row_label, num_k))


rule merge_segment_stats:
    input:
        stats = expand('output/score_stats/{sample}_{hap}_stats.tsv', sample=samples, hap=['H1', 'H2'])
    output:
        'output/score_stats/merged_stats.tsv'
    run:
        import pandas as pd

        merged = []
        merged_cols = []

        for stat_file in input:
            filename = stat_file.split('/')[-1].rsplit('_', 1)[0]
            merged_cols.append(filename)
            df = pd.read_csv(stat_file, sep='\t', header=None, index_col=0)
            merged.append(df)

        merged = pd.concat(
            merged,
            axis=1,
            ignore_index=False
        )
        merged.fillna(0., inplace=True)
        merged.columns = merged_cols

        merged.to_csv(merged, sep='\t', header=True, index=True)


rule master:
    input:
        output_files
