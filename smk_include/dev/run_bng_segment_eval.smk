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
pattern = 'output/segment_align/{}_{}.CHM13.wmap-k19.track.bed'
for s in samples:
    tmp = pattern.format(s, 'H1')
    output_files.append(tmp)
    tmp = pattern.format(s, 'H2')
    output_files.append(tmp)
output_files.append('output/score_stats/CHM13_reftoassm_stats.mindist.tsv')
output_files.append('output/score_stats/CHM13_reftoassm_stats.maxscore.tsv')
output_files.append('output/score_stats/CHM13_assmtoref_stats.mindist.tsv')
output_files.append('output/score_stats/CHM13_assmtoref_stats.maxscore.tsv')


BNG_SEGMENT_COLORS = [
    'pink',
    'green',
    'cyan',
    'unknown',
    'blue',
    'orange',
    'purple',
    'red',
    'yellow'
]

BNG_SEGMENT_COLORS_RGB = {
    'pink': '255-192-203',
    'green': '0-255-0',
    'cyan': '0-255-255',
    'unknown': '0-0-0',
    'blue': '0-0-255',
    'orange': '255-165-0',
    'purple': '128-0-128',
    'red': '255-0-0',
    'yellow': '255-255-0'
}


def extract_segment_color(seg_name):
    for c in BNG_SEGMENT_COLORS:
        if c in seg_name:
            return c
    return 'blank'


def split_segments(segment_chain):
    
    split_segments = []
    
    segments = segment_chain.split(';')
    for segment in segments:
        color, start, end = segment.split('-')
        try:
            end = int(end)
        except ValueError:
            start, end, color = segment.split('-')
            start = int(start)
            end = int(end)
        else:
            start = int(start)
        if start > end:
            start, end = end, start
        if end - start > 100000:
            raise ValueError('Likely annotation error: {} - {} ({})'.format(start, end, color))
        color = color.strip()
        regular_color = extract_segment_color(color)
        try:
            color_rgb = BNG_SEGMENT_COLORS_RGB[regular_color]
        except KeyError:
            raise ValueError('Cannot process segment name: {} / {} / {}'.format(segment, color, regular_color))
        split_segments.append((start, end, color, color_rgb))

    return split_segments


rule flatten_bng_segments_table:
    """
    '/beeond/projects/bng_eval/annotation/20210327_1p36_HiFiAsm_SegmentInfo.xlsx'
    """
    input:
        '/beeond/projects/bng_eval/annotation/20210412_1p36_HiFiAsm_SegmentInfo.xlsx'
    output:
        'input/table/20210412_1p36_HiFiAsm_SegmentInfo.flat.tsv'
    run:
        import pandas as pd

        table_header = ['sample', 'cluster_id', 'cluster_length', 'sample_id',
            'haplotype', 'contig_id', 'orientation', 'segments']
        
        df = pd.read_excel(
            input[0],
            header=0,
            names=table_header,
        )
        df.drop('sample_id', axis=1, inplace=True)

        # specific for T2T reference
        df['cluster_id'].fillna('chr1', inplace=True)
        df['cluster_length'].fillna('248387497', inplace=True)
        df['cluster_length'] = df['cluster_length'].astype(int)
        df['haplotype'].fillna('H0', inplace=True)
        df['haplotype'].replace({'h1': 'H1', 'h2': 'H2'}, inplace=True)
        df['contig_id'].fillna(1, inplace=True)
        df['orientation'].fillna(0, inplace=True)
        df['orientation'].replace({'direct': 1, 'inverted': -1}, inplace=True)

        new_rows = []
        for idx, row in df.iterrows():
            assembly_segments = split_segments(row['segments'])
            for (start, end, color, color_rgb) in assembly_segments:
                new_rows.append(
                    (
                        row['sample'].strip(),
                        row['haplotype'].strip(),
                        row['cluster_id'].strip(),
                        row['cluster_length'],
                        row['contig_id'],
                        row['orientation'],
                        start,
                        end,
                        color,
                        color_rgb
                    )
                )

        df = pd.DataFrame(
            new_rows,
            columns=[
                'sample',
                'haplotype',
                'cluster_id',
                'cluster_length',
                'contig_id',
                'orientation',
                'start',
                'end',
                'color',
                'color_rgb'
            ]
        )
        df['contig_id'] = df['contig_id'].astype(int)

        # red segments just indicate region boundaries
        df = df.loc[df['color'] != 'red', :].copy()

        df.to_csv(
            output[0],
            sep='\t',
            header=True,
            index=False,
            encoding='ascii'
        )
    # END OF RUN BLOCK


rule extract_reference_segments:
    input:
        table = 'input/table/20210412_1p36_HiFiAsm_SegmentInfo.flat.tsv'
    output:
        tig_names = 'output/tig_names/T2Tv1_38p13Y_chm13.tigs.txt',
        segment_bed = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.segments.bed',
        track_bed = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.track.bed'
    run:
        import pandas as pd

        df = pd.read_csv(input.table, sep='\t', header=0, index_col=None)

        with open(output.tig_names, 'w') as tig_names:
            _ = tig_names.write('chr1\n')

        segments = df.loc[df['sample'] == 'chm13', ['cluster_id', 'start', 'end', 'color', 'color_rgb', 'orientation']].copy()
        segments['orientation'] = 'plus'
        segments['name'] = 'CHM13_1p3613_' + \
            segments[['start', 'end', 'color', 'color_rgb', 'orientation']].astype(str).apply(lambda row: '_'.join(row), axis=1)
        segments.sort_values(['start', 'end'], inplace=True)
        segments.to_csv(
            output.segment_bed,
            sep='\t',
            header=False,
            index=False,
            columns=['cluster_id', 'start', 'end', 'name']
        )

        track = df.loc[df['sample'] == 'chm13', ['cluster_id', 'start', 'end', 'color', 'color_rgb']].copy()
        track['name'] = track['color']
        track['score'] = (track['end'] - track['start']).astype(int)
        track['score'].clip(0, 1000, inplace=True)
        track['strand'] = '.'
        track['thickStart'] = track['start']
        track['thickEnd'] = track['end']
        track['color'] = track['color_rgb'].map(lambda x: ','.join(x.split('-')))

        with open(output.track_bed, 'w') as track_bed:
            _ = track_bed.write('track name="CHM13_1p3613_segments" itemRgb="On"\n')
            track.to_csv(
                track_bed,
                sep='\t',
                header=False,
                index=False,
                columns=['cluster_id', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'color']
            )


rule extract_reference_segment_sequences:
    input:
        segment_bed = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.segments.bed',
        reference_fasta = '/beeond/data/hifiasm_v13_assemblies/T2Tv1_38p13Y_chm13.fasta',
    output:
        'output/segment_sequences/T2Tv1_38p13Y_chm13.segments.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'bedtools getfasta -fi {input.reference_fasta} -bed {input.segment_bed} -nameOnly > {output}'


rule extract_assembly_segments:
    input:
        table = 'input/table/20210412_1p36_HiFiAsm_SegmentInfo.flat.tsv'
    output:
        tig_names = 'output/tig_names/{sample}_{hap}_tigs.txt',
        segment_bed = 'output/segment_coordinates/{sample}_{hap}.segments.bed',
        track_bed = 'output/segment_coordinates/{sample}_{hap}.track.bed'
    run:
        import pandas as pd

        df = pd.read_csv(input.table, sep='\t', header=0, index_col=None)

        row_indexer = (df['sample'] == wildcards.sample) & (df['haplotype'] == wildcards.hap)

        segments = df.loc[row_indexer, ['cluster_id', 'start', 'end', 'color', 'color_rgb', 'orientation']].copy()
        segments['orientation'].replace({1: 'plus', -1: 'minus'}, inplace=True)

        tigs = set(segments['cluster_id'].values)
        with open(output.tig_names, 'w') as tig_names:
            _ = tig_names.write('\n'.join(sorted(tigs)) + '\n')

        segments['name'] = '{}_{}_'.format(wildcards.sample, wildcards.hap) + \
            segments[['cluster_id', 'start', 'end', 'color', 'color_rgb', 'orientation']].astype(str).apply(lambda row: '_'.join(row), axis=1)
        segments.sort_values(['start', 'end'], inplace=True)
        segments.to_csv(
            output.segment_bed,
            sep='\t',
            header=False,
            index=False,
            columns=['cluster_id', 'start', 'end', 'name']
        )

        segments['score'] = (segments['end'] - segments['start']).astype(int)
        segments['score'].clip(0, 1000, inplace=True)
        segments['strand'] = '.'
        segments['thickStart'] = segments['start']
        segments['thickEnd'] = segments['end']
        segments['color'] = segments['color_rgb'].map(lambda x: ','.join(x.split('-')))

        with open(output.track_bed, 'w') as track_bed:
            _ = track_bed.write('track name="{}_{}_segments" itemRgb="On"\n'.format(wildcards.sample, wildcards.hap))
            segments.to_csv(
                track_bed,
                sep='\t',
                header=False,
                index=False,
                columns=['cluster_id', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'color']
            )


rule extract_hap1_segment_sequences:
    input:
        segment_bed = 'output/segment_coordinates/{sample}_H1.segments.bed',
        reference_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta'
    output:
        'output/segment_sequences/{sample}_H1.segments.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'bedtools getfasta -fi {input.reference_fasta} -bed {input.segment_bed} -nameOnly > {output}'


rule extract_hap2_segment_sequences:
    input:
        segment_bed = 'output/segment_coordinates/{sample}_H2.segments.bed',
        reference_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta'
    output:
        'output/segment_sequences/{sample}_H2.segments.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'bedtools getfasta -fi {input.reference_fasta} -bed {input.segment_bed} -nameOnly > {output}'


rule extract_reference_contigs:
    input:
        tig_names = 'output/tig_names/T2Tv1_38p13Y_chm13.tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/T2Tv1_38p13Y_chm13.fasta',
    output:
        assm_tigs = 'output/tig_sequences/T2Tv1_38p13Y_chm13_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'


rule extract_hap1_contigs:
    input:
        tig_names = 'output/tig_names/{sample}_H1_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta',
    output:
        assm_tigs = 'output/tig_sequences/{sample}_H1_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'


rule extract_hap2_contigs:
    input:
        tig_names = 'output/tig_names/{sample}_H2_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta',
    output:
        assm_tigs = 'output/tig_sequences/{sample}_H2_tigs.fasta'
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


rule count_reference_kmers:
    input:
        fasta = '/beeond/data/hifiasm_v13_assemblies/T2Tv1_38p13Y_chm13.fasta'
    output:
        kmer_db = directory('output/kmer/T2Tv1_38p13Y_chm13.k19.db/'),
        rep_kmer = 'output/kmer/T2Tv1_38p13Y_chm13.k19.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer/T2Tv1_38p13Y_chm13.count-dump.rsrc'
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


rule align_assembly_segments_to_reference:
    input:
        tigs = 'output/tig_sequences/T2Tv1_38p13Y_chm13_tigs.fasta',
        segments = 'output/segment_sequences/{sample}_{hap}.segments.fasta',
        kmers = 'output/kmer/T2Tv1_38p13Y_chm13.k19.rep-grt09998.txt'
    output:
        'output/segment_align/{sample}_{hap}.CHM13.wmap-k19.bam'
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


rule align_reference_segments_to_assembly:
    input:
        tigs = 'output/tig_sequences/{sample}_{hap}_tigs.fasta',
        segments = 'output/segment_sequences/T2Tv1_38p13Y_chm13.segments.fasta',
        kmers = 'output/kmer/{sample}_{hap}.k19.rep-grt09998.txt'
    output:
        'output/segment_align/CHM13.{sample}_{hap}.wmap-k19.bam'
    benchmark:
        'rsrc/output/segment_align/CHM13.{sample}_{hap}.wmap.rsrc'
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
        'output/segment_align/{filename}.wmap-k19.bam'
    output:
        'output/segment_align/{filename}.wmap-k19.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools bamtobed -ed -i {input} > {output}'


rule make_bed_dump_to_track:
    input:
        'output/segment_align/{filename}.wmap-k19.bed'
    output:
        'output/segment_align/{filename}.wmap-k19.track.bed'
    run:
        import pandas as pd
        header = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        df = pd.read_csv(input[0], sep='\t', header=None, names=header)
        df['thickStart'] = df['start']
        df['thickEnd'] = df['end']
        df['color'] = df['name'].map(lambda x: ','.join(x.split('_')[-2].split('-')))

        with open(output[0], 'w') as dump:
            _ = dump.write('track name="{}_segment_alignments" itemRgb="On"\n'.format(wildcards.filename))
            df.to_csv(
                dump,
                sep='\t',
                header=False,
                index=False
            )


rule intersect_annotation_ref_to_assm:
    input:
        aln = 'output/segment_align/CHM13.{sample}_{hap}.wmap-k19.bed',
        annotation = 'output/segment_coordinates/{sample}_{hap}.segments.bed'
    output:
        'output/intersect/CHM13.{sample}_{hap}.isect.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -wao -a {input.annotation} -b {input.aln} > {output}'


rule intersect_annotation_assm_to_ref:
    input:
        aln = 'output/segment_align/{sample}_{hap}.CHM13.wmap-k19.bed',
        annotation = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.segments.bed'
    output:
        'output/intersect/{sample}_{hap}.CHM13.isect.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -wao -a {input.annotation} -b {input.aln} > {output}'


def compute_match_score(row):
    score = 0
    if row['assm_color'] == row['seg_color']:
        score += 1
        if row['assm_strand'] == row['strand']:
            score += 1
    return score


def extract_scored_color_pair(alignments, scoring, select_columns):

    if scoring == 'mindist':
        row_indexer = alignments['edist'] == alignments['edist'].min()    
        assm_color, seg_color, score = alignments.loc[row_indexer, select_columns].values[0]
    elif scoring == 'maxscore':
        tmp = alignments.sort_values(['score', 'edist'], ascending=[False, True], inplace=False)
        row_indexer = tmp.index[0]    
        assm_color, seg_color, score = alignments.loc[row_indexer, select_columns].values
    else:
        raise ValueError('unexpected scoring: {}'.format(scoring))
    return assm_color, seg_color, score
 

rule score_bng_segments:
    input:
        'output/intersect/{assembly_1}.{assembly_2}.isect.tsv'
    output:
        'output/score_stats/{assembly_1}.{assembly_2}.stats.{scoring}.tsv'
    run:
        import collections as col
        import pandas as pd    

        BNG_SEGMENT_INTERSECT_HEADER = ['assm_cluster', 'assm_start', 'assm_end',
            'assm_name', 'seg_cluster', 'seg_start', 'seg_end', 'seg_name',
            'edist', 'strand', 'overlap']

        df = pd.read_csv(input[0], sep='\t', names=BNG_SEGMENT_INTERSECT_HEADER,
            usecols=['assm_name', 'seg_name', 'edist', 'strand', 'overlap'])

        df['assm_color'] = df['assm_name'].apply(extract_segment_color)
        if wildcards.assembly_1 == 'CHM13':
            df['assm_strand'] = df['assm_name'].apply(lambda x: x.split('_')[-1])
        else:
            df['assm_strand'] = df['seg_name'].apply(lambda x: x.split('_')[-1])
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
            
            assm_color, seg_color, score = extract_scored_color_pair(
                region_alns,
                wildcards.scoring,
                select_columns
            )
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
    # END OF RUN BLOCK

rule merge_segment_stats_reference_to_assembly:
    input:
        stats = expand('output/score_stats/CHM13.{sample}_{hap}.stats.{{stattype}}.tsv', sample=samples, hap=['H1', 'H2'])
    output:
        'output/score_stats/CHM13_reftoassm_stats.{stattype}.tsv'
    run:
        import pandas as pd

        merged = []
        merged_cols = []

        for stat_file in input:
            filename = stat_file.split('/')[-1].split('.')[1]
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

        merged.to_csv(output[0], sep='\t', header=True, index=True, index_label='statistic')


rule merge_segment_stats_assembly_to_reference:
    input:
        stats = expand('output/score_stats/{sample}_{hap}.CHM13.stats.{{stattype}}.tsv', sample=samples, hap=['H1', 'H2'])
    output:
        'output/score_stats/CHM13_assmtoref_stats.{stattype}.tsv'
    run:
        import pandas as pd

        merged = []
        merged_cols = []

        for stat_file in input:
            filename = stat_file.split('/')[-1].split('.')[0]
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

        merged.to_csv(output[0], sep='\t', header=True, index=True, index_label='statistic')


rule merge_reference_segments:
    input:
        'output/segment_coordinates/T2Tv1_38p13Y_chm13.segments.bed'
    output:
        temp('output/segment_coordinates/T2Tv1_38p13Y_chm13.merged.bed')
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools merge -i {input} | cut -f 1-3 > {output}'


rule define_extended_region:
    input:
        'output/segment_coordinates/T2Tv1_38p13Y_chm13.segments.bed'
    output:
        temp('output/segment_coordinates/T2Tv1_38p13Y_chm13.ext-region.bed')
    run:
        import pandas as pd
        df = pd.read_csv(
            input[0],
            sep='\t',
            names=['chrom', 'start', 'end', 'name'],
            usecols=['chrom', 'start', 'end']
        )

        chrom = df['chrom'].values[0]
        start = df['start'].min() - 100000
        end = df['end'].max() + 100000

        with open(output[0], 'w') as dump:
            _ = dump.write('{}\t{}\t{}\n'.format(chrom, int(start), int(end)))


rule subtract_reference_segments:
    input:
        region = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.ext-region.bed',
        segments = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.merged.bed'
    output:
        temp('output/segment_coordinates/T2Tv1_38p13Y_chm13.subtract.bed')
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools subtract -a {input.region} -b {input.segments} > {output}'


rule build_segment_complements:
    input:
        'output/segment_coordinates/T2Tv1_38p13Y_chm13.subtract.bed'
    output:
        bed = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.complements.bed',
        track = 'output/segment_coordinates/T2Tv1_38p13Y_chm13.complements.track.bed',
    run:
        track_line = 'track name="CHM13_1p3613_complements" itemRgb="On"'
        name_template = 'CHM13_1p3613_{}_{}_grey-{}_128-128-128_plus'
        bed_template = '{}\t{}\t{}\t{}'
        track_template = '{}\t{}\t{}\t{}\t1000\t.\t{}\t{}\t128,128,128'

        bed_rows = []
        track_rows = []

        with open(input[0], 'r') as bedfile:
            for ln, line in enumerate(bedfile, start=1):
                chrom, start, end = line.strip().split()
                row_name = name_template.format(start, end, ln)

                bed_row = bed_template.format(chrom, start, end, row_name)
                bed_rows.append(bed_row)

                track_row = track_template.format(chrom, start, end, row_name, start, end)
                track_rows.append(track_row)
        
        with open(output.bed, 'w') as dump:
            _ = dump.write('\n'.join(bed_rows) + '\n')

        with open(input.track, 'w') as dump:
            _ = dump.write(track_line + '\n')
            _ = dump.write('\n'.join(track_rows) + '\n')
        



rule master:
    input:
        output_files
