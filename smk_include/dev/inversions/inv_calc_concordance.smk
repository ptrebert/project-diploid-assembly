
import pathlib as pl

ANNOTATION_REFERENCE_SEGMENTS = '/home/local/work/code/github/project-diploid-assembly/annotation/bng_segments/20211006_chm13_1p36.13_allSegmentsV2.tsv'
ANNOTATION_ASSEMBLY_SEGMENTS = '/home/local/work/code/github/project-diploid-assembly/annotation/bng_segments/20211006_hifiasmv13dev_1p36.13_allSegmentsV2.tsv'
REFERENCE_FASTA = '/home/local/work/data/references/T2Tv11_38p13Y_chm13.fasta'
ASSEMBLY_FASTA_FOLDER = '/home/local/work/data/hgsvc_2021/v13dev/phased_assm'


rule run_all:
    input:
        'output/concordance/segments_to_segments.1p3613.tsv',
        'output/concordance/segments_to_reference.1p3613.RO30.tsv',
        'output/concordance/segments_to_reference.1p3613.RO50.tsv'


def collect_sample_names():

    fasta = pl.Path(ASSEMBLY_FASTA_FOLDER).glob('*.fasta')
    samples = set()
    for f in fasta:
        smp = f.name.split('_')[0]
        samples.add(smp)
    assert samples
    return sorted(samples)

SAMPLES = collect_sample_names()


def load_sequence(fasta, chromosome):

    seq = ''
    buffer = False
    with open(fasta, 'r') as fa:
        for line in fa:
            if line.startswith('>'):
                if line.strip() == f'>{chromosome}':
                    buffer = True
                    continue
                elif seq:
                    break
                else:
                    buffer = False
            elif buffer:
                seq += line.strip()
            else:
                continue
    if not seq:
        raise ValueError('Reference sequence empty')
    return seq


def extract_segment_info(row):
    segment = row['segment']
    order_lut = {
        'one': '1',
        'two': '2',
        'three': '3',
        'four': '4'
    }
    if segment.endswith('inv'):
        orientation = 'invert'
    else:
        orientation = 'direct'
    segment = segment.split('inv')[0]
    parts = segment.split('-')
    if len(parts) == 1:
        color = parts[0]
        order = '1'
    else:
        color = parts[0]
        order = order_lut[parts[1]]
    return color, order, orientation


def compute_region_endpoints(min_start, max_end):
    region_span = max_end - min_start
    assert 1e5 < region_span < 1e6
    blunt_factor = 1e5
    extended_span = region_span // blunt_factor * blunt_factor + blunt_factor
    region_mid = min_start + region_span // 2
    region_left_end = int(region_mid - extended_span // 2)
    region_right_end = int(region_mid + extended_span // 2)
    return region_left_end, region_right_end


rule extract_roi_from_reference:
    """
    This rule assumes that BNG coordinates
    are right-inclusive
    """
    input:
        ref_fasta = REFERENCE_FASTA,
        ref_segments = ANNOTATION_REFERENCE_SEGMENTS
    output:
        ref_region = 'output/regions/T2Tv11_roi.{region}.fasta',
        ref_segments = 'output/segments/T2Tv11_segments.{region}.fasta',
        ref_colors = 'output/segments/T2Tv11_segment-colors.{region}.txt'
    run:
        import pandas as pd
        right_offset = 1

        assert wildcards.region == '1p3613'
        ref_chrom = 'chr' + wildcards.region[0]  # breaks for chr10+

        segments = pd.read_csv(
            input.ref_segments,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end', 'segment']
        )
        segments[['color', 'order', 'orientation']] = segments.apply(extract_segment_info, axis=1, result_type='expand')
        with open(output.ref_colors, 'w') as dump:
            colors = set(segments['color'].values)
            _ = dump.write('\n'.join(sorted(colors)) + '\n')

        chrom_seq = load_sequence(input.ref_fasta, ref_chrom)
        ref_left_end, ref_right_end = compute_region_endpoints(segments['start'].min(), segments['end'].max() + right_offset)
        roi_header = f'>T2Tv11_{wildcards.region}_{ref_chrom}_{ref_left_end}-{ref_right_end}'
        with open(output.ref_region, 'w') as fasta:
            _ = fasta.write(roi_header + '\n')
            _ = fasta.write(chrom_seq[ref_left_end:ref_right_end] + '\n')
        
        for row in segments.itertuples(index=False):
            seg_header = '_'.join([
                'T2Tv11',
                wildcards.region,
                ref_chrom,
                str(row.start) + '-' + str(row.end),
                row.color + '-' + row.order + '-' + row.orientation + '-plus'
            ])
            seg_sequence = chrom_seq[row.start:row.end+right_offset]
            with open(output.ref_segments, 'a') as fasta:
                _ = fasta.write(f'>{seg_header}\n{seg_sequence}\n')
    # END OF RUN BLOCK


def select_assembly_input(wildcards):

    if wildcards.hap == 'h1':
        ext_hap = 'h1-un'
    elif wildcards.hap == 'h2':
        ext_hap = 'h2-un'
    else:
        raise

    assm_fasta = list(pl.Path(ASSEMBLY_FASTA_FOLDER).glob(f'{wildcards.sample}*{ext_hap}*.fasta'))
    if len(assm_fasta) != 1:
        raise ValueError(assm_fasta)
    return str(assm_fasta[0])


rule extract_segments_from_assemblies:
    """
    This rule assumes that BNG coordinates
    are right-inclusive
    """
    input:
        assm_segments = ANNOTATION_ASSEMBLY_SEGMENTS,
        assm_fasta = select_assembly_input
    output:
        assm_segments = 'output/segments/{sample}_segments.{hap}.fasta',
        assm_colors = 'output/segments/{sample}_segment-colors.{hap}.txt',
        segment_list = 'output/segments/{sample}_segments.{hap}.txt',
    run:
        import pandas as pd
        right_offset = 1

        # short header for version 2021-08-30
        # table_header = ['sample', 'hap', 'contig', 'ctg_orient', 'start', 'end', 'color', 'seg_orient']
        # use_columns = list(range(8))
        
        #  short header for version 2021-10-06
        table_header = ['sample', 'hap', 'contig', 'ctg_orient', 'start', 'end', 'color', 'seg_orient', 'seg_state', 'ctg_length']
        use_columns = list(range(10))

        segments = pd.read_csv(
            input.assm_segments,
            sep='\t',
            usecols=use_columns,
            header=0,
            names=table_header
        )
        segments['ctg_orient'] = segments['ctg_orient'].replace({'+': 'plus', '-': 'minus', '.': 'unknown'})
        segments['seg_orient'] = segments['seg_orient'].apply(lambda x: x.rstrip('ed'))
        segments = segments.loc[(segments['sample'] == wildcards.sample) & (segments['hap'] == wildcards.hap), :]
        segments.sort_values(['contig', 'start', 'end'], ascending=True, inplace=True)
        segments['order'] = (segments.groupby(['contig', 'color']).cumcount()+1).astype(str)
        with open(output.assm_colors, 'w') as dump:
            colors = set(segments['color'].values)
            _ = dump.write('\n'.join(sorted(colors)) + '\n')
        
        segment_names = []
        for ctg, ctg_segments in segments.groupby('contig'):
            ctg_sequence = load_sequence(input.assm_fasta, ctg)
            for row in ctg_segments.itertuples(index=False):
                seg_header = '_'.join([
                    wildcards.sample,
                    wildcards.hap,
                    ctg.replace('_', '-'),
                    str(row.start) + '-' + str(row.end),
                    row.color + '-' + row.order + '-' + row.seg_orient + '-' + row.ctg_orient
                ])
                segment_names.append(seg_header)
                seg_sequence = ctg_sequence[row.start:row.end+right_offset]
                with open(output.assm_segments, 'a') as fasta:
                    _ = fasta.write(f'>{seg_header}\n{seg_sequence}\n')
        with open(output.segment_list, 'w') as listing:
            _ = listing.write('\n'.join(sorted(segment_names)) + '\n')
    # END OF RUN BLOCK


rule merge_segments:
    input:
        ref_segments = 'output/segments/T2Tv11_segments.{region}.fasta',
        assm_segments = expand(
            'output/segments/{sample}_segments.{hap}.fasta',
            sample=SAMPLES,
            hap=['h1', 'h2']
        ),
        name_lists = expand(
            'output/segments/{sample}_segments.{hap}.txt',
            sample=SAMPLES,
            hap=['h1', 'h2']
        ),
    output:
        fasta = 'output/segments/merged_segments.{region}.fasta',
        stats = 'output/segments/merged_segments.{region}.stats.tsv',
        lengths = 'output/segments/merged_segments.{region}.lengths.tsv',
        names = 'output/segments/merged_segments.{region}.names.txt',
    run:
        import collections as col
        import io
        import math
        import numpy as np

        def nck(n, k):
            return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))
        
        out_buffer = io.StringIO()
        colored_segments = col.Counter()
        segment_lengths = col.defaultdict(list)
        for seg_file in [input.ref_segments] + input.assm_segments:
            with open(seg_file, 'r') as fasta:
                while 1:
                    header = fasta.readline().strip()
                    if not header:
                        break
                    seg_color = header.split('_')[-1].split('-')[0]
                    seg_start, seg_end = header.split('_')[-2].split('-')
                    seg_length = int(seg_end) - int(seg_start)
                    colored_segments[seg_color] += 1
                    segment_lengths[seg_color].append(seg_length)
                    _ = out_buffer.write(header + '\n')
                    _ = out_buffer.write(fasta.readline())
        with open(output.fasta, 'w') as fasta:
            _ = fasta.write(out_buffer.getvalue())
        with open(output.stats, 'w') as table:
            for color in sorted(colored_segments.keys()):
                count = colored_segments[color]
                try:
                    _ = table.write(f'{color}\t{count}\t{nck(count, 2)}\n')
                except ValueError:
                    assert count < 2
                    _ = table.write(f'{color}\t{count}\t-1\n')
        with open(output.lengths, 'w') as table:
            for color in sorted(segment_lengths.keys()):
                values = np.sort(np.array(segment_lengths[color], dtype=np.int32))
                _ = table.write(f'{color}_NUM_total\t{values.size}\n')
                _ = table.write(f'{color}_LEN_minimum\t{values.min()}\n')
                _ = table.write(f'{color}_LEN_median\t{values[values.size//2]}\n')
                _ = table.write(f'{color}_LEN_mean\t{round(values.mean(), 0)}\n')
                _ = table.write(f'{color}_LEN_stddev\t{round(values.std(), 0)}\n')
                _ = table.write(f'{color}_LEN_maximum\t{values.max()}\n')

        out_buffer = io.StringIO()
        for listing in input.name_lists:
            with open(listing, 'r') as fd:
                content = fd.read()
                out_buffer.write(content)
        
        with open(output.names, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())
    # END OF RUN BLOCK


rule mmap_all_vs_all_segments:
    input:
        fasta = 'output/segments/merged_segments.{region}.fasta',
    output:
        paf = 'output/all_vs_all/merged_segments.{region}.ava.{secsim}.paf'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    threads: 4
    params:
        p = lambda wildcards: str(round(int(wildcards.secsim)/100,2))
    shell:
        'minimap2 -o {output.paf} -t {threads} -X -x ava-pb -p {params.p} -N 200 {input.fasta} {input.fasta}'


rule mmap_segments_to_reference:
    input:
        ref_region = 'output/regions/T2Tv11_roi.{region}.fasta',
        segments = 'output/segments/merged_segments.{region}.fasta',
    output:
        paf = 'output/seg_to_ref/merged_segments.{region}.ref.paf'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    threads: 4
    shell:
        'minimap2 -o {output.paf} -t {threads} -x map-hifi --secondary=no {input.ref_region} {input.segments}'


PAF_HEADER = [
    'qname',
    'qlength',
    'qstart',
    'qend',
    'map_orient',
    'tname',
    'tlength',
    'tstart',
    'tend',
    'matches',
    'block_length',
    'mapq',
    'map_type',
    'foo',
    'bar',
    'baz',
    'seqdiv',
    'boo'
]


rule compute_reference_concordance:
    """
    for some weird reason, intevalindex.overlaps
    throws an NotImplementedError, so the below
    uses brute force...
    """
    input:
        paf = 'output/seg_to_ref/merged_segments.{region}.ref.paf',
        segments = 'output/segments/merged_segments.{region}.stats.tsv',
        lengths = 'output/segments/merged_segments.{region}.lengths.tsv',
        names = 'output/segments/merged_segments.{region}.names.txt',
    output:
        table = 'output/concordance/segments_to_reference.{region}.RO{ro}.tsv'
    run:
        import pandas as pd
        import collections as col

        ro_overlap_fraction = round(float(wildcards.ro) / 100, 2)

        segment_names = set(open(input.names, 'r').read().strip().split())

        stats = col.Counter()
        with open(input.segments, 'r') as counts:
            for line in counts:
                _, c, _ = line.split()
                stats['total_segments'] += int(c)


        paf = pd.read_csv(
            input.paf, sep='\t', header=None,
            names=PAF_HEADER,
            usecols=[0,1,2,3,4,7,8,16]
        )

        paf['seqdiv'] = paf['seqdiv'].apply(lambda x: float(x.split(':')[-1]))
        paf['sample'] = paf['qname'].apply(lambda x: x.split('_')[0])
        paf['qcolor'] = paf['qname'].apply(lambda x: x.split('_')[-1].split('-')[0])
        paf['qorient'] = paf['qname'].apply(lambda x: x.split('_')[-1].split('-')[-2])
        paf['qcontig'] = paf['qname'].apply(lambda x: x.split('_')[-1].split('-')[-1])
        paf['hap'] = paf['qname'].apply(lambda x: x.split('_')[1])
        paf['map_orient'] = paf['map_orient'].replace({'+': 'plus', '-': 'minus'})
        
        paf.sort_values(['tstart', 'tend'], inplace=True)

        t2t_segments = paf.loc[paf['sample'] == 'T2Tv11', :].copy()
        stats['total_segments'] -= t2t_segments.shape[0]
        assm_segments = paf.loc[paf['sample'] != 'T2Tv11', :].copy()

        aligned = set(assm_segments['qname'].values)
        unaligned = segment_names - aligned
        
        split_align = []
        score_orientation = {
            ('invert', 'minus'): 1,
            ('direct', 'plus'): 1
        }
        for segment, alignments in assm_segments.groupby('qname'):
            color, _, expected_orientation, ctg_orientation = segment.split('_')[-1].split('-')
            is_split = False
            if alignments.shape[0] > 1:
                split_align.append(segment)
                # split alignment - discordant
                # stats['discordant_total'] += 1
                # stats['discordant_split'] += 1
                # continue
                # select alignment with smallest seq. div.
                alignments = alignments.loc[alignments['seqdiv'] == alignments['seqdiv'].min(), :]
                if alignments.shape[0] > 1:
                    # two alignments with identical seq. div.
                    stats['discordant_total'] += 1
                    stats['discordant_split'] += 1
                    continue
                stats['split_align'] += 1
                is_split = True

            map_orient = alignments['map_orient'].values[0]
            select_start = alignments['tstart'].values[0] < t2t_segments['tend']
            select_end = t2t_segments['tstart'] < alignments['tend'].values[0]
            selector = select_start & select_end
            t2t_match = t2t_segments.loc[selector, :]
            if t2t_match.empty:
                stats['discordant_total'] += 1
                stats['discordant_no-match'] += 1
                continue
            if t2t_match.shape[0] > 1:
                stats['multi_ref_hit'] += 1
                # T2T segments are close to each other
                # check if overlap is unambigous
                assm_start = alignments['tstart'].values[0]
                assm_end = alignments['tend'].values[0]
                assm_length = assm_end - assm_start
                overlaps = []
                for row in t2t_match.itertuples(index=False):
                    t2t_length = row.tend - row.tstart
                    ovl = min(assm_end, row.tend) - max(assm_start, row.tstart)
                    color = row.qcolor
                    overlaps.append(
                        (
                            ovl,
                            ovl/assm_length > ro_overlap_fraction,
                            ovl/t2t_length > ro_overlap_fraction,
                            color
                        )
                    )
                overlaps = [ovl for ovl in overlaps if ovl[1] & ovl[2]]
                if len(overlaps) > 1:
                    stats['discordant_total'] += 1
                    stats['multi_ref_hit_ambig'] += 1
                    continue
                elif len(overlaps) == 0:
                    stats['discordant_total'] += 1
                    stats['multi_ref_hit_weak'] += 1
                    continue
                else:
                    stats['multi_ref_hit_uniq'] += 1
                    t2t_color = overlaps[0][3]
            else:
                t2t_color = t2t_match['qcolor'].values[0]
            color_match = t2t_color == color
            orient_match = score_orientation.get((expected_orientation, map_orient), 0) > 0
            if color_match and orient_match and not is_split:
                stats[('HIT', expected_orientation, ctg_orientation, map_orient, 'contig')] += 1
                stats['concordant_total'] += 1
                stats['concordant_full_ctg_total'] += 1
                stats['concordant_ctg_total'] += 1
                stats['concordant_full_total'] += 1
            elif color_match and orient_match:
                stats[('HIT', expected_orientation, ctg_orientation, map_orient, 'split')] += 1
                stats['concordant_total'] += 1
                stats['concordant_full_split_total'] += 1
                stats['concordant_split_total'] += 1
                stats['concordant_full_total'] += 1
            elif color_match and not is_split:
                # collect orientation statistics to check for biases
                stats[('MISO', expected_orientation, ctg_orientation, map_orient, 'contig')] += 1
                stats['concordant_total'] += 1
                stats['concordant_miso_ctg_total'] += 1
                stats['concordant_ctg_total'] += 1
                stats['concordant_miso_total'] += 1
            elif color_match and is_split:
                # collect orientation statistics to check for biases
                stats[('MISO', expected_orientation, ctg_orientation, map_orient, 'split')] += 1
                stats['concordant_total'] += 1
                stats['concordant_miso_split_total'] += 1
                stats['concordant_split_total'] += 1
                stats['concordant_miso_total'] += 1
            else:
                stats['discordant_total'] += 1
        stats['aligned_segments'] = stats['concordant_total'] + stats['discordant_total']
        pct_values = [
            'concordant_total',
            'concordant_full_ctg_total',
            'concordant_full_split_total',
            'concordant_full_total',
            'concordant_miso_ctg_total',
            'concordant_miso_split_total',
            'concordant_miso_total',
            'concordant_ctg_total',
            'concordant_split_total',
            'discordant_total',
            'aligned_segments',
            'split_align'
        ]
        with open(output[0], 'w') as table:
            for key in sorted([k for k in stats.keys() if isinstance(k, str)]):
                value = stats[key]
                _ = table.write(f'{key}\t{value}\n')
                if key in pct_values:
                    pct_value = round(value/stats['total_segments']*100,1)
                    _ = table.write(f'{key}_pct\t{pct_value}\n')
            for key in sorted([k for k in stats.keys() if isinstance(k, tuple)]):
                out_key = f'{key[0]}_SEG-{key[1]}_CTG-{key[2]}_ALN-{key[3]}_STATE-{key[4]}'
                value = stats[key]
                _ = table.write(f'{out_key}\t{value}\n')
            for idx, split in enumerate(split_align, start=1):
                _ = table.write(f'split_{idx}\t{split}\n')
            for idx, segname in enumerate(sorted(unaligned), start=1):
                _ = table.write(f'unaligned_{idx}\t{segname}\n')
    # END OF RUN BLOCK


rule compute_all_vs_all_concordance:
    input:
        paf = 'output/all_vs_all/merged_segments.{region}.ava.080.paf',
        segments = 'output/segments/merged_segments.{region}.stats.tsv',
        lengths = 'output/segments/merged_segments.{region}.lengths.tsv',
    output:
        table = 'output/concordance/segments_to_segments.{region}.tsv'
    run:
        import pandas as pd
        import numpy as np
        import collections as col

        pairs = col.Counter()
        with open(input.segments, 'r') as counts:
            for line in counts:
                color, _, nck = line.split()
                pairs[color] = int(nck)

        orientation_scoring = {
            ('direct', 'direct', 'plus'): 1,
            ('direct', 'invert', 'minus'): 1,
            ('invert', 'direct', 'minus'): 1,
            ('direct', 'direct', 'minus'): 0,
            ('direct', 'invert', 'plus'): 0,
            ('invert', 'direct', 'plus'): 0,
            ('invert', 'invert', 'minus'): 0,
            ('invert', 'invert', 'plus'): 1
        }

        def score_orientation(row, scoring):
            qorient = row['qorient']
            torient = row['torient']
            map_orient = row['map_orient']
            try:
                score = scoring[(qorient, torient, map_orient)]
            except KeyError:
                print(row)
                raise
            return score
            

        paf = pd.read_csv(
            input.paf, sep='\t', header=None,
            names=PAF_HEADER[:15] + PAF_HEADER[16:],
            usecols=[0,1,2,3,4,5,6,7,8,9,10,15]
        )
        paf['seqdiv'] = paf['seqdiv'].apply(lambda x: float(x.split(':')[-1]))
        paf['qsample'] = paf['qname'].apply(lambda x: x.split('_')[0])
        paf['tsample'] = paf['tname'].apply(lambda x: x.split('_')[0])
        paf['map_orient'] = paf['map_orient'].replace({'+': 'plus', '-': 'minus'})
        paf['qcolor'] = paf['qname'].apply(lambda x: x.split('_')[-1].split('-')[0])
        paf['qorient'] = paf['qname'].apply(lambda x: x.split('_')[-1].split('-')[-2])
        paf['qfrac'] = paf['block_length'] / paf['qlength']        
        paf['tcolor'] = paf['tname'].apply(lambda x: x.split('_')[-1].split('-')[0])
        paf['torient'] = paf['tname'].apply(lambda x: x.split('_')[-1].split('-')[-2])
        paf['tfrac'] = paf['block_length'] / paf['tlength']
        paf['score_orientation'] = paf.apply(score_orientation, axis=1, scoring=orientation_scoring)

        stats = col.Counter()
        # drop remaining (partial) self-mappings
        drop_rows = paf['qname'] == paf['tname']
        stats['partial_self_mappings'] = drop_rows.sum()
        paf = paf.loc[~drop_rows, :].copy()
        
        for (qcolor, tcolor), mappings in paf.groupby(['qcolor', 'tcolor']):
            key_prefix = f'{qcolor}_{tcolor}'
            if qcolor == tcolor:
                all_comb = pairs[qcolor]
                stats[f'{qcolor}_{tcolor}_pairs_total'] = all_comb
                stats[f'{qcolor}_{tcolor}_pairs_total_pct'] = round(mappings.shape[0] / all_comb * 100, 1)
            stats[key_prefix + '_mappings'] = mappings.shape[0]
            aln_fractions = np.sort(np.concatenate((mappings['qfrac'].values, mappings['tfrac'].values), dtype=np.float32))
            stats[key_prefix + '_lenfrac_mean'] = aln_fractions.mean()
            stats[key_prefix + '_lenfrac_stddev'] = aln_fractions.std()
            stats[key_prefix + '_lenfrac_50pct'] = aln_fractions[aln_fractions.size//2]
            stats[key_prefix + '_lenfrac_25pct'] = aln_fractions[int(aln_fractions.size * 0.25)]
            stats[key_prefix + '_lenfrac_80pct'] = aln_fractions[int(aln_fractions.size * 0.8)]
            stats[key_prefix + '_lenfrac_99pct'] = aln_fractions[int(aln_fractions.size * 0.99)]

        with open(output.table, 'w') as tsv:
            for key in sorted(stats.keys()):
                value = stats[key]
                _ = tsv.write(f'{key}\t{value}\n')
