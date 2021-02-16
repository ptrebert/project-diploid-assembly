
BREAK_CONFIG = {
    'bng_ro_overlap': 50,
    'reference': 'GRCh38_HGSVC2_noalt'
}

PS_ASSEMBLY_ALIGNMENT_PATH = os.path.join(os.getcwd(), 'output/alignments/contigs_to_reference/evaluation/phased_assemblies')
PS_ASSEMBLY_READ_COVERAGE = os.path.join(os.getcwd(), 'output/evaluation/hap_read_coverage')

BNG_UNIQ_CLUSTERS = 'Table_S13-4_Bionano_unique_clusters'

CHILD_SAMPLES = ['HG00733', 'HG00514', 'NA19240']

H64_ASSEMBLIES = sorted(
    set(
        [x.split('_map-to_')[0] for x in os.listdir(PS_ASSEMBLY_ALIGNMENT_PATH) if x.endswith('.bed') and x.split('_')[0] not in CHILD_SAMPLES]
    )
)

H70_ASSEMBLIES = sorted(
    set(
        [x.split('_map-to_')[0] for x in os.listdir(PS_ASSEMBLY_ALIGNMENT_PATH) if x.endswith('.bed')]
    )
)

REPEAT_CLASSES = [
    'SimpleRep',
    'Satellite',
    'LTR',
    'LINE',
    'SINE',
    'DNA',
    'LowComplex',
    'RNA',
    'Retroposon',
    'RC'
]

SIMPLE_REPEAT_CLASSES = [
    'homopolymer',
    'dinucleotide',
    'trinucleotide',
    'other'
]

REFERENCE_CHROMOSOMES = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']


rule split_repeatmasker_annotation:
    input:
        'references/downloads/{known_ref}.rmsk.tsv.gz'
    output:
        expand('references/annotation/{{known_ref}}.rmsk_{repclass}.bed',
                repclass=REPEAT_CLASSES),
        'references/annotation/{known_ref}.rmsk_all.bed'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 if attempt < 2 else 3072 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 if attempt < 2 else 3072 * attempt
    run:
        import os
        import pandas as pd

        out_columns = [
            '#chrom',
            'chromStart',
            'chromEnd',
            'name',
            'score',
            'strand',
            'repClass',
            'repFamily'
        ]
        df = pd.read_csv(input[0], sep='\t', header=0)
        
        ucsc_columns = dict([
            ('#swScore', 'score'),
            ('genoName', '#chrom'),
            ('genoStart', 'chromStart'),
            ('genoEnd', 'chromEnd'),
            ('repName', 'name')
        ])

        try:
            nonref_chrom = ~df['#chrom'].isin(REFERENCE_CHROMOSOMES)
        except KeyError:
            new_columns = [ucsc_columns.get(c, c) for c in df.columns]
            df.columns = new_columns
            nonref_chrom = ~df['#chrom'].isin(REFERENCE_CHROMOSOMES)

        uncertain_class = df['repClass'].str.contains('?', regex=False)
        uncertain_family = df['repFamily'].str.contains('?', regex=False)
        unknown_class = df['repClass'].str.contains('Unknown', regex=False)
        unspec_class = df['repClass'].str.contains('Unspecified', regex=False)

        drop_rows = (nonref_chrom | uncertain_class | uncertain_family | unknown_class | unspec_class)

        df = df.loc[~drop_rows, :].copy()

        df['out_class'] = df['repClass']
        df.loc[df['repClass'].str.contains('RNA', regex=False), 'out_class'] = 'RNA'
        df.loc[df['repClass'].str.contains('Simple_repeat', regex=False), 'out_class'] = 'SimpleRep'
        df.loc[df['repClass'].str.contains('Low_complexity', regex=False), 'out_class'] = 'LowComplex'

        out_classes = df['out_class'].unique().tolist()

        if not all([x in REPEAT_CLASSES for x in out_classes]):
            raise ValueError('Unhandled repeat class: {}'.format(sorted(out_classes)))
        
        out_path = os.path.dirname(output[0])
        for out_class in df['out_class'].unique():
            out_file = os.path.join(out_path, '{}.rmsk_{}.bed'.format(wildcards.known_ref, out_class))
            sub = df.loc[df['out_class'] == out_class, out_columns]
            sub.to_csv(out_file, sep='\t', header=True, index=False)
        
        out_path = output[-1]
        df[out_columns].to_csv(out_path, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule split_simple_repeat_annotation:
    input:
        'references/downloads/{known_ref}.simplerepeats.tsv.gz'
    output:
        expand('references/annotation/{{known_ref}}.smprep_{repclass}.bed',
                repclass=SIMPLE_REPEAT_CLASSES),
        'references/annotation/{known_ref}.smprep_all.bed'
    run:
        import os
        import pandas as pd

        out_columns = [
            '#chrom',
            'chromStart',
            'chromEnd',
            'name',
            'score',
            'strand',
            'motif',
            'copyNum',
            'perMatch',
            'perIndel'
        ]
        df = pd.read_csv(input[0], sep='\t', header=0)

        if 'motif' not in df.columns:
            df['motif'] = df['sequence']

        nonref_chrom = ~df['#chrom'].isin(REFERENCE_CHROMOSOMES)
        
        df = df.loc[~nonref_chrom, :].copy()

        df['length'] = df['chromEnd'] - df['chromStart']
        df['strand'] = '.'
        df['name'] = 'TRF_' + df['#chrom'] + '_' + df['chromStart'].astype(str) + '_' + df['length'].astype(str)
        df['out_class'] = 'other'
        df.loc[df['motif'].str.len() == 1, 'out_class'] = 'homopolymer'
        df.loc[df['motif'].str.len() == 2, 'out_class'] = 'dinucleotide'
        df.loc[df['motif'].str.len() == 3, 'out_class'] = 'trinucleotide'

        out_path = os.path.dirname(output[0])
        for out_class in df['out_class'].unique():
            out_file = os.path.join(out_path, '{}.smprep_{}.bed'.format(wildcards.known_ref, out_class))
            sub = df.loc[df['out_class'] == out_class, out_columns]
            sub.to_csv(out_file, sep='\t', header=True, index=False)

        out_file = output[-1]
        df[out_columns].to_csv(out_file, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


def map_segdup_color(frac_match):
    """
    Follows roughly Mitchell Vollger's color coding for T2T,
    is different from UCSC / hg38 color coding
    """

    if frac_match < 0.9:
        color_rgb, color_name = '147,112,219', 'purple'
    elif 0.9 <= frac_match < 0.98:
        color_rgb, color_name = '128,128,128', 'gray'
    elif 0.98 <= frac_match < 0.99:
        color_rgb, color_name = '204,204,0', 'yellow'
    elif 0.99 <= frac_match:
        color_rgb, color_name = '255,103,0', 'orange'
    else:
        raise ValueError('fracMatch value out of range: {}'.format(frac_match))

    return color_rgb, color_name


rule process_segdup_annotation:
    input:
        'references/downloads/{known_ref}.segdups.tsv.gz'
    output:
        'references/annotation/{known_ref}.segdups.bed'
    run:
        import os
        import pandas as pd

        out_columns = [
            '#chrom',
            'chromStart',
            'chromEnd',
            'name',
            'score',
            'strand',
            'color_rgb',
            'color_name'
        ]

        df = pd.read_csv(input[0], sep='\t', header=0)

        try:
            nonref_chrom = ~df['#chrom'].isin(REFERENCE_CHROMOSOMES)
        except KeyError:
            out_columns[0] = 'chrom'
            nonref_chrom = ~df['chrom'].isin(REFERENCE_CHROMOSOMES)

        df = df.loc[~nonref_chrom, :].copy()

        df['score'] = df['fracMatch'].apply(lambda x: int(min(round(x * 1000, 0), 1000)))
        assigned_colors = pd.DataFrame.from_records(
            df['fracMatch'].apply(map_segdup_color),
            index=df.index,
            columns=['color_rgb', 'color_name']
        )

        df = pd.concat([df, assigned_colors], axis=1, ignore_index=False)

        df.sort_values(['#chrom', 'chromStart', 'chromEnd'], axis=0, inplace=True)
        #sort_order = dict(('chr' + str(i), i) for i in range(1, 23))
        #sort_order['chrX'] = 23
        #sort_order['chrY'] = 24
        #df['sort_order'] = df['#chrom'].map(sort_order)

        #df.sort_values(['sort_order', 'chromStart', 'chromEnd'], axis=0, inplace=True)

        df[out_columns].to_csv(
            output[0],
            sep='\t',
            header=True,
            index=False
        )
    # END OF RUN BLOCK


rule extract_pav_no_bng_support:
    """
    Download source is the Zenodo upload (aka: submitted status)
    MD5: fe229fc1d66283211032019bd2cdff9b
    DOI: https://doi.org/10.5281/zenodo.4268827
    https://zenodo.org/record/4268828/files/variants_freeze3_sv_insdel.tsv.gz
    """
    input:
        'references/downloads/variants_freeze3_sv_insdel.tsv.gz'
    output:
        'references/annotation/GRCh38_HGSVC2_noalt.SVnobngV3.bed',
        'references/annotation/GRCh38_HGSVC2_noalt.DELnobngV3.bed',
        'references/annotation/GRCh38_HGSVC2_noalt.INSnobngV3.bed'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
    
        use_columns = [
            'ID',
            '#CHROM',
            'POS',
            'END',
            'SVTYPE',
            'SVLEN',
            'DISC_CLASS',
            'BIONANO'
        ]

        df = pd.read_csv(
            input[0],
            sep='\t',
            header=0,
            usecols=use_columns,
            encoding='ascii'
        )

        # select all variants with no BIONANO support
        df = df.loc[df['BIONANO'].isna(), :].copy()

        df['#chrom'] = df['#CHROM']
        df['start'] = df['POS'].astype('int64')
        df['end'] = df['start'] + df['SVLEN'].astype('int64')  # this affects INS
        df['name'] = df['ID']
        df['svtype'] = df['SVTYPE']
        df['svclass'] = df['DISC_CLASS']

        df.drop(use_columns, axis=1, inplace=True)

        df.sort_values(['#chrom', 'start', 'end'], ascending=True, inplace=True)

        dump_columns = [
            '#chrom',
            'start',
            'end',
            'name',
            'svtype',
            'svclass'
        ]

        selectors = [['DEL', 'INS'], ['DEL'], ['INS']]

        for selector, outfile in zip(selectors, output):
            sub = df.loc[df['svtype'].isin(selector), dump_columns]
            sub.to_csv(
                outfile,
                sep='\t',
                header=True,
                index=False
            )
    # END OF RUN BLOCK


###################################
###
### BELOW
### Process BNG unique clusters
###
###################################


BNG_TABLE_HEADER = [
    'name',
    'svtype',
    'chrom',
    'start',
    'end',
    'clusterSize',
    'num_samples',
    'InsCls',
    'DelCls',
    'NClustAtRegion',
    'InsMem',
    'DelMem',
    'NMemAtRegion',
    'SDOverlap',
    'GeneOverlap',
    'AFR',
    'AMR',
    'EAS',
    'EUR',
    'SAS'
]

BNG_REDUCED_HEADER = ['chrom', 'start', 'end', 'name', 'svtype', 'num_samples', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

rule preprocess_bng_uniq_table:
    input:
        'references/downloads/{}.tsv'.format(BNG_UNIQ_CLUSTERS)
    output:
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniqSV.sep.tsv',
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniqDEL.sep.tsv',
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniqINS.sep.tsv'
    log:
        'log/references/annotation/GRCh38_HGSVC2_noalt.BNGuniq.sep.log'
    run:
        import pandas as pd

        df = pd.read_csv(
            input[0],
            sep='\t',
            header=None,
            names=BNG_TABLE_HEADER,
            usecols=BNG_REDUCED_HEADER,
            comment='#'
        )

        df['name'] = 'BNG' + df['name'].astype(str)
        df['chrom'] = 'chr' + df['chrom']
        df['svtype'] = df['svtype'].map({'deletion': 'DEL', 'insertion': 'INS'})
        df = df[BNG_REDUCED_HEADER]
        
        # drop zero-length events - apparently, some type of error in the BNG data
        select_errors = (df['end'] - df['start']) < 10
        errors = df.loc[select_errors, :].copy()

        if not errors.empty:
            with open(log[0], 'w') as error_dump:
                _ = error_dump.write('## Zero-length clusters\n#')
                errors.to_csv(error_dump, sep='\t', header=True, index=False)

        df = df.loc[~select_errors, :].copy()

        with open(output[0], 'w') as dump:
            _ = dump.write('#')
            df.to_csv(
                dump,
                sep='\t',
                header=True,
                index=False
            )

        with open(output[1], 'w') as dump:
            _ = dump.write('#')
            df.loc[df['svtype'] == 'DEL', :].to_csv(
                dump,
                sep='\t',
                header=True,
                index=False
            )
        
        with open(output[2], 'w') as dump:
            _ = dump.write('#')
            df.loc[df['svtype'] == 'INS', :].to_csv(
                dump,
                sep='\t',
                header=True,
                index=False
            )
    # END OF RUN BLOCK


rule bng_clusters_self_overlap:
    """
    The merging strategy for the BNG unique clusters can be done in different ways:

    - do NOT merge if less than 80% size concordance, i.e. reciprocal overlap (comment by Joyce Lee)
    - do NOT merge if less than 50% size concordance, i.e. reciprocal overlap (this roughly follows PAV,
    although PAV variant merging strategy is more complex, see SI 9.2)

    """
    input:
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{svtype}.sep.tsv'
    output:
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}.ovlRO{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        svtype = '(SV|DEL|INS)'
    params:
        ro_overlap = round(BREAK_CONFIG['bng_ro_overlap'] / 100, 2)
    shell:
        'bedtools intersect -f {params.ro_overlap} -r -wo -a {input} -b {input} > {output}'


rule bng_clusters_merge:
    input:
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}.ovlRO{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap'])
    output:
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}-ro{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap']),
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}-ro{}.bed'.format(BREAK_CONFIG['bng_ro_overlap']),
        'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}-ro{}.stats'.format(BREAK_CONFIG['bng_ro_overlap'])
    wildcard_constraints:
        svtype = '(SV|DEL|INS)'
    run:
        import pandas as pd
        import itertools as itt
        import collections as col
        import numpy as np

        bng_columns = BNG_REDUCED_HEADER

        ovl_columns = ['_'.join(x) for x in itt.product(['a', 'b'], bng_columns)]
        ovl_columns.append('bp_overlap')
        int_columns = ['a_' + x for x in bng_columns[5:]]

        df = pd.read_csv(input[0], sep='\t', header=None, names=ovl_columns)
        self_overlaps = df.loc[df['a_name'] == df['b_name'], :]
        other_overlaps = df.loc[df['a_name'] != df['b_name'], :]

        clusters_to_merge = col.defaultdict(set)
        for pair, subset in other_overlaps.groupby(['a_name', 'b_name']):
            clusters_to_merge[pair[0]].add(pair[0])
            clusters_to_merge[pair[0]].add(pair[1])
            clusters_to_merge[pair[1]].add(pair[0])
            clusters_to_merge[pair[1]].add(pair[1])
        
        clusters_to_merge = col.deque(clusters_to_merge.values())
        clusters_to_merge.append(None)

        merged_clusters = []
        while 1:
            set_a = clusters_to_merge.popleft()
            if set_a is None:
                break
            while 1:
                set_b = clusters_to_merge.popleft()
                if set_b is None:
                    if len(clusters_to_merge) == 0:
                        merged_clusters.append(set_a)
                        clusters_to_merge.append(None)
                        break
                    else:
                        clusters_to_merge.append(None)
                        merged_clusters.append(set_a)
                        break
                ovl = set_a.intersection(set_b)
                if len(ovl) > 0:
                    set_a = set_a.union(set_b)
                else:
                    clusters_to_merge.append(set_b)
                continue    

        stats_clusters = col.Counter([len(x) for x in merged_clusters])

        merged_indices = []
        merged_rows = []
        for merge_set in merged_clusters:
            
            sub = self_overlaps.loc[self_overlaps['a_name'].isin(merge_set), :].copy()
            merged_indices.extend(sub.index.tolist())
            
            mrg_start = sub['a_start'].min()
            mrg_end = sub['a_end'].max()
            mrg_stats = sub.loc[:, int_columns].sum(axis=0)
            sub.sort_values('a_name', inplace=True, ascending=True)
            mrg_name = '|'.join(sub['a_name'].values)
            if sub['a_svtype'].nunique() == 1:
                mrg_type = sub.loc[sub.index[0], 'a_svtype']
            else:
                mrg_type = '|'.join(sub['a_svtype'].values)
            mrg_stats['name'] = mrg_name
            mrg_stats['svtype'] = mrg_type
            mrg_stats['start'] = mrg_start
            mrg_stats['end'] = mrg_end
            mrg_stats['chrom'] = sub.loc[sub.index[0], 'a_chrom']
            merged_rows.append(mrg_stats)

        if not len(merged_indices) == len(set(merged_indices)):
            raise ValueError('non-unique merge operation')
    
        merged_df = pd.DataFrame(merged_rows)
        merged_df.rename(dict((k, k.split('_', 1)[-1]) for k in int_columns), inplace=True, axis=1)

        stats_clean_samples = col.Counter(merged_df['num_samples'].values)

        select_not_merged = np.logical_not(self_overlaps.index.isin(merged_indices))
        unmerged_df = self_overlaps.loc[select_not_merged, :].copy()
        unmerged_df.rename(dict((k, k.split('_', 1)[-1]) for k in ovl_columns if k.startswith('a')), inplace=True, axis=1)

        final_df = pd.concat([unmerged_df, merged_df], axis=0, ignore_index=False)
        final_df.drop([c for c in final_df.columns if c.startswith('b')], axis=1, inplace=True)
        final_df.sort_values(['chrom', 'start', 'end'], inplace=True, ascending=True)

        with open(output[0], 'w') as dump:
            _ = dump.write('#')
            final_df[bng_columns].to_csv(
                dump,
                sep='\t',
                header=True,
                index=False
            )
        
        with open(output[1], 'w') as dump:
            _ = dump.write('#')
            final_df[['chrom', 'start', 'end', 'name']].to_csv(
                dump,
                sep='\t',
                header=True,
                index=False
            )
        
        with open(output[2], 'w') as dump:
            _ = dump.write('# BNG unique clusters merge info\n')
            _ = dump.write('num_clusters\t{}\n'.format(self_overlaps.shape[0]))
            _ = dump.write('num_overlap_pairs\t{}\n'.format(other_overlaps.shape[0] // 2))  # div by two because of reciprocality
            _ = dump.write('num_overlap_merged\t{}\n'.format(merged_df.shape[0]))
            _ = dump.write('num_clusters_final\t{}\n'.format(final_df.shape[0]))
            _ = dump.write('# statistic: merged clusters\n')
            _ = dump.write('# NUM members --- count\n')
            for num_members, count in stats_clusters.most_common():
                _ = dump.write('{}\t{}\n'.format(num_members, count))
            _ = dump.write('# statistic: clean samples per merged cluster\n')
            _ = dump.write('# NUM samples --- count\n')
            for num_samples, count in stats_clean_samples.most_common():
                _ = dump.write('{}\t{}\n'.format(num_samples, count))

    # END OF RUN BLOCK


#################################
###
### BELOW
### Process assembly breaks
###
#################################


rule merge_unfiltered_alignments:
    input:
        'output/alignments/contigs_to_reference/evaluation/phased_assemblies/{assembly}_map-to_{known_ref}.bed'
    output:
        'output/evaluation/break_analysis/ref_alignments/merged/{assembly}_map-to_{known_ref}.merged.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools merge -d 0 -c 4 -o distinct -i {input} > {output}'


rule complement_merged_alignments:
    input:
        'output/evaluation/break_analysis/ref_alignments/merged/{assembly}_map-to_{known_ref}.merged.bed',
        'references/assemblies/{known_ref}.fasta.fai'
    output:
        'output/evaluation/break_analysis/ref_alignments/complement/{assembly}_map-to_{known_ref}.complement.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools complement -L -i {input[0]} -g {input[1]} > {output}'


rule make_reference_windows:
    input:
        'references/assemblies/{known_ref}.fasta.fai'
    output:
        'references/assemblies/{known_ref}.win-100kbp.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools makewindows -g {input} -w 100000 > {output}'


rule compute_break_coverage_over_windows:
    input:
        'references/assemblies/{known_ref}.win-100kbp.bed',
        'output/evaluation/break_analysis/ref_alignments/complement/{assembly}_map-to_{known_ref}.complement.bed'
    output:
        'output/evaluation/break_analysis/ref_alignments/window_cov/{assembly}_map-to_{known_ref}.break-cov-100kbp.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools coverage -counts -a {input[0]} -b {input[1]} > {output}'


rule compute_break_windows:
    input:
        expand('output/evaluation/break_analysis/ref_alignments/window_cov/{assembly}_map-to_{{known_ref}}.break-cov-100kbp.bed',
                assembly=H64_ASSEMBLIES)
    output:
        'references/annotation/{known_ref}.H64.win-100kbp-breaks.bed',
        'references/annotation/{known_ref}.H64.win-100kbp-stable.bed',
    run:
        import pandas as pd

        hap_cov = None
        for cov_file in input:
            df = pd.read_csv(cov_file, sep='\t', header=None, index_col=[0, 1, 2], names=['#chrom', 'start', 'end', 'break_cov'])
            df['break_cov'] = df['break_cov'].clip(0, 1)  # just need to know if the haploid assembly is breaking
            if hap_cov is None:
                hap_cov = df
            else:
                hap_cov += df

        # consider windows with cumulative coverage > 1 as "break windows"
        hap_cov.loc[hap_cov['break_cov'] > 1, :].to_csv(output[0], sep='\t', index=True, header=True)

        # all other windows are considered stable
        hap_cov.loc[hap_cov['break_cov'] <= 1, :].to_csv(output[1], sep='\t', index=True, header=True)

    # END OF RUN BLOCK       


rule subtract_stable_regions_from_haploid_breaks:
    input:
        'output/evaluation/break_analysis/ref_alignments/complement/{assembly}_map-to_{known_ref}.complement.bed',
        'references/annotation/{known_ref}.H64.win-100kbp-stable.bed',
    output:
        'output/evaluation/break_analysis/ref_alignments/subtract/{assembly}_map-to_{known_ref}.minus-stable.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools subtract -a {input[0]} -b {input[1]} > {output}'


rule concat_haploid_breaks:
    """
    For whatever reason, bedtools merge does not eat the sort output here from stdin...
    """
    input:
        expand('output/evaluation/break_analysis/ref_alignments/subtract/{assembly}_map-to_{{known_ref}}.minus-stable.bed',
                assembly=H64_ASSEMBLIES)
    output:
        temp('references/annotation/{known_ref}.H64.raw-breaks.bed')
    shell:
        'cut -f 1,2,3 {input} | egrep "chr[0-9XY]+\s" | sort -k1,1 -k2,2n > {output}'


rule merge_haploid_breaks:
    input:
        'references/annotation/{known_ref}.H64.raw-breaks.bed'
    output:
        'references/annotation/{known_ref}.H64breaks.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools merge -i {input} | awk \'{{ OFS="\t" ; print ($1,$2,$3,$1"-"$2"-"$3-$2"-H64") }}\' > {output}'


rule intersect_breaks_alignments:
    input:
        'references/annotation/{known_ref}.H64breaks.bed',
        'output/evaluation/break_analysis/ref_alignments/complement/{assembly}_map-to_{known_ref}.complement.bed',
    output:
        'output/evaluation/break_analysis/break_coverages/assemblies/{known_ref}_cov_{assembly}.H64breaks.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools coverage -a {input[0]} -b {input[1]} > {output}'


rule intersect_region_set_annotations:
    input:
        'references/annotation/{known_ref}.{region_set}.bed',
        'references/annotation/{known_ref}.{annotation}.bed'
    output:
        'output/evaluation/break_analysis/break_coverages/annotation/{known_ref}_cov_{annotation}.{region_set}.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    wildcard_constraints:
        known_ref = BREAK_CONFIG['reference'],
        annotation = '[a-z0-9A-Z\_]+',
    shell:
        'bedtools coverage -a {input[0]} -b {input[1]} > {output}'


rule haploid_read_coverage_in_regions:
    input:
        'references/annotation/GRCh38_HGSVC2_noalt.{region_set}.bed',
        'output/evaluation/hap_read_coverage/{readset}_map-to_hg38_GCA_p13.{hap}.bigWig'
    output:
        'output/evaluation/break_analysis/hap_read_coverage/{region_set}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 if attempt < 2 else 49152 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 if attempt < 2 else 49152 * attempt
    shell:
        'bigWigAverageOverBed {input[1]} {input[0]} {output}'


rule map_scores_region_set:
    input:
        'references/annotation/{known_ref}.{region_set}.bed',
        'references/annotation/{known_ref}.{annotation}.bed'
    output:
        'output/evaluation/break_analysis/break_coverages/annotation/{known_ref}_score_{annotation}.{region_set}.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        annotation = 'segdups'
    shell:
        'bedtools map -prec 3 -c 5 -o mean -f 0.25 -a {input[0]} -b {input[1]} > {output}'


rule profile_region_set_sequence:
    input:
        'references/annotation/{known_ref}.{region_set}.bed',
        'references/assemblies/{known_ref}.fasta',
    output:
        'references/annotation/{known_ref}.{region_set}.nucprof.tsv',
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        known_ref = BREAK_CONFIG['reference']
    shell:
        'bedtools nuc -bed {input[0]} -fi {input[1]} > {output}'


TABLE_HEADERS = {
    'assembly_cov': ['chrom', 'start', 'end', 'name', 'num_overlaps', 'bp_coverage', 'total_length', 'pct_coverage'],
    'annotation_cov': ['chrom', 'start', 'end', 'name', 'num_overlaps', 'bp_coverage', 'total_length', 'pct_coverage'],
    'score_cov': ['chrom', 'start', 'end', 'name', 'mean_seqid'],
    'nuc_profile': [
        'chrom', 'start', 'end', 'name', 'pct_nucAT', 'pct_nucGC',
        'bp_nucA', 'bp_nucC', 'bp_nucG', 'bp_nucT', 'bp_nucN',
        'bp_nucX', 'bp_length'
        ],
    'read_cov': ['name', 'total_length', 'bp_coverage', 'total_sum', 'avg_depth', 'avg_nzdepth']
}

FEATURE_COLUMNS = ['bp', 'pct', 'avg']


def load_break_annotation(file_path, header, col_names, sample, platform, haplotype, name):

    import pandas as pd

    df = pd.read_csv(
        file_path,
        sep='\t',
        names=col_names,
        header=header,
        index_col=['chrom', 'start', 'end', 'name']
    )

    if 'mean_seqid' in df:
        df.loc[df['mean_seqid'] == '.', 'mean_seqid'] = '0'
        df['pct_seqid'] = df['mean_seqid'].astype(float) / 1000

    if 'bp_nucX' in df:
        # compute percentages for all nuc counts
        columns = df.columns
        for c in columns:
            if not c.startswith('bp_nuc'):
                continue
            pct_col = c.replace('bp_', 'pct_')
            df[pct_col] = (df[c] / df['bp_length']).round(3)

    col_index = pd.MultiIndex.from_tuples(
        [(sample, platform, haplotype, name, c.split('_')[0], c.split('_')[1]) for c in df.columns],
        names=['sample', 'platform', 'haplotype', 'name', 'unit', 'statistic']
    )

    df.columns = col_index
    select_columns = [u in FEATURE_COLUMNS for u in df.columns.get_level_values('unit')]
    df = df.loc[:, select_columns].copy()

    return df


def load_break_read_coverage(read_cov_file, break_regions):

    import pandas as pd

    components = os.path.basename(read_cov_file).split('_')
    sample = components[2]
    platform = {'pbsq2-clr': 'CLR', 'pbsq2-ccs': 'HiFi'}[components[4]]
    hap = read_cov_file.split('.')[-2].upper()
    if hap == 'UN':
        hap = 0
    elif hap == 'H1':
        hap = 1
    elif hap == 'H2':
        hap = 2
    else:
        raise ValueError('Cannot parse haplotype: {}'.format(hap))
   
    df = pd.read_csv(
        read_cov_file,
        sep='\t',
        names=TABLE_HEADERS['read_cov'],
        index_col=['name'],
        header=None,
    )

    df.drop(['total_length', 'total_sum'], axis=1, inplace=True)

    # this join adds chrom-start-end to the read coverage file data
    df = break_regions.join(df, how='outer')
   
    name_tuples = [(sample, platform, hap, 'read_cov', c.split('_')[0], c.split('_')[1]) for c in df.columns]

    col_index = pd.MultiIndex.from_tuples(
        name_tuples,
        names=['sample', 'platform', 'haplotype', 'name', 'unit', 'statistic']
    )
    df.columns = col_index
    
    df.drop([c for c in df.columns if 'BNG' in c], axis=1, inplace=True)

    return df


if 'GRCh38' not in BREAK_CONFIG['reference']:
    read_cov_regions_breaks = []
    read_cov_regions_bng = []
else:
    # HG00513_hgsvc_pbsq2-ccs_1000_map-to_hg38_GCA_p13.h2.bigWig
    # {region_set}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab

    out_path = 'output/evaluation/break_analysis/hap_read_coverage'

    read_cov_regions_breaks = []
    read_cov_regions_bng = []
    for bigwig in os.listdir(PS_ASSEMBLY_READ_COVERAGE):
        if '.h1.' not in bigwig:
            continue
        readset = bigwig.split('_map-to_')[0]
        for h in ['h1', 'h2', 'un']:
            break_file = '{region_set}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab'.format(
                **{
                    'region_set': 'H64breaks',
                    'readset': readset,
                    'hap': h
                }
            )
            read_cov_regions_breaks.append(os.path.join(out_path, break_file))

            bng_file = '{region_set}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab'.format(
                **{
                    'region_set': 'BNGuniq{{svtype}}-ro{}'.format(BREAK_CONFIG['bng_ro_overlap']),
                    'readset': readset,
                    'hap': h
                }
            )
            read_cov_regions_bng.append(os.path.join(out_path, bng_file))


rule merge_break_coverage_annotation:
    input:
        assm_cov = expand('output/evaluation/break_analysis/break_coverages/assemblies/{{known_ref}}_cov_{assembly}.H64breaks.tsv',
                            assembly=H70_ASSEMBLIES),
        annot_cov = expand('output/evaluation/break_analysis/break_coverages/annotation/{{known_ref}}_cov_{annotation}.H64breaks.tsv',
                            annotation=['rmsk_' + x for x in REPEAT_CLASSES] + ['smprep_' + i for i in SIMPLE_REPEAT_CLASSES] + ['segdups']),
        score_cov = ['output/evaluation/break_analysis/break_coverages/annotation/{known_ref}_score_segdups.H64breaks.tsv'],
        nuc_prof = 'references/annotation/{known_ref}.H64breaks.nucprof.tsv',
        read_cov = read_cov_regions_breaks,
        break_info = 'references/annotation/{known_ref}.H64breaks.bed'
    output:
        tsv = 'output/evaluation/break_analysis/tables/{known_ref}.H64breaks.tsv',
        hdf = 'output/evaluation/break_analysis/tables/{known_ref}.H64breaks.h5',
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        annotations = []

        breaks = pd.read_csv(
            input.break_info,
            sep='\t',
            names=['chrom', 'start', 'end', 'name'],
            header=None,
            index_col=['chrom', 'start', 'end', 'name']
        )

        for read_cov_file in input.read_cov:
            df = load_break_read_coverage(
                read_cov_file,
                breaks
            )
            annotations.append(df)
            
 
        for assm_cov_file in input.assm_cov:
            assembly = assm_cov_file.split('_cov_')[1]
            sample = assembly.split('_')[0]
            platform = 'CLR' if '-clr' in assembly else 'HiFi'
            hap = 10 if 'h1-un' in assembly else 20
            name = '{}_{}_H{}'.format(sample, platform, int(hap/10))
            df = load_break_annotation(
                assm_cov_file,
                None,
                TABLE_HEADERS['assembly_cov'],
                sample,
                platform,
                hap,
                'contig_cov'  #name
            )
            annotations.append(df)

        for annot_cov_file in input.annot_cov:
            annot = annot_cov_file.split('_cov_')[1]
            annot = annot.split('.')[0]
            ref_assm = wildcards.known_ref.split('_')[0]
            df = load_break_annotation(
                annot_cov_file,
                None,
                TABLE_HEADERS['annotation_cov'],
                ref_assm,
                'REF',
                0,
                annot
            )
            annotations.append(df)
        
        for score_cov_file in input.score_cov:
            annot = score_cov_file.split('_score_')[1]
            annot = annot.split('.')[0]
            ref_assm = wildcards.known_ref.split('_')[0]
            df = load_break_annotation(
                score_cov_file,
                None,
                TABLE_HEADERS['score_cov'],
                ref_assm,
                'REF',
                0,
                annot
            )
            annotations.append(df)

        df = load_break_annotation(
            input.nuc_prof,
            0,
            TABLE_HEADERS['nuc_profile'],
            ref_assm,
            'REF',
            0,
            'nucprof'
        )
        annotations.append(df)

        df = pd.concat(annotations, axis=1, ignore_index=False)
        df.sort_index(axis=0, inplace=True, ignore_index=False)
        df.sort_index(axis=1, inplace=True, ignore_index=False)
        
        df.to_csv(output.tsv, sep='\t', header=True, index=True)
        df.to_hdf(output.hdf, 'annotation', format='fixed')


def load_bng_annotation(file_path, header, col_names, name):

    import pandas as pd

    df = pd.read_csv(
            file_path,
            sep='\t',
            names=col_names,
            header=header,
            index_col=['chrom', 'start', 'end', 'name']
        )

    if 'mean_seqid' in df:
        df.loc[df['mean_seqid'] == '.', 'mean_seqid'] = '0'
        df['pct_seqid'] = df['mean_seqid'].astype(float) / 1000

    if 'bp_nucX' in df:
        # compute percentages for all nuc counts
        columns = df.columns
        for c in columns:
            if not c.startswith('bp_nuc'):
                continue
            pct_col = c.replace('bp_', 'pct_')
            df[pct_col] = (df[c] / df['bp_length']).round(3)

    try:
        col_index = pd.MultiIndex.from_tuples(
            [(name, c.split('_')[0], c.split('_')[1]) for c in df.columns],
            names=['name', 'unit', 'statistic']
        )
    except IndexError:
        print(name)
        print(col_names)
        print(df.head())
        print(df.columns)
        raise

    df.columns = col_index
    select_columns = [u in FEATURE_COLUMNS for u in df.columns.get_level_values('unit')]
    df = df.loc[:, select_columns].copy()

    return df


def load_bng_infos(file_path):

    import pandas as pd

    df = pd.read_csv(
        file_path,
        sep='\t',
        names=BNG_REDUCED_HEADER,  # ['chrom', 'start', 'end', 'clusterID', 'clusterSV', 'cleanNumSamples', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
        header=None,
        index_col=['chrom', 'start', 'end', 'name'],
        comment='#'
    )

    new_columns = {
        'name': 'category_clusterID',
        'svtype': 'category_SVtype',
        'num_samples': 'num_samples',
        'AFR': 'num_AFR',
        'AMR': 'num_AMR',
        'EAS': 'num_EAS',
        'EUR': 'num_EUR',
        'SAS': 'num_SAS'
    }

    df.rename(columns=new_columns, inplace=True)

    try:
        col_index = pd.MultiIndex.from_tuples(
            [('BNG', c.split('_')[0], c.split('_')[1]) for c in df.columns],
            names=['name', 'unit', 'statistic']
        )
    except IndexError:
        print(df.head())
        print(df.columns)
        raise

    df.columns = col_index
    
    return df


def load_bng_read_coverage(read_cov_file, bng_regions):

    import pandas as pd

    components = os.path.basename(read_cov_file).split('_')
    sample = components[2]
    platform = {'pbsq2-clr': 'CLR', 'pbsq2-ccs': 'HiFi'}[components[4]]
    hap = read_cov_file.split('.')[-2].upper()
    if hap == 'UN':
        hap = 'H0'

    data_column_name = '-'.join([sample, platform, hap])
    
    df = pd.read_csv(
        read_cov_file,
        sep='\t',
        names=TABLE_HEADERS['read_cov'],
        index_col=['name'],
        header=None,
    )

    df.drop(['total_length', 'total_sum'], axis=1, inplace=True)
    
    name_tuples = [(data_column_name, c.split('_')[0], c.split('_')[1]) for c in df.columns]

    col_index = pd.MultiIndex.from_tuples(
        name_tuples,
        names=['name', 'unit', 'statistic']
    )
    df.columns = col_index

    df = bng_regions.join(df, how='outer')
    
    df.drop([c for c in df.columns if 'BNG' in c], axis=1, inplace=True)

    return df


rule merge_bng_cluster_coverage_annotation:
    """
    Adapt overlap operation for reduced upset plot.
    Change from:
    ['segdups', 'H64breaks', '{svtype}nobngV3', 'INVv3', 'issues']
    to
    ['segdups', 'issues']
    """
    input:
        annot_cov = expand('output/evaluation/break_analysis/break_coverages/annotation/{{known_ref}}_cov_{annotation}.BNGuniq{{svtype}}-ro{ro}.tsv',
                            annotation=['rmsk_' + x for x in REPEAT_CLASSES] + ['smprep_' + i for i in SIMPLE_REPEAT_CLASSES] + ['segdups', 'issues'],
                            ro=BREAK_CONFIG['bng_ro_overlap']
                    ),
        score_cov = ['output/evaluation/break_analysis/break_coverages/annotation/{{known_ref}}_score_segdups.BNGuniq{{svtype}}-ro{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap'])],
        nuc_prof = 'references/annotation/{{known_ref}}.BNGuniq{{svtype}}-ro{}.nucprof.tsv'.format(BREAK_CONFIG['bng_ro_overlap']),
        bng_annotation = 'references/annotation/GRCh38_HGSVC2_noalt.BNGuniq{{svtype}}-ro{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap']),
        read_cov = []  #read_cov_regions_bng
    output:
        tsv = 'output/evaluation/break_analysis/tables/{{known_ref}}.BNGuniq{{svtype}}-ro{}.tsv'.format(BREAK_CONFIG['bng_ro_overlap']),
        hdf = 'output/evaluation/break_analysis/tables/{{known_ref}}.BNGuniq{{svtype}}-ro{}.h5'.format(BREAK_CONFIG['bng_ro_overlap']),
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        annotations = []
        
        bng_regions = load_bng_infos(input.bng_annotation)
        annotations.append(bng_regions)

        for read_cov_file in input.read_cov:
            df = load_bng_read_coverage(
                read_cov_file,
                bng_regions
            )
            annotations.append(df)

        for annot_cov_file in input.annot_cov:
            annot = annot_cov_file.split('_cov_')[1]
            annot = annot.split('.')[0]
            if annot == 'H64':
                annot = 'H64breaks'
            ref_assm = wildcards.known_ref.split('_')[0]
            df = load_bng_annotation(
                annot_cov_file,
                None,
                TABLE_HEADERS['annotation_cov'],
                annot
            )
            annotations.append(df)
        
        for score_cov_file in input.score_cov:
            annot = score_cov_file.split('_score_')[1]
            annot = annot.split('.')[0]
            ref_assm = wildcards.known_ref.split('_')[0]
            df = load_bng_annotation(
                score_cov_file,
                None,
                TABLE_HEADERS['score_cov'],
                annot
            )
            annotations.append(df)

        df = load_bng_annotation(
            input.nuc_prof,
            0,
            TABLE_HEADERS['nuc_profile'],
            'nucprof'
        )
        annotations.append(df)

        df = pd.concat(annotations, axis=1, ignore_index=False)
        df.sort_index(axis=0, inplace=True, ignore_index=False)
        df.sort_index(axis=1, inplace=True, ignore_index=False)
        
        df.to_csv(output.tsv, sep='\t', header=True, index=True)
        df.to_hdf(output.hdf, 'annotation', format='fixed')

