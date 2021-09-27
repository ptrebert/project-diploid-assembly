
localrules: dump_contig_name, run_all

PGAS_PATH = '/gpfs/project/projects/medbioinf/projects/hgsvc/2021_pgas/run_folder'
ONTQC_PATH = '/gpfs/project/projects/medbioinf/projects/hgsvc/ontqc/run_folder'
REFERENCE = 'T2Tv11_hg002Yv2_chm13'

rule run_all:
    input:
        'output/read_ctg_align/NA18989_HIFIEC-SUB_MAP-TO_chimeric-contigs.mmap.paf.gz',
        'output/read_ctg_align/HG02666_HIFIEC-SUB_MAP-TO_chimeric-contigs.mmap.paf.gz',
        'output/read_ctg_align/NA18989_HIFIEC_MAP-TO_chimeric-contigs.mmap.paf.gz',
        'output/read_ctg_align/HG02666_HIFIEC_MAP-TO_chimeric-contigs.mmap.paf.gz',
        'output/read_ctg_align/NA18989_ONTUL_MAP-TO_chimeric-contigs.mmap.paf.gz',
        'output/read_ctg_align/HG02666_ONTUL_MAP-TO_chimeric-contigs.mmap.paf.gz',
        expand('output/repmask/{sample}.repmask.tsv', sample=['NA18989', 'HG02666']),
        expand(
            'output/contig_ref_align/{sample}_CTG_MAP-TO_{reference}.mmap.paf.gz',
            sample=['NA18989', 'HG02666'],
            reference=[REFERENCE]
        ),
        expand(
            'output/cache/{sample}_{readset}.h5',
            sample=['HG02666', 'NA18989'],
            readset=['ONTUL', 'HIFIEC', 'HIFIEC-SUB']
        ),
        expand(
            'output/cache/{sample}.repmask.h5',
            sample=['HG02666', 'NA18989']
        )



rule dump_contig_name:
    input:
        asm_error = ancient(os.path.join(
            PGAS_PATH,
            'output/reference_assembly/clustered/temp/saarclust/results',
            '{sample}_hgsvc_pbsq2-ccs_1000_nhr-hifiasm/{sample}_hgsvc_ilnxs-80pe_sseq/clustered_assembly',
            'asmErrorsReport_2e+05bp_dynamic.tsv'
        )),
        ctg_order = ancient(os.path.join(
            PGAS_PATH,
            'output/reference_assembly/clustered/temp/saarclust/results',
            '{sample}_hgsvc_pbsq2-ccs_1000_nhr-hifiasm/{sample}_hgsvc_ilnxs-80pe_sseq/clustered_assembly',
            'ordered\&oriented_2e+05bp_dynamic.tsv'
        )),
    output:
        ctg_names = 'output/{sample}.chimeric-contigs.names.txt',
        ctg_info = 'output/{sample}.chimeric-contigs.info.tsv',
    run:
        import pandas as pd

        error_header = [
            'ctg_name',
            'start',
            'end',
            'width',
            'strand',
            'error_type',
            'misasm_bases'
        ]
        df_error = pd.read_csv(input.asm_error, sep='\t', header=0, names=error_header, usecols=[0,5,6])
        df_error = df_error.loc[df_error['error_type'] == 'chimerism', :].copy()

        df_error = df_error.loc[df_error['misasm_bases'] >= int(1e6), :].copy()
        assert not df_error.empty

        order_header = [
            'ctg_name',
            'start',
            'end',
            'length',
            'strand',
            'direction',
            'order_num',
            'cluster_id'
        ]

        df_order = pd.read_csv(input.ctg_order, sep='\t', header=0, names=order_header, usecols=[0,1,2,3,5,6,7])
        df_order['direction'] = (df_order['direction'].replace({'dir': 1, 'revcomp': -1})).astype('int8')
        df_order = df_order.loc[df_order['ctg_name'].isin(df_error['ctg_name']), :].copy()

        df_order = df_order.merge(df_error, left_on='ctg_name', right_on='ctg_name')

        contig_names = df_order['ctg_name'].unique().tolist()
        with open(output.ctg_names, 'w') as dump:
            _ = dump.write('\n'.join(contig_names) + '\n')

        df_order.to_csv(output.ctg_info, sep='\t', header=True, index=False)

    # END OF RUN BLOCK


rule extract_contig_sequences:
    input:
        ctg_names = 'output/{sample}.chimeric-contigs.names.txt',
        assembly = ancient(os.path.join(PGAS_PATH, 'output/reference_assembly/non-hap-res', '{sample}_hgsvc_pbsq2-ccs_1000_nhr-hifiasm.fasta'))
    output:
        'output/{sample}.chimeric-contigs.fasta'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    shell:
        'seqtk subseq {input.assembly} {input.ctg_names} > {output}'


rule generate_fasta_index:
    input:
        fasta = 'output/{sample}.chimeric-contigs.fasta'
    output:
        'output/{sample}.chimeric-contigs.fasta.fai'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input.fasta}'


rule extract_contig_read_names:
    input:
        ctg_names = 'output/{sample}.chimeric-contigs.names.txt',
        graph = ancient(
            os.path.join(
                PGAS_PATH,
                'output/reference_assembly/non-hap-res/layout/hifiasm',
                '{sample}_hgsvc_pbsq2-ccs_1000/{sample}_hgsvc_pbsq2-ccs_1000.p_ctg.noseq.gfa'
            )
        ),
        reads = ancient(
            os.path.join(
                PGAS_PATH,
                'output/reference_assembly/non-hap-res/layout/hifiasm',
                '{sample}_hgsvc_pbsq2-ccs_1000.ec-reads.fasta.gz'
            )
        )
    output:
        read_info = 'output/{sample}.chimeric-contigs.reads.tsv',
        renamed_reads = 'output/read_subsets/{sample}.HIFIEC.fasta.gz'
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    run:
        import pandas as pd
        import gzip as gzip
        import io as io
        import collections as col
        import sys

        contig_names = set(cn.strip() for cn in open(input.ctg_names, 'r').readlines())

        records = []
        name_lut = dict()
        dir_map = {
            '+': 'forward',
            '-': 'reverse'
        }
        read_mult = col.Counter()
        with open(input.graph, 'r') as graph:
            for line in graph:
                if not line.startswith('A'):
                    continue
                _, contig, ctg_start, read_dir, read_name, read_start, read_end, _, _ = line.strip().split()
                if contig not in contig_names:
                    continue
                new_name = f'{read_name}|CTGNAME_{contig}|CTGSTART_{ctg_start}|RDDIR_{dir_map[read_dir]}|RDSTART_{read_start}|RDEND_{read_end}'
                records.append((contig, ctg_start, read_dir, read_name, read_start, read_end, new_name))
                # NB: it indeed happens that a read is used more than once for assembly,
                # just keep last record as new name for now
                #assert read_name not in name_lut, f'exists: {read_name} to {name_lut[read_name]} / collision: {new_name}'
                read_mult[read_name] += 1
                name_lut[read_name] = new_name

        df = pd.DataFrame.from_records(records, columns=['contig', 'ctg_start', 'read_dir', 'read_name', 'read_start', 'read_end', 'read_new_name'])

        found = 0
        limit = len(name_lut)
        sys.stdout.write(f'Total reads to extract: {limit}\n')
        out_buffer = io.StringIO()
        read_lengths = {}
        with gzip.open(input.reads, 'rt') as fasta:
            while found < limit:
                read_name = fasta.readline().strip().strip('>')
                if read_name not in name_lut:
                    _ = fasta.readline()
                    continue
                new_name = name_lut[read_name]
                read_sequence = fasta.readline().strip()
                read_lengths[read_name] = len(read_sequence)
                out_buffer.write(f'>{new_name}\n{read_sequence}\n')
                found += 1
                if found % 10000 == 0:
                    sys.stdout.write(f'Found {found} reads\n')
        assert found == limit
        sys.stdout.write(f'Dumping buffer\n')
        with gzip.open(output.renamed_reads, 'wt') as fasta:
            _ = fasta.write(out_buffer.getvalue())

        sys.stdout.write(f'Dumping metadata\n')
        df['read_length'] = (df['read_name'].replace(read_lengths)).astype('int32')
        df['read_multiplicity'] = (df['read_name'].replace(read_mult)).astype('int8')
        df.to_csv(output.read_info, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule align_chimeric_contigs:
    input:
        contigs = 'output/{sample}.chimeric-contigs.fasta',
        reference = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        'output/contig_ref_align/{sample}_CTG_MAP-TO_{reference}.mmap.paf.gz'
    benchmark:
        'rsrc/output/contig_ref_align/{sample}_CTG_MAP-TO_{reference}.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 24
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt,
    shell:
        'minimap2 -t 22 -x asm20 --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.contigs} | pigz --best -p 4 > {output}'


rule align_ont_reads:
    """
    """
    input:
        reads = ancient(os.path.join(ONTQC_PATH, 'input/ONTUL', '{sample}_ONTUL_guppy-5.0.11-sup-prom.fasta.gz')),
        reference = 'output/{sample}.chimeric-contigs.fasta'
    output:
        paf = 'output/read_ctg_align/{sample}_ONTUL_MAP-TO_chimeric-contigs.mmap.paf.gz',
    benchmark:
        'rsrc/output/read_ctg_align/{sample}_ONTUL_MAP-TO_chimeric-contigs.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 24
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt,
    shell:
        'minimap2 -t 22 -x map-ont -k17 --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p 4 > {output}'


rule align_hifiec_reads:
    """
    """
    input:
        reads = ancient(ancient(
            os.path.join(
                PGAS_PATH,
                'output/reference_assembly/non-hap-res/layout/hifiasm',
                '{sample}_hgsvc_pbsq2-ccs_1000.ec-reads.fasta.gz'
            )
        )),
        reference = 'output/{sample}.chimeric-contigs.fasta'
    output:
        paf = 'output/read_ctg_align/{sample}_HIFIEC_MAP-TO_chimeric-contigs.mmap.paf.gz',
    benchmark:
        'rsrc/output/read_ctg_align/{sample}_HIFIEC_MAP-TO_chimeric-contigs.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 24
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt,
    shell:
        'minimap2 -t 22 -x map-hifi --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p 4 > {output}'


rule align_hifiec_subset_reads:
    """
    """
    input:
        reads = ancient('output/read_subsets/{sample}.HIFIEC.fasta.gz'),
        reference = 'output/{sample}.chimeric-contigs.fasta'
    output:
        paf = 'output/read_ctg_align/{sample}_HIFIEC-SUB_MAP-TO_chimeric-contigs.mmap.paf.gz',
    benchmark:
        'rsrc/output/read_ctg_align/{sample}_HIFIEC-SUB_MAP-TO_chimeric-contigs.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 24
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'minimap2 -t 22 -x map-hifi --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p 4 > {output}'


rule mask_repeats_in_chimeric_contigs:
    input:
        contigs = 'output/{sample}.chimeric-contigs.fasta'
    output:
        repmsk_dir = directory('output/repmask/{sample}')
    threads: 36
    resources:
        mem_total_mb = lambda wildcards, attempt: 524288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt
    shell:
        'module load RepeatMasker && rm -rfd {output.repmsk_dir} && mkdir -p {output.repmsk_dir} && '
        'RepeatMasker -species human -no_is -pa {threads} -s -dir {output.repmsk_dir} {input.contigs}'


rule normalize_repmask_output:
    input:
        rm_folder = 'output/repmask/{sample}'
    output:
        'output/repmask/{sample}.repmask.tsv'
    params:
        input_file = lambda wildcards, input: os.path.join(input.rm_folder, f'{wildcards.sample}.chimeric-contigs.fasta.out')
    shell:
        'tail -n +4 {params.input_file} | tr -s " " "\\t" | sed "s/^\\t//" > {output}'


def load_sequence_sizes(file_path):
    seq_sizes = {}
    with open(file_path, 'r') as fa_idx:
        for line in fa_idx:
            c, s = line.split()[:2]
            seq_sizes[c] = int(s)
    return seq_sizes


RM_HEADER = [
    'SW_score',
    'sub_pct',
    'del_pct',
    'ins_pct',
    'seq_name',
    'start',
    'end',
    'remain',
    'orientation',
    'repeat_name',
    'repeat_class',
    'column_end1',
    'column_end2',
    'column_end3',
    'column_end4',
    'column_end5'
]

RM_USE = [
    'seq_name',
    'start',
    'end',
    'repeat_class'
]


rule cache_repmask_output:
    input:
        table = 'output/repmask/{sample}.repmask.tsv',
        faidx = 'output/{sample}.chimeric-contigs.fasta.fai'
    output:
        'output/cache/{sample}.repmask.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import pandas as pd
        import numpy as np
        seq_sizes = load_sequence_sizes(input.faidx)

        with pd.HDFStore(output[0], 'w', complevel=9, complib='blosc'):
            pass

        rm = pd.read_csv(input.table, sep='\t', header=None, index_col=None, names=RM_HEADER, usecols=RM_USE)
        rm['is_all'] = 1
        rm['is_simple'] = rm['repeat_class'].apply(lambda x: 1 if 'Simple_repeat' in x else 0)
        rm['is_lowcomplex'] = rm['repeat_class'].apply(lambda x: 1 if 'Low_complexity' in x else 0)
        rm['is_satellite'] = rm['repeat_class'].apply(lambda x: 1 if 'Satellite' in x else 0)

        for sequence in rm['seq_name'].unique():
            for rep_select in ['is_all', 'is_simple', 'is_lowcomplex', 'is_satellite']:
                selector = np.logical_and(rm['seq_name'] == sequence, rm[rep_select] == 1)
                subset = rm.loc[selector, :]
                bases = np.zeros(seq_sizes[sequence], dtype=np.int8)
                for row in subset.itertuples(index=False):
                    bases[row[1]:row[2]+1] = 1
                with pd.HDFStore(output[0], 'a', complevel=9, complib='blosc') as hdf:
                    hdf.put(f'{wildcards.sample}/{sequence}/repmask/{rep_select}', pd.Series(bases), format='fixed')
    # END OF RUN BLOCK


PAF_HEADER = [
    'read_name',
    'read_length',
    'read_aln_start',
    'read_aln_end',
    'orientation',
    'ref_name',
    'ref_length',
    'ref_aln_start',
    'ref_aln_end',
    'residue_matches',
    'block_length',
    'mapq',
    'aln_type',
    'tag_cm',
    'tag_s1',
    'tag_s2',
    'divergence',
    'tag_rl'
]

PAF_USE = [
    'read_length',
    'read_aln_start',
    'read_aln_end',
    'ref_name',
    'ref_aln_start',  # idx 4
    'ref_aln_end',  # idx 5
    'residue_matches',
    'block_length',
    'mapq',  # idx 8
    'divergence',  # idx 9
]


rule cache_read_coverage_output:
    input:
        paf = 'output/read_ctg_align/{sample}_{readset}_MAP-TO_chimeric-contigs.mmap.paf.gz',
        faidx = 'output/{sample}.chimeric-contigs.fasta.fai'
    output:
        'output/cache/{sample}_{readset}.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import pandas as pd
        import numpy as np
        seq_sizes = load_sequence_sizes(input.faidx)

        aln = pd.read_csv(input.paf, sep='\t', header=None, names=PAF_HEADER, usecols=PAF_USE)
        aln['divergence'] = aln['divergence'].apply(lambda x: float(x.split(':')[-1]))

        idx_start = 4
        idx_end = 5
        idx_mapq = 8
        idx_div = 9

        with pd.HDFStore(output[0], 'w', complevel=9, complib='blosc'):
            pass

        for sequence in aln['ref_name'].unique():
            coverage = np.zeros(seq_sizes[sequence], dtype=np.int32)
            mapq = np.zeros(seq_sizes[sequence], dtype=np.float32)
            divergence = np.zeros(seq_sizes[sequence], dtype=np.float32)
            subset = aln.loc[aln['ref_name'] == sequence, :]
            for row in subset.itertuples(index=False):
                coverage[row[idx_start]:row[idx_end]] += 1
                mapq[row[idx_start]:row[idx_end]] += row[idx_mapq]
                divergence[row[idx_start]:row[idx_end]] += row[idx_div]
            nz_values = coverage > 0
            mapq[nz_values] /= coverage[nz_values]
            divergence[nz_values] /= coverage[nz_values]
            with pd.HDFStore(output[0], 'a', complevel=9, complib='blosc') as hdf:
                hdf.put(f'{wildcards.sample}/read_cov/{wildcards.readset.replace("-", "")}', pd.Series(coverage), format='fixed')
                hdf.put(f'{wildcards.sample}/aln_mapq/{wildcards.readset.replace("-", "")}', pd.Series(mapq), format='fixed')
                hdf.put(f'{wildcards.sample}/aln_div/{wildcards.readset.replace("-", "")}', pd.Series(divergence), format='fixed')
            
        if wildcards.readset == 'ONTUL':
            aln = aln.loc[aln['read_length'] >= 1e5, :].copy()
            aln['read_aln_fraction'] = (aln['read_aln_end'] - aln['read_aln_end']) / aln['read_length']
            aln = aln.loc[aln['read_aln_fraction'] >= 0.8, :].copy()
            aln.drop('read_aln_fraction', axis=1, inplace=True)
            for sequence in aln['ref_name'].unique():
                subset = aln.loc[aln['ref_name'] == sequence, :]
                coverage = np.zeros(seq_sizes[sequence], dtype=np.int32)
                mapq = np.zeros(seq_sizes[sequence], dtype=np.float32)
                divergence = np.zeros(seq_sizes[sequence], dtype=np.float32)
                for row in subset.itertuples(index=False):
                    coverage[row[idx_start]:row[idx_end]] += 1
                    mapq[row[idx_start]:row[idx_end]] += row[idx_mapq]
                    divergence[row[idx_start]:row[idx_end]] += row[idx_div]
                nz_values = coverage > 0
                mapq[nz_values] /= coverage[nz_values]
                divergence[nz_values] /= coverage[nz_values]
                with pd.HDFStore(output[0], 'a', complevel=9, complib='blosc') as hdf:
                    hdf.put(f'{wildcards.sample}/read_cov_UL', pd.Series(coverage), format='fixed')
                    hdf.put(f'{wildcards.sample}/aln_mapq_UL', pd.Series(mapq), format='fixed')
                    hdf.put(f'{wildcards.sample}/aln_div_UL', pd.Series(divergence), format='fixed')
        # END OF RUN BLOCK
