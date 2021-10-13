
rule dump_tigs_to_fasta:
    input:
        tigs = lambda wildcards: str(SAMPLE_INFOS[wildcards.sample][wildcards.tigs]),
    output:
        fasta = 'output/tigs/{sample_info}_{sample}.{tigs}.fasta'
    benchmark:
        'rsrc/output/tigs/{sample_info}_{sample}.{tigs}.dump.rsrc'
    wildcard_constraints:
        tigs = '(TIGPRI|TIGALT|TIGRAW)',
        sample = CONSTRAINT_REGULAR_SAMPLES
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        header_prefix = f'>{wildcards.sample_info}_{wildcards.sample}.'
        with open(output.fasta, 'w') as fasta:
            with open(input.tigs, 'r') as gfa:
                for ln, line in enumerate(gfa, start=1):
                    if not line.startswith('S'):
                        continue
                    try:
                        parts = line.split()
                        node_id = parts[1].strip()
                        sequence = parts[2].strip()
                    except (ValueError, IndexError) as err:
                        if not line.strip():
                            continue
                        raise ValueError(f'LN {ln} : FILE {gfa_file}: {str(err)}')
                    _ = fasta.write(f'{header_prefix}{node_id}\n{sequence}\n')
    # END OF RUN BLOCK


rule align_tigs_to_reference:
    input:
        fa_tigs = 'output/tigs/{sample_info}_{sample}.{tigs}.fasta',
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        bam = 'output/tig_aln/split/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.psort.bam',
        bai = 'output/tig_aln/split/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.psort.bam.bai',
        temp_dir = temp(
            directory('temp/minimap/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}')
        )
    log:
        mm = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.mm.log',
        view = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.view.log',
        sort = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.sort.log',
    benchmark:
        'rsrc/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.rsrc',
    wildcard_constraints:
        tigs = '(TIGPRI|TIGALT|TIGRAW)'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: (config['num_cpu_medium'] + config['num_cpu_low'])
    resources:
        mem_total_mb = lambda wildcards, attempt: 36864 + 36864 * attempt,
        runtime_hrs = lambda wildcards, attempt: 16 * attempt,
        mem_sort_mb = 4096
    params:
        align_threads = config['num_cpu_medium'],
        sort_threads = config['num_cpu_low'],
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: wildcards.tigs,
        discard_flag = 260  # unmapped or not primary
    shell:
        'mkdir -p {output.temp_dir} && '
        'minimap2 -t {params.align_threads} '
            '--secondary=no --eqx -L -Y -ax asm20 -m 10000 '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}\\tCN:HHU" '
            '{input.fa_ref} {input.fa_tigs} 2> {log.mm} | '
        'samtools view -u -F {params.discard_flag} /dev/stdin 2> {log.view} | '
        'samtools sort -@ {params.sort_threads} -m {resources.mem_sort_mb}M -T {output.temp_dir}/tmp_part -O BAM -l 6 /dev/stdin > {output.bam} 2> {log.sort} '
        ' && '
        'samtools index {output.bam}'


rule merge_tigs_to_reference_bams:
    input:
        bams = expand(
            'output/tig_aln/split/{{sample_info}}_{{sample}}_{tigs}_MAP-TO_{{reference}}.psort.bam',
            tigs=['TIGPRI', 'TIGALT']
        ),
        bais = expand(
            'output/tig_aln/split/{{sample_info}}_{{sample}}_{tigs}_MAP-TO_{{reference}}.psort.bam.bai',
            tigs=['TIGPRI', 'TIGALT']
        )
    output:
        bam = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam',
        bai = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam.bai',
    wildcard_constraints:
        sample = CONSTRAINT_REGULAR_SAMPLES
    log:
        'log/output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.merge.log'
    benchmark:
        'rsrc/output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.merge.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools merge -l 6 -O BAM --no-PG -@ {threads} {output.bam} {input.bams} 2> {log} '
            ' && '
            'samtools index {output.bam}'


rule dump_tigs_to_reference_to_bed:
    input:
        bam = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam',
        bai = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam.bai',
    output:
        bed = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bed',
    wildcard_constraints:
        sample = CONSTRAINT_REGULAR_SAMPLES
    benchmark:
        'rsrc/output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.dumpbed.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'bedtools bamtobed -i {input.bam} > {output.bed}'


def load_sequence_names(file_path):
    seq_names = set()
    with open(file_path, 'r') as fasta:
        for line in fasta:
            if not line.startswith('>'):
                continue
            seq_name = line.strip()[1:]
            seq_names.add(seq_name)
    return seq_names


rule classify_tigs:
    input:
        tigs_pri = 'output/tigs/{sample_info}_{sample}.TIGPRI.fasta',
        tigs_alt = 'output/tigs/{sample_info}_{sample}.TIGALT.fasta',
        bed = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bed',
    output:
        multiext(
            'output/tig_aln/chrom_groups/{sample_info}_{sample}_MAP-TO_{reference}.tigs',
            '.AM.pass.txt', '.AM.fail.txt',
            '.XY.pass.txt', '.XY.fail.txt',
            '.AMXY.fail.txt', '.UN.fail.txt'
        )
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: max(0, attempt - 1)
    run:
        import pandas as pd
        import collections as col
        bed_columns = ['chrom', 'start', 'end', 'tig', 'mapq', 'orientation']
        df = pd.read_csv(input.bed, sep='\t', header=None, names=bed_columns)
        auto_mito = [f'chr{i}' for i in range(1, 23)] + ['chrM']
        df['karyo_group'] = df['chrom'].apply(lambda x: 'AM' if x in auto_mito else 'XY')
        all_tigs = set().union(load_sequence_names(input.tigs_pri), load_sequence_names(input.tigs_alt))

        # init output to ensure that all output files are always created
        groups = {
            ('AM', 'pass'): [],
            ('AM', 'fail'): [],
            ('XY', 'pass'): [],
            ('XY', 'fail'): [],
            ('AMXY', 'fail'): [],
            ('UN', 'fail'): []
        }

        aligned_tigs = set()
        for contig, alignments in df.groupby('tig'):
            short_tig_name = contig.split('.')[1]
            ext_tig_name = contig + '.'
            if alignments['chrom'].nunique() == 1:
                chrom = alignments['chrom'].values[0]
                ext_tig_name += chrom
                group_key = (alignments['karyo_group'].values[0], 'pass')
            elif alignments['karyo_group'].nunique() == 1:
                chrom = 'chr' + alignments['karyo_group'].values[0]
                ext_tig_name += chrom
                group_key = (alignments['karyo_group'].values[0], 'fail')
            elif alignments['chrom'].nunique() > 1 and alignments['karyo_group'].nunique() > 1:
                chrom = 'chrNN'
                ext_tig_name += chrom
                group_key = ('AMXY', 'fail')
            else:
                raise ValueError(f'Cannot process alignments: {alignments}')
            groups[group_key].append(
                (
                    contig,
                    short_tig_name,
                    ext_tig_name
                )
            )
            aligned_tigs.add(contig)
        
        for contig in (all_tigs - aligned_tigs):
            groups[('UN', 'fail')].append(
                (
                    contig,
                    contig.split('.')[1],
                    contig + '.chrUN'
                )
            )
        
        base_out = output[0].rsplit('.', 3)[0]
        for (karyo_group, label), tig_names in groups.items():
            out_file = base_out + f'.{karyo_group}.{label}.txt'
            with open(out_file, 'w') as dump:
                _ = dump.write('\n'.join(['\t'.join(t) for t in sorted(tig_names)]) + '\n')
    # END OF RUN BLOCK


rule generate_gonosomal_reference:
    input:
        tigs_pri = 'output/tigs/{sample_long}.TIGPRI.fasta',
        tigs_alt = 'output/tigs/{sample_long}.TIGALT.fasta',
        tig_groups = multiext(
            'output/tig_aln/chrom_groups/{sample_long}_MAP-TO_{reference}.tigs',
            '.XY.pass.txt', '.XY.fail.txt',
            '.AMXY.fail.txt', '.UN.fail.txt'
        )
    output:
        fasta = 'output/gonosomal_reference/fasta/{sample_long}.{reference}.AMXYUN.tigs.fasta',
        gfa = 'output/gonosomal_reference/graph/{sample_long}.{reference}.AMXYUN.tigs.gfa',
        stats = 'output/gonosomal_reference/{sample_long}.{reference}.AMXYUN.tigs.stats.tsv',
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import pandas as pd
        import io as io
        import collections as col

        tig_table = []
        for tig_group in input.tig_groups:
            df = pd.read_csv(tig_group, sep='\t', names=['tig_long', 'tig', 'tig_ext_group'])
            if df.empty:
                continue
            tig_table.append(df)
        tig_table = pd.concat(tig_table, axis=0, ignore_index=False)

        fasta_buffer = io.StringIO()
        graph_buffer = io.StringIO('H\tVN:Z:1\n')
        length_stats = col.Counter()
        buffer = False
        seq_name = None
        for fasta_file in [input.tigs_pri, input.tigs_alt]:
            with open(fasta_file, 'r') as fasta:
                for line in fasta:
                    if line.startswith('>'):
                        seq_name = line.strip().strip('>')
                        seq_name = tig_table.loc[tig_table['tig_long'] == seq_name, 'tig_ext_group']
                        if seq_name.empty:
                            buffer = False
                            continue
                        assert seq_name.size == 1
                        seq_name = seq_name.values[0]
                        assert seq_name not in length_stats
                        buffer = True
                    elif buffer:
                        seq_len = len(line.strip())
                        length_stats[seq_name] += seq_len
                        length_stats['total_bp'] += seq_len
                        chrom_group = seq_name.rsplit('.',1)[-1]
                        length_stats[chrom_group] += seq_len

                        fasta_buffer.write(f'>{seq_name}\n{line}')
                        graph_buffer.write(f'S\t{seq_name}\t{line.strip()}\tLN:i:{seq_len}\n')
                        buffer = False
                        seq_name = None
                    else:
                        continue
        with open(output.fasta, 'w') as dump:
            _ = dump.write(fasta_buffer.getvalue())
        with open(output.gfa, 'w') as dump:
            _ = dump.write(graph_buffer.getvalue())
        with open(output.stats, 'w') as dump:
            for n, c in length_stats.most_common():
                _ = dump.write(f'{n}\t{c}\n')
    # END OF RUN BLOCK

