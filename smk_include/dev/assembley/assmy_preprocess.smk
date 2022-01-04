import pathlib as pl

def find_script_path(script_name, subfolder=''):
    """
    Find full path to script to be executed. Function exists
    to avoid config parameter "script_dir"

    :param script_name:
    :param subfolder:
    :return:
    """
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


rule merge_hg002_chry_draft:
    input:
        genome = '/gpfs/project/projects/medbioinf/data/references/T2Tv11_T2TC_chm13.fasta',
        chry = '/gpfs/project/projects/medbioinf/data/references/hg002.chrY.v2.fasta'
    output:
        genome = '/gpfs/project/projects/medbioinf/data/references/T2Tv11_hg002Yv2_chm13.fasta'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import shutil
        import io

        shutil.copy(input.genome, output.genome)

        out_buffer = io.StringIO()
        _ = out_buffer.write('>chrY\n')
        with open(input.chry, 'r') as fasta_in:
            _ = fasta_in.readline()
            _ = out_buffer.write(fasta_in.read())

        with open(output.genome, 'a') as fasta_out:
            _ = fasta_out.write(out_buffer.getvalue())


def select_input_graph(wildcards):
    assemblers = {
        'HAS': 'hifiasm',
        'MBG': 'mbg',
        'LJA': 'lja'
    }
    graphs_paths = {
        'hifiasm': 'output/target_assembly/{chrom}/hifiasm/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.{tigs}.gfa',
        'mbg': 'output/target_assembly/{chrom}/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.gfa',
        'lja': 'no_uncompressed_graph_available'
    }
    try:
        # for hifiasm whole-genome graphs
        gfa = SAMPLE_INFOS[wildcards.sample][wildcards.tigs]
    except KeyError:
        if 'OHEC' in wildcards.tigs:
            read_type = 'OHEC'
        elif 'OEC' in wildcards.tigs:
            read_type = 'ONTEC'
        elif 'EC' in wildcards.tigs:
            read_type = 'HIFIEC'
        elif 'AF' in wildcards.tigs:
            read_type = 'HIFIAF'
        else:
            raise
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise
        assembler = assemblers[wildcards.tigs[:3]]
        if assembler == 'lja':
            raise ValueError(f'Not supported at the moment - assembler LJA / {str(wildcards)}')
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'read_type': read_type,
            'mapq': mapq,
            'chrom': chrom
        }
        if assembler == 'hifiasm':
            formatter['tigs'] = get_hifiasm_tigs(wildcards.tigs)
        elif assembler == 'mbg':
            k, w, r = get_mbg_param(wildcards.tigs)
            formatter['kmer'] = k
            formatter['window'] = w
            formatter['resolvek'] = r
        else:
            raise
        gfa = graphs_paths[assembler].format(**formatter)
    return gfa


rule check_overlong_edges:
    input:
        gfa = select_input_graph
    output:
        discard = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.discard.links'
    conda: '../../../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('gfa_check_ovl.py'),
        stats_out = lambda wildcards, output: output.discard.replace('.discard.links', '.stats')
    resources:
        mem_total_mb = lambda wildcards, attempt: 512 * attempt,
        runtime_hrs = 0,
        runtime_min = lambda wildcards, attempt: 10 * attempt,
    shell:
        '{params.script_exec} --out-discard {output.discard} -g {input.gfa} > {params.stats_out}'


rule clean_input_gfa:
    input:
        gfa = select_input_graph,
        discard = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.discard.links'
    output:
        gfa = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import os
        import io
        import shutil

        skip_lines = []
        with open(input.discard, 'r') as table:
            for line in table:
                if line.startswith('#') or not line.strip():
                    continue
                ln = line.split('\t')[0]
                skip_lines.append(int(ln))
        
        if not skip_lines:
            shutil.copy(input.gfa, output.gfa)
        else:
            gfa_buffer = io.StringIO()
            with open(input.gfa, 'r') as gfa:
                for ln, line in enumerate(gfa, start=1):
                    if ln in skip_lines:
                        continue
                    gfa_buffer.write(line)
            with open(output.gfa, 'w') as gfa:
                _ = gfa.write(gfa_buffer.getvalue())
    # END OF RUN BLOCK


rule merge_reference_confirmed_y_contigs:
    input:
        contigs = lambda wildcards: [str(p) for p in pl.Path(config['path_chry_contigs']).glob('*.fasta')],
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        fa = 'output/references/{reference}.augY.fasta',
        faidx = 'output/references/{reference}.augY.fasta.fai'
    message: 'DEPRECATED'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    shell:
        'cat {input.fa_ref} {input.contigs} > {output.fa}'
            ' && '
            'samtools faidx {output.fa}'


rule extract_confirmed_y_contig_names:
    input:
        contigs = pl.Path(config['path_chry_contigs']) / pl.Path('{sample_long}.TIGPRI_TIGALT_Y_X-PAR_contigs.fasta')
    output:
        'output/references/{sample_long}.XYPAR.contigs.txt'
    shell:
        'egrep "^>" {input.contigs} | egrep -o "[a-z0-9]+$" > {output}'


rule extract_confirmed_y_contig_read_names:
    input:
        tig_names = 'output/references/{sample_info}_{sample}.XYPAR.contigs.txt',
        tig_pri = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['TIGPRI'],
        tig_alt = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['TIGALT'],
    output:
        'output/references/{sample_info}_{sample}.XYPAR.reads.txt'
    message: 'DEPRECATED'
    shell:
        'egrep "^A" {input.tig_pri} | egrep -f {input.tig_names} | cut -f 5 > {output}'
            ' && '
        'egrep "^A" {input.tig_alt} | egrep -f {input.tig_names} | cut -f 5 >> {output}'


rule extract_confirmed_y_contig_reads:
    input:
        read_names = 'output/references/{sample_info}_{sample}.XYPAR.reads.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['HIFIAF'],
    output:
        fq = 'output/references/{sample_info}_{sample}.XYPAR.reads.fastq.gz'
    message: 'DEPRECATED'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt ** 3
    shell:
        'seqtk subseq {input.reads} {input.read_names} | pigz -p 4 --best > {output}'


def select_alignment_cache_file(wildcards):

    cache_dir = pl.Path(config['path_alignment_cache'])
    cache_files = list(cache_dir.glob(f'{wildcards.sample}*{wildcards.read_type}*T2Tv11_hg002Yv2_chm13.cov.cache.h5'))
    if not cache_files:
        raise ValueError('too few ' + str(wildcards))
    if len(cache_files) > 1:
        raise ValueError('too many ' + str(wildcards))
    return cache_files[0]


rule extract_aligned_chry_read_names:
    input:
        cached_aln = select_alignment_cache_file
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.mq60.txt',
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.mq00.txt'
    run:
        import pandas as pd
        with pd.HDFStore(input.cached_aln, 'r') as hdf:
            chrom_aln = hdf[wildcards.chrom]
            highq_reads = chrom_aln.loc[chrom_aln['mapq'] == 60, 'read_name'].unique().tolist()
            with open(output[0], 'w') as dump:
                _ = dump.write('\n'.join(sorted(highq_reads)) + '\n')
            all_reads = chrom_aln['read_name'].unique().tolist()
            with open(output[1], 'w') as dump:
                _ = dump.write('\n'.join(sorted(all_reads)) + '\n')


rule extract_aligned_chrom_read_sequences:
    input:
        names = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type],
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.{seq_type}.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = CONSTRAINT_REGULAR_SAMPLES,
        read_type = '(HIFIEC|HIFIAF|ONTUL|ONTEC)',
        chrom = '(chrX|chrY)'
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 4
    shell:
        'seqtk subseq {input.reads} {input.names} | pigz -p 4 --best > {output}'


rule merge_sex_chrom_reads:
    input:
        chrx = 'output/read_subsets/chrX/{sample_info}_{sample}_{read_type}.chrX-reads.{mapq}.{seq_type}.gz',
        chry = 'output/read_subsets/chrY/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.{seq_type}.gz',
    output:
        'output/read_subsets/chrXY/{sample_info}_{sample}_{read_type}.chrXY-reads.{mapq}.{seq_type}.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF|ONTUL|ONTEC)',
        sample = CONSTRAINT_REGULAR_SAMPLES,
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'pigz -d -c {input.chrx} {input.chry} | pigz --best -p 4 > {output}'


rule merge_read_types:
    input:
        hifiec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_HIFIEC.{chrom}-reads.{mapq}.{seq_type}.gz',
        ontec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTEC.{chrom}-reads.{mapq}.{seq_type}.gz',
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_OHEC.{chrom}-reads.{mapq}.{seq_type}.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        chrom = '(chrY|chrX|chrXY)',
        sample = CONSTRAINT_REGULAR_SAMPLES,
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'pigz -d -c {input.hifiec} {input.ontec} | pigz --best -p 4 > {output}'


def select_afr_mix_subsets(wildcards):

    assert 'mq0' in wildcards.mapq
    if wildcards.read_type == 'HIFIAF':
        seq_type = 'fastq'
    else:
        seq_type = 'fasta'
    merge_samples = SAMPLE_INFOS[wildcards.sample]['merge']
    template = 'output/read_subsets/{chrom}/{sample_long}_{read_type}.{chrom}-reads.{mapq}.{seq_type}.gz'
    selected_readsets = []
    for sample in merge_samples:
        sample_long = SAMPLE_INFOS[sample]['long_id']
        formatter = {
            'sample_long': sample_long,
            'read_type': wildcards.read_type,
            'mapq': wildcards.mapq,
            'seq_type': seq_type,
            'chrom': wildcards.chrom
        }
        selected_readsets.append(template.format(**formatter))
    return sorted(selected_readsets)


rule merge_afr_mix_subsets:
    input:
        reads = select_afr_mix_subsets
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.{seq_type}.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(NA193N7|NA193NN|AFR4MIX)'
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'pigz -d -c {input.reads} | pigz --best -p 4 > {output}'
