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

    try:
        gfa = SAMPLE_INFOS[wildcards.sample][wildcards.tigs]
    except KeyError:
        template = 'output/target_assembly/chry_reads/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.r_utg.gfa'
        if 'EC' in wildcards.tigs:
            read_type = 'HIFIEC'
        elif 'AF' in wildcards.tigs:
            read_type = 'HIFIAF'
        else:
            raise
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        assert 'RAW' in wildcards.tigs
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'read_type': read_type,
            'mapq': mapq
        }
        gfa = template.format(**formatter)
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
        'output/read_subsets/chry/{sample_info}_{sample}_{read_type}.chrY-reads.mq60.txt',
        'output/read_subsets/chry/{sample_info}_{sample}_{read_type}.chrY-reads.mq00.txt'
    run:
        import pandas as pd
        with pd.HDFStore(input.cached_aln, 'r') as hdf:
            chry_aln = hdf['chrY']
            highq_reads = chry_aln.loc[chry_aln['mapq'] == 60, 'read_name'].unique().tolist()
            with open(output[0], 'w') as dump:
                _ = dump.write('\n'.join(sorted(highq_reads)) + '\n')
            all_reads = chry_aln['read_name'].unique().tolist()
            with open(output[1], 'w') as dump:
                _ = dump.write('\n'.join(sorted(all_reads)) + '\n')


rule extract_aligned_chry_read_sequences:
    input:
        names = 'output/read_subsets/chry/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type],
    output:
        'output/read_subsets/chry/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.{seq_type}.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 4
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 4
    shell:
        'seqtk subseq {input.reads} {input.names} | pigz -p 4 --best > {output}'
