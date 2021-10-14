
rule create_conda_environment_ashleys:
    output:
        'output/check_files/environment/conda_ashleys.ok'
    log:
        'log/output/check_files/environment/conda_ashleys.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_ashleys.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule install_source_ashleys_qc:
    input:
        'output/check_files/environment/conda_ashleys.ok'
    output:
        touch('output/check_files/src_build/install_ashleys.ok')
    log:
       'log/output/check_files/src_build/install_ashleys.log'
    conda:
        '../environment/conda/conda_ashleys.yml'
    params:
        repo_folder = 'output/repositories'
    shell:
         '( rm -rf {params.repo_folder}/ashleys-qc && '
         'mkdir -p {params.repo_folder} && '
         'cd {params.repo_folder} && '
         'git clone https://github.com/friendsofstrandseq/ashleys-qc.git && '
         'cd ashleys-qc && '
         'git checkout develop && '
         'python setup.py install ; ) > {log} 2>&1'


def select_raw_strandseq_input_mate1(wildcards):
    return _select_raw_strandseq_input_mate(wildcards.sseq_reads, wildcards.lib_id, 0)


def select_raw_strandseq_input_mate2(wildcards):
    return _select_raw_strandseq_input_mate(wildcards.sseq_reads, wildcards.lib_id, 1)


def _select_raw_strandseq_input_mate(sseq_reads, lib_id, mate_num):
    import os
    import fnmatch
    import json

    source_file = os.path.join('input', 'data_sources', sseq_reads + '.json')
    try:
        with open(source_file, 'r') as source_info:
            dump = json.load(source_info)
        raw_input_files = [(d['local_path'], d['remote_path']) for d in dump.values()]
    except FileNotFoundError:
        sseq_data_source = read_local_data_sources()[sseq_reads]
        raw_input_files = [(d['local_path'], d['remote_path']) for d in sseq_data_source.values()]

    matched_input_files = sorted(filter(lambda x: lib_id in x[0], raw_input_files))
    if len(matched_input_files) != 2:
        raise ValueError('Could not extract two mates for Strand-seq input: {} / {} / {}'.format(sseq_reads, lib_id, matched_input_files))

    # return the remote path for the PGAS-internal Strand-seq library ID
    mate_file = matched_input_files[mate_num][1]
    return mate_file


SSEQ_QC_REFERENCE = 'no-reference' if 'sseq_qc_reference' not in config else config['sseq_qc_reference'] 

ruleorder: run_strandseq_qc_alignment > samtools_index_bam_alignment

rule run_strandseq_qc_alignment:
    input:
        data_source = 'input/data_sources/{sseq_reads}.json',
        ref_idx = 'references/assemblies/{}/bwa_index/{}.bwt'.format(SSEQ_QC_REFERENCE, SSEQ_QC_REFERENCE),
        mate1 = select_raw_strandseq_input_mate1,
        mate2 = select_raw_strandseq_input_mate2,
    output:
        bam_sort = 'temp/ashleys_qc/{sseq_reads}/{lib_id}.psort.bam',
        bam_mdup = 'temp/ashleys_qc/{sseq_reads}/{lib_id}.psort.mdup.bam',
        bai = 'temp/ashleys_qc/{sseq_reads}/{lib_id}.psort.mdup.bam.bai',
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(12288 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt
    params:
        idx_prefix = lambda wildcards, input: input.ref_idx.rsplit('.', 1)[0],
        individual = lambda wildcards: wildcards.sseq_reads.split('_')[0]
    shell:
        'bwa mem -t {threads}'
            ' -R "@RG\\tID:{wildcards.sseq_reads}_{wildcards.lib_id}\\tPL:Illumina\\tSM:{params.individual}"'
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} | '
            ' samtools sort -@ 2 -m 1024M -o {output.bam_sort} -O BAM /dev/stdin '
            ' && '
            ' sambamba markdup -t {threads} --overflow-list-size 600000 {output.bam_sort} {output.bam_mdup} '
            ' && '
            ' samtools index {output.bam_mdup}'


def select_raw_strandseq_input_bam(wildcards):
    return _select_raw_strandseq_input_all(wildcards.sseq_reads, '')

def select_raw_strandseq_input_bai(wildcards):
    return _select_raw_strandseq_input_all(wildcards.sseq_reads, '.bai')


def _select_raw_strandseq_input_all(sseq_reads, file_ext):
    import os
    import json

    source_file = os.path.join('input', 'data_sources', sseq_reads + '.json')
    try:
        with open(source_file, 'r') as source_info:
            dump = json.load(source_info)
        raw_input_files = [(d['local_path'], d['remote_path']) for d in dump.values()]
    except FileNotFoundError:
        sseq_data_source = read_local_data_sources()[sseq_reads]
        raw_input_files = [(d['local_path'], d['remote_path']) for d in sseq_data_source.values()]

    all_libs = set([os.path.basename(t[0]).rsplit('_', 1)[0] for t in raw_input_files])

    output_files = sorted(['temp/ashleys_qc/{}/{}.psort.mdup.bam{}'.format(sseq_reads, l, file_ext) for l in all_libs])

    return output_files


rule compute_strandseq_library_features:
    input:
        setup_ok = 'output/check_files/src_build/install_ashleys.ok',
        bam_files = select_raw_strandseq_input_bam,
        bai_files = select_raw_strandseq_input_bai,
    output:
        table = 'temp/ashleys_qc/{sseq_reads}.features.tsv'
    log:
        'log/temp/ashleys_qc/{sseq_reads}.features.log'
    benchmark:
        'rsrc/temp/ashleys_qc/{{sseq_reads}}.features.t{}.rsrc'.format(config['num_cpu_medium'])
    conda:
        '../environment/conda/conda_ashleys.yml'
    threads: config['num_cpu_medium']
    params:
        bam_folder = lambda wildcards: os.path.join('temp/ashleys_qc', wildcards.sseq_reads),
        bam_ext = '.psort.mdup.bam',
        windows = '5000000 2000000 1000000 800000 600000 400000 200000',
    shell:
        'ashleys.py -j {threads} --verbose features -f {params.bam_folder} '
        ' -w {params.windows} --recursive_collect -e {params.bam_ext} -o {output.table} &> {log}'


rule predict_strandseq_library_quality:
    input:
        setup_ok = 'output/check_files/src_build/install_ashleys.ok',
        features = 'temp/ashleys_qc/{sseq_reads}.features.tsv'
    output:
        table = 'temp/ashleys_qc/{sseq_reads}/prediction.tsv',
        subset = temp('temp/ashleys_qc/{sseq_reads}/prediction_critical.tsv'),
    conda:
        '../environment/conda/conda_ashleys.yml'
    params:
        model_file = os.path.abspath('output/repositories/ashleys-qc/models/svc_default.pkl')
    shell:
        'ashleys.py predict -p {input.features} -o {output.table} -m {params.model_file}'


def determine_suffix_to_trim(file_names):
    """
    Guess what suffix can be trimmed off of all file names
    Faster with itt.takewhile / zip on strings?
    """
    import difflib
    import collections as col

    rev_first = file_names[0][::-1]
    seqmatch = difflib.SequenceMatcher()
    seqmatch.set_seq2(rev_first)

    suffix_counts = col.Counter()

    for fn in file_names[1:]:
        seqmatch.set_seq1(fn[::-1])
        for block in seqmatch.get_matching_blocks():
            assert block.a == 0, 'Block match does not start at beginning: {} / {} / {}'.format(rev_first, fn[::-1], block)
            suffix_counts[rev_first[:block.size][::-1]] += 1
            break
    (top_suffix, top_count), (next_suffix, next_count) = suffix_counts.most_common(2)
    if top_count == next_count:
        raise ValueError('Cannot estimate suffix to trim from file names: '
                        '{}/#{} vs {}/#{}'.format(top_suffix, top_count, next_suffix, next_count))

    return top_suffix


def clean_table_index(index_names):
    """
    """
    trim_suffix = determine_suffix_to_trim(index_names)
    new_index = []
    for name in index_names:
        # in Python 3.9: str.removesuffix  :-)
        # here: make sure only suffix is removed
        # no matter how weird the file name is...
        new_name = name.rsplit(trim_suffix, 1)[0]
        new_index.append(new_name)
    return new_index


def find_source_library_ids(data_source_file, library_ids):
    """
    """
    import os
    import json
    import pandas as pd
    with open(data_source_file, 'r') as dump:
        metadata = json.load(dump)

    # select only first mate
    matched_names = dict()
    for lib_path, path_info in metadata.items():
        if lib_path.endswith('_2'):
            continue
        lib_id = os.path.basename(lib_path).rsplit('_1', 1)[0]
        remote_name = os.path.basename(path_info['remote_path'])
        matched_names[lib_id] = remote_name

    num_intersect = len(set(matched_names.keys()).intersection(set(library_ids)))
    num_extract = len(matched_names)

    if num_intersect != num_extract:
        raise ValueError('PGAS library ID mismatch: {} vs {}'.format(num_intersect, num_extract))

    trim_suffix = determine_suffix_to_trim(sorted(matched_names.values()))
    source_ids = []
    for lib_id in library_ids:
        matched_name = matched_names[lib_id]
        matched_name = matched_name.rsplit(trim_suffix, 1)[0]
        source_ids.append(matched_name)
    
    source_ids = pd.Series(source_ids, index=pd.Index(library_ids), name='library_source')
    return source_ids
      

rule exclude_low_quality_libraries:
    input:
        data_source = 'input/data_sources/{sseq_reads}.json',
        quality_labels = 'temp/ashleys_qc/{sseq_reads}/prediction.tsv',
        library_features = 'temp/ashleys_qc/{sseq_reads}.features.tsv'
    output:
        table = 'output/sseq_qc/{sseq_reads}.ashleys-qc.tsv',
        listing = 'output/sseq_qc/{sseq_reads}.exclude.txt',
    log:
        'log/output/sseq_qc/{sseq_reads}.ashleys-qc.log'
    params:
        t_conf_low = 0.3,  # follows hard-coded values in ASHLEYS
        t_conf_high = 0.7,
        error_lowq = config.get('sseq_qc_fail_threshold', 50)
    run:
        import collections as col
        import pandas as pd

        assert 1 < params.error_lowq, f'Invalid value [1...+INF]: {params.error_lowq}'

        with open(log[0], 'w') as logfile:
            _ = logfile.write(f'strandseq_sample\t{wildcards.sseq_reads}\n')
    
            qc_labels = pd.read_csv(
                input.quality_labels,
                sep='\t',
                header=0,
                # library here: PGAS-internal ID
                names=['library_id', 'label', 'probability'],
                index_col='library_id'
            )
            total_num_libs = qc_labels.shape[0]
            _ = logfile.write(f'total_num_libs\t{total_num_libs}\n')

            qc_features = pd.read_csv(
                input.library_features,
                sep='\t',
                header=0,
                index_col='sample_name'
            )
            num_feat_libs = qc_features.shape[0]
            _ = logfile.write(f'total_num_features\t{qc_features.shape[1]}\n')
            if num_feat_libs != total_num_libs:
                _ = logfile.write(f'ERROR - number of libraries in feature table is inconsistent: {num_feat_libs} / {total_num_libs}\n')
                raise ValueError(f'Cannot process ASHLEYS output for SSEQ reads {wildcards.sseq_reads}')

            qc_info = pd.concat([qc_labels, qc_features], axis=1, ignore_index=False)
            qc_info['label_confidence'] = 'high'
            selector = (qc_info['probability'] > params.t_conf_low) & (qc_info['probability'] < params.t_conf_high)
            qc_info.loc[selector, 'label_confidence'] = 'low'

            _ = logfile.write(f'LowConfBracket_lower_bound\t{params.t_conf_low}\n')
            _ = logfile.write(f'LowConfBracket_upper_bound\t{params.t_conf_high}\n')

            num_low_conf = (qc_info['label_confidence'] == 'low').sum()
            _ = logfile.write(f'total_num_libs_LCB\t{num_low_conf}\n')
            pct_low_conf = round(num_low_conf / total_num_libs * 100, 1)
            _ = logfile.write(f'total_pct_libs_LCB\t{pct_low_conf}\n')

            qc_info.sort_index(inplace=True)
            new_index_names = clean_table_index(qc_info.index.tolist())
            qc_info.index = pd.Index(new_index_names, name='library_id')
            _ = logfile.write('reset_index_success\tyes\n')

            # add info of original library name
            source_ids = find_source_library_ids(input.data_source, new_index_names)
            qc_info = pd.concat([qc_info, source_ids], axis=1, ignore_index=False)
            _ = logfile.write('library_source_added\tyes\n')

            column_sort_order = ['library_source', 'label', 'probability', 'label_confidence'] + qc_features.columns.tolist()
            qc_info = qc_info[column_sort_order]

            lowq_libraries = (qc_info['label'] == 0).sum()
            usable_libraries = total_num_libs - lowq_libraries
            _ = logfile.write(f'fail_threshold_lowQ\t{params.error_lowq}\n')
            _ = logfile.write(f'total_lowQ_libs\t{lowq_libraries}\n')
            _ = logfile.write(f'total_usable_libs\t{usable_libraries}\n')
            frac_lowq_libs = round(lowq_libraries / total_num_libs, 2)
            pct_lowq_libs = round(lowq_libraries / total_num_libs * 100, 1)
            _ = logfile.write(f'fraction_lowQ_libs\t{frac_lowq_libs}\n')
            _ = logfile.write(f'pct_lowQ_libs\t{pct_lowq_libs}\n')

            if usable_libraries < params.error_lowq:
                _ = logfile.write('ERROR\tnot_enough_usable_libraries\n')
                # create an untracked error output file for simpler
                # summary of what samples failed and why
                with open(f'ERROR_sseq-quality_{wildcards.sseq_reads}.err', 'w'):
                    pass
                ignore_fail = bool(config.get('sseq_qc_ignore_fail', False))
                if not ignore_fail:
                    raise ValueError(f'Cannot process ASHLEYS output for SSEQ reads (lowQ) {wildcards.sseq_reads}')

            qc_info.to_csv(
                output.table,
                sep='\t',
                header=True,
                index=True,
                index_label='library_id'
            )
            _ = logfile.write(f'output_table_saved\t{output.table}\n')
            
            qc_info.loc[qc_info['label'] == 0, 'library_source'].to_csv(
                output.listing,
                sep='\t',
                header=False,
                index=False
            )
            _ = logfile.write(f'output_list_saved\t{output.listing}\n')
            _ = logfile.write('operation_success\tyes\n')
    # END OF RUN BLOCK
