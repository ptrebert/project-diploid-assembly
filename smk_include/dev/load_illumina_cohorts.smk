
import os
import pandas as pd

localrules: master


def collect_existing_samples(path):
    status_files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.status')]

    existing_samples = []

    for sf in status_files:
        with open(sf, 'r') as dump:
            status = dump.read().strip()
            if status != 'PASS':
                continue
            fastq_file = os.path.basename(sf).split('.')[0]
            fastq_file += '.fastq.gz'
            existing_samples.append(fastq_file)
    return set(existing_samples)
        

def determine_698_cohort_downloads(wildcards):

    JANA_LOCATION = '/gpfs/project/ebler/hgsvc/data/high_cov_samples/fastq/related/'
    METADATA = 'metadata/PRJEB36890.698.ready.tsv'

    if not os.path.isfile(METADATA):
        # just trigger generation of clean metadata table
        return []
    
    existing_samples = collect_existing_samples(JANA_LOCATION)

    md = pd.read_csv(METADATA, sep='\t')

    fastq1 = set(md['fastq1_name'].values) - existing_samples
    fastq2 = set(md['fastq2_name'].values) - existing_samples
    fastq_any = fastq1.union(fastq2)

    missing_fastqs = sorted([os.path.join('output', 'PRJEB36890.698', f.replace('.fastq.gz', '.done')) for f in fastq_any])
    return missing_fastqs


rule master:
    input:
        'metadata/PRJEB36890.698.ready.tsv',
        determine_698_cohort_downloads


rule prepare_metadata:
    input:
        'metadata/{cohort}.tsv'
    output:
        ready = 'metadata/{cohort}.ready.tsv',
        incomplete = 'metadata/{cohort}.incomplete.tsv'
    run:
        import pandas as pd

        keep_columns = [
            'study_accession',
            'sample_accession',
            'sample_alias',
            'run_accession',
            'fastq1_md5',
            'fastq2_md5',
            'fastq1_bytes',
            'fastq2_bytes',
            'fastq1_path',
            'fastq2_path'
        ]

        md = pd.read_csv(input[0], sep='\t')

        extract_md5 = '(?P<fastq1_md5>[a-z0-9]+);(?P<fastq2_md5>[a-z0-9]+)'
        md5_sums = md['fastq_md5'].str.extract(extract_md5, expand=True)

        extract_bytes = '(?P<fastq1_bytes>[0-9]+);(?P<fastq2_bytes>[0-9]+)'
        fastq_sizes = md['fastq_bytes'].str.extract(extract_bytes, expand=True)

        extract_fastq = '(?P<fastq1_path>[A-Za-z0-9/\._]+);(?P<fastq2_path>[A-Za-z0-9/\._]+)'
        fastq_paths = md['fastq_ftp'].str.extract(extract_fastq, expand=True)

        md = pd.concat([md, md5_sums, fastq_sizes, fastq_paths], axis=1)

        cram_based_fastq = md['submitted_ftp'].str.endswith('.cram')
        missing_values = md.isna().any(axis=1)

        drop_out_select = missing_values | ~cram_based_fastq

        drop_outs = md.loc[drop_out_select, :].copy()

        md = md.loc[~drop_out_select, keep_columns].copy()

        # add for simple cross checking with existing files
        md['fastq1_name'] = md['sample_alias'] + '_' + md['fastq1_path'].map(os.path.basename)
        md['fastq2_name'] = md['sample_alias'] + '_' + md['fastq2_path'].map(os.path.basename)

        with open(output.error, 'w') as dump:
            drop_outs.to_csv(dump, sep='\t', header=True, index=False)

        with open(output.ready, 'w') as dump:
            md.to_csv(dump, sep='\t', header=True, index=False)


rule generate_download_request:
    input:
        'metadata/{cohort}.ready.tsv'
    output:
        'requests/{cohort}/{sample}_{accession}_1.request',
        'requests/{cohort}/{sample}_{accession}_2.request',
    run:
        import pandas as pd

        table = pd.read_csv(input[0], sep='\t', header=True)

        subset = table.loc[table['sample_alias'] == wildcards.sample, :]
        assert subset.shape[0] == 1, 'Too many read files for sample: {}'.format(wildcards.sample)

        ebi_ftp_host = 'ftp.sra.ebi.ac.uk'

        for idx in [1,2]:
            with open(output[idx-1], 'w') as request_file:
                remote_path = subset['fastq{}_path'.format(idx)].value.replace(ebi_ftp_host, '')
                assert remote_path.startswith('/'), 'FTP host replacement failed: {} / {}'.format(wildcards.sample, remote_path)
                local_name = '{}_{}_{}.fastq.gz'.format(wildcards.sample, wildcards.accession, idx)
                assert local_name == subset['fastq{}_name'.format(idx)], 'Name mismatch: {} / {}'.format(local_name, subset['fastq{}_name'.format(idx)])
                local_path = os.path.join('output', wildcards.cohort, local_name)
                file_md5 = str(subset['fastq{}_md5'.format(idx)].value)
                file_size = str(subset['fastq{}_bytes'.format(idx)].value)

                _ = request_file.write('\n'.join([
                    remote_path,
                    local_path,
                    file_md5,
                    file_size
                ]) + '\n')


rule run_aspera_download:
    input:
        'requests/{cohort}/{fastq}.request'
    output:
        touch('output/{cohort}/{fastq}.done')
    log:
        'log/output/{fastq}.ascp.log'
    benchmark:
        'rsrc/output{fastq}.ascp.rsrc'
    params:
        remote_path = lambda wildcards, input: open(input[0]).readlines()[0] if os.path.isfile(input[0]) else 'DRY-RUN',
        local_path = lambda wildcards, input: open(input[0]).readlines()[1] if os.path.isfile(input[0]) else 'DRY-RUN',
    shell:
        'ascp -i /home/ebertp/.aspera/connect/etc/asperaweb_id_dsa.openssh '
            '-T -Q -l 200M -P33001 -L- -k2  '
            'era-fasp@fasp.sra.ebi.ac.uk:{params.remote_path}
            '{params.local_path} &> {log}'
