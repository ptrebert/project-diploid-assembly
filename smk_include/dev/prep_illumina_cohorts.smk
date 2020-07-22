
import os
import pandas as pd

localrules: master,
            check_file_consistency,
            delete_original_fastq


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


def determine_conversion_targets(wildcards):

    import collections as col

    input_paths = [
        '/gpfs/project/ebertp/data/share/jana/run_folder/output/PRJEB36890.698'
    ]

    file_pairs = col.defaultdict(list)

    for p in input_paths:
        cohort = os.path.split(p)[-1]
        done_files = [f for f in os.listdir(p) if f.endswith('done')]
        [file_pairs[(cohort, f.rsplit('_', 1)[0])].append(f) for f in done_files]


    # 'output/reduced/{cohort}/{sample}_{accession}.fastq.del'
    conv_targets = set()
    for (cohort, sample_id), done_files in file_pairs.items():
        if len(done_files) != 2:
            continue
        #conv_targets.add(os.path.join('output', 'reduced', cohort, sample_id + '.fastq.del'))
        conv_targets.add(os.path.join('output', 'reduced', cohort, sample_id + '.fasta.count'))

    return sorted(conv_targets)


rule master:
    input:
        determine_conversion_targets


rule check_file_consistency:
    input:
        req1 = 'requests/{cohort}/{sample}_{accession}_1.request',
        req2 = 'requests/{cohort}/{sample}_{accession}_2.request',
        done1 = 'output/{cohort}/{sample}_{accession}_1.done',
        done2 = 'output/{cohort}/{sample}_{accession}_2.done',
    output:
        touch('output/{cohort}/{sample}_{accession}.check.ok')
    log:
        'log/output/{cohort}/{sample}_{accession}.check.log'
    run:
        import os
        import time
        import io

        error_buffer = io.StringIO()
        log_buffer = io.StringIO()

        log_buffer.write('{}\n'.format(time.ctime()))

        for req_file_path in [input.req1, input.req2]:

            with open(req_file_path, 'r') as reqfile:
                _ = reqfile.readline()
                local_fastq = reqfile.readline().strip()
                fastq_md5 = reqfile.readline().strip()
                fastq_bytes = int(reqfile.readline().strip())
                if not os.path.isfile(local_fastq):
                    error_buffer.write('Invalid path to local FASTQ: {}\n'.format(local_fastq))
                    log_buffer.write('Error: invalid path to local FASTQ: {}\n'.format(local_fastq))
                else:
                    statinfo = os.stat(local_fastq)
                    size_on_fs = int(statinfo.st_size)
                    if size_on_fs != fastq_bytes:
                        error_buffer.write('Error: local FASTQ "{}" has unexpected size: {} vs {}\n'.format(
                            local_fastq, fastq_bytes, size_on_fs)
                            )
                        log_buffer.write('Error: local FASTQ "{}" has unexpected size: {} vs {}\n'.format(
                            local_fastq, fastq_bytes, size_on_fs)
                            )
                    else:
                        log_buffer.write('File check ok: {}\n'.format(local_fastq))
        
        with open(log[0], 'w') as log_dump:
            _ = log_dump.write(log_buffer.getvalue())
        
        errmsg = error_buffer.getvalue()
        if len(errmsg) > 1:
            raise ValueError(errmsg)


rule convert_fastq_to_gzip_fasta:
    input:
        done1 = 'output/{cohort}/{sample}_{accession}_1.done',
        done2 = 'output/{cohort}/{sample}_{accession}_2.done',
        check_ok = 'output/{cohort}/{sample}_{accession}.check.ok'
    output:
        check = touch('output/reduced/{cohort}/{sample}_{accession}.convert.ok')
    log:
        'log/output/reduced/{cohort}/{sample}_{accession}.convert.log'
    benchmark:
        'rsrc/output/reduced/{cohort}/{sample}_{accession}.convert.rsrc'
    conda:
        '../../environment/conda/conda_convert.yml'
    threads: 2
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((4096 * attempt) / 2),
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        fastq1 = lambda wildcards, input: input.done1.replace('.done', '.fastq.gz'),
        fastq2 = lambda wildcards, input: input.done2.replace('.done', '.fastq.gz'),
        gzip_fasta = lambda wildcards, output: output.check.replace('.convert.ok', '.fasta.gz')
    shell:
        '( ( seqtk seq -A -C {params.fastq1} ; seqtk seq -A -C {params.fastq2} ; ) | gzip > {params.gzip_fasta} ; ) 2> {log}'


rule count_reads_in_fasta:
    input:
        done1 = 'output/{cohort}/{sample}_{accession}_1.done',
        done2 = 'output/{cohort}/{sample}_{accession}_2.done',
        convert_ok = 'output/reduced/{cohort}/{sample}_{accession}.convert.ok',
    output:
        'output/reduced/{cohort}/{sample}_{accession}.fasta.count'
    params:
        gzip_fasta = lambda wildcards, input: input.convert_ok.replace('.convert.ok', '.fasta.gz')
    shell:
        'zgrep -E -c "^>" {params.gzip_fasta} > {output}'


rule delete_original_fastq:
    input:
        metadata = 'metadata/PRJEB36890.698.ready.tsv',
        done1 = 'output/{cohort}/{sample}_{accession}_1.done',
        done2 = 'output/{cohort}/{sample}_{accession}_2.done',
        convert_ok = 'output/reduced/{cohort}/{sample}_{accession}.convert.ok',
        read_count = 'output/reduced/{cohort}/{sample}_{accession}.fasta.count'
    output:
        'output/reduced/{cohort}/{sample}_{accession}.fastq.del'
    run:
        import os
        import sys
        import pandas as pd

        fasta_path = input.convert_ok.replace('.convert.ok', '.fasta.gz')
        if not os.path.isfile(fasta_path):
            raise ValueError('FASTA conversion failed / no file at: {}'.format(fasta_path))
        
        df = pd.read_csv(input.metadata, sep='\t')
        expected_count = int(df.loc[df['run_accession'] == wildcards.accession, 'read_count'])
        assert expected_count > 0, 'Expected read count is zero: {} / {}'.format(wildcards.sample, wildcards.accession)

        with open(input.read_count, 'r') as dump:
            observed_count = int(dump.read())
        
        if observed_count != expected_count:
            error_msg = '\nRead count mismatch: {} / Obs. {} / Exp. {}\n'.format(
                wildcards.sample, observed_count, expected_count
                )
            sys.stderr.write(error_msg)
            raise RuntimeError(error_msg)

        with open(output[0], 'w') as out:
            fastq1_path = input.done1.replace('.done', '.fastq.gz')
            if os.path.isfile(fastq1_path):
                out.write('FASTQ will be deleted: {}\n'.format(fastq1_path))
                os.remove(fastq1_path)
            else:
                out.write('FASTQ no longer exists: {}\n'.format(fastq1_path))
            
            fastq2_path = input.done2.replace('.done', '.fastq.gz')
            if os.path.isfile(fastq2_path):
                out.write('FASTQ will be deleted: {}\n'.format(fastq2_path))
                os.remove(fastq2_path)
            else:
                out.write('FASTQ no longer exists: {}\n'.format(fastq2_path))
