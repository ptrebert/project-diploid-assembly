
localrules: master_eval_known_reference

rule master_eval_known_reference:
    input:
        []


rule compute_delta_assembly_reference:
    input:
        known_ref = 'references/assemblies/{known_ref}.fasta',
        assembly = 'output/{folder_path}/{file_name}.fasta'
    output:
        delta = 'output/evaluation/mummer_delta/{known_ref}/{folder_path}/{file_name}.delta'
    log:
        'log/output/evaluation/mummer_delta/{known_ref}/{folder_path}/{file_name}.mummer.log'
    benchmark:
        'run/output/evaluation/mummer_delta/{{known_ref}}/{{folder_path}}/{{file_name}}.mummer.t{}.rsrc'.format(config['num_cpu_medium'])
    conda:
         '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((24576 + attempt * 24576) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 24576 + attempt * 24576,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.custom_ref} --delta={output.delta} &> {log}'


rule quast_analysis_assembly:
    input:
        dl_chk = 'output/check_files/quast-lg/busco_db_download.ok',
        known_ref = 'references/assemblies/{known_ref}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz',
        assembly = 'output/{folder_path}/{file_name}.fasta',
    output:
        pdf_report = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/{folder_path}/{file_name}/report.pdf',
        html_icarus = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/{folder_path}/{file_name}/icarus.html',
    log:
        'log/output/evaluation/quastlg_busco/{known_ref}-{genemodel}/{folder_path}/{file_name}/quast_run.log',
    benchmark:
        'run/output/evaluation/quastlg_busco/{{known_ref}}-{{genemodel}}/{{folder_path}}/{{file_name}}/quast_run.t{}.rsrc'.format(config['num_cpu_medium'])
    conda:
         '../environment/conda/conda_rtools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((36864 + attempt * 36864) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 36864 + attempt * 36864,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    shell:
        'quast-lg.py --threads {threads} -r {input.known_ref}'
            ' --features gene:{input.genes} --conserved-genes-finding'
            ' --output-dir {params.output_dir} {input.assembly}'
            ' &> {log}'
