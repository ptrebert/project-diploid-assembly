
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'racon_polishing.smk'
include: 'arrow_polishing.smk'

localrules: master_eval_known_reference

rule master_eval_known_reference:
    input:


rule download_quast_busco_databases:
    output:
        'output/check_files/quast-lg/busco_db_download.ok'
    resources:
        runtime_hrs = 0,
        runtime_min = 30
    shell:
        'quast-download-busco &> {output}'


rule compute_delta_assembly_reference:
    input:
        known_ref = 'references/assemblies/{known_ref}.fasta',
        assembly = 'output/{folder_path}/{file_name}.fasta'
    output:
        delta = 'output/evaluation/mummer_delta/{known_ref}/{folder_path}/{file_name}.delta'
    log:
        'log/output/evaluation/mummer_delta/{known_ref}/{folder_path}/{file_name}.mummer.log'
    benchmark:
        'run/output/evaluation/mummer_delta/{known_ref}/{folder_path}/{file_name}.mummer.rsrc'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((24576 + attempt * 24576) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 24576 + attempt * 24576,
        runtime_hrs = 4
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
        'run/output/evaluation/quastlg_busco/{known_ref}-{genemodel}/{folder_path}/{file_name}/quast_run.rsrc',
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((24576 + attempt * 24576) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 24576 + attempt * 24576,
        runtime_hrs = 4
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    shell:
        'quast-lg.py --threads {threads} -r {input.known_ref}' \
            ' --features gene:{input.genes} --conserved-genes-finding' \
            ' --output-dir {params.output_dir} {input.assembly}' \
            ' &> {log}'
