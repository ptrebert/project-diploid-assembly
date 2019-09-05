
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_eval_known_reference

rule master_eval_known_reference:
    input:


rule compute_assembly_delta_reference_custom:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        custom_ref = 'references/assemblies/{custom}.fasta'
    output:
        delta = 'output/evaluation/known_reference/mummer/{known}_vs_{custom}.delta'
    log:
        'log/output/evaluation/known_reference/mummer/{known}_vs_{custom}.log'
    benchmark:
        'run/output/evaluation/known_reference/mummer/{known}_vs_{custom}.delta.rsrc'
    threads: 12
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.custom_ref} --delta={output.delta} &> {log}'


rule gzip_delta_for_assemblytics:
    input:
        'output/evaluation/known_reference/mummer/{ref1}_vs_{ref2}.delta'
    output:
        'output/evaluation/known_reference/assemblytics/delta_upload/{ref1}_vs_{ref2}.delta.gz'
    shell:
        "gzip -c {input} > {output}"


rule quast_analysis_reference_custom:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        custom_ref = 'references/assemblies/{custom}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/icarus.html',
    threads: 12
    log:
        'log/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}'
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.custom_ref}'
        exec += ' &> {log}'
        shell(exec)


rule compute_assembly_delta_reference_haploid_assembly:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/{approach}/{var_callset}/{reference}/{readset}/consensus/{readset}.{hap}.fasta'
    output:
        delta = 'output/evaluation/known_reference/mummer/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}/{known}_vs_{readset}.{hap}.delta',
    log:
        'log/output/evaluation/known_reference/mummer/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}/{known}_vs_{readset}.{hap}.log',
    benchmark:
        'run/output/evaluation/known_reference/mummer/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}/{known}_vs_{readset}.{hap}.rsrc',
    threads: 12
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.haploid_assm} --delta={output.delta} &> {log}'


rule quast_analysis_reference_haploid_assembly:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/{approach}/{var_callset}/{reference}/{readset}/consensus/{readset}.{hap}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}.{genemodel}/icarus.html'
    threads: 12
    log:
        'log/output/evaluation/known_reference/quastlg_busco/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/{approach}.{var_callset}.{reference}.{readset}.{hap}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.haploid_assm}'
        exec += ' &> {log}'
        shell(exec)