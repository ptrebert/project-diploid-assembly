
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'

localrules: master_eval_known_reference

rule master_eval_known_reference:
    input:
        rules.master_preprocess_references.input,
        rules.master_prepare_custom_references.input,

        expand('output/evaluation/known_reference/mummer/{reference}_vs_{squashed}.delta',
                reference=['GRCh38_UCSC_noalt', 'GRCh38_ENSv97_primary'],
                squashed=[s + '_sqa' for s in config['samples']]),

        expand('output/evaluation/known_reference/quastlg_busco/{squashed}.{reference}.{genemodel}/report.pdf',
                reference=['GRCh38_UCSC_noalt', 'GRCh38_ENSv97_primary'],
                squashed=[s + '_sqa' for s in config['samples']],
                genemodel=['GRCh38_GENCODEv31_basic'])


rule compute_assembly_delta_reference_custom:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        custom_ref = 'references/assemblies/{custom}.fasta'
    output:
        delta = 'output/evaluation/known_reference/mummer/{known}_vs_{custom}.delta'
    log: 'log/output/evaluation/known_reference/mummer/{known}_vs_{custom}.log'
    benchmark: 'run/output/evaluation/known_reference/mummer/{known}_vs_{custom}.delta.rsrc'
    threads: 24
    run:
        exec = 'nucmer --maxmatch -l 100 -c 500'
        exec += ' --threads={threads}'
        exec += ' {input.known_ref} {input.custom_ref}'
        exec += ' --delta={output.delta}'
        exec += ' &> {log}'
        shell(exec)


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
    log: 'log/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.log',
    benchmark: 'run/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}'
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'  # this is the BUSCO option
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.custom_ref}'
        exec += ' &> {log}'
        shell(exec)