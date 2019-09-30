
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga.smk'
include: 'racon_polishing.smk'

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
    wildcard_constraints:
        known = 'GRCh38\w+',
        custom = '[\-\w]+'
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
    wildcard_constraints:
        known = 'GRCh38\w+',
        custom = '[\-\w]+'
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.custom_ref}'
        exec += ' &> {log}'
        shell(exec)


rule compute_assembly_delta_reference_canonical_haploid_assembly:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/canonical/{var_callset}/{reference}/{vc_reads}/consensus/{hap_reads}.{hap}.fasta'
    output:
        delta = 'output/evaluation/known_reference/mummer/canonical.{var_callset}.{reference}.{vc_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.delta',
    log:
        'log/output/evaluation/known_reference/mummer/canonical.{var_callset}.{reference}.{vc_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.log',
    benchmark:
        'run/output/evaluation/known_reference/mummer/canonical.{var_callset}.{reference}.{vc_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.rsrc',
    threads: 12
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.haploid_assm} --delta={output.delta} &> {log}'


rule quast_analysis_reference_canonical_haploid_assembly:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/canonical/{var_callset}/{reference}/{vc_reads}/consensus/{hap_reads}.{hap}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/icarus.html'
    threads: 12
    log:
        'log/output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.haploid_assm}'
        exec += ' &> {log}'
        shell(exec)


rule compute_assembly_delta_reference_strandseq_haploid_assembly:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/strandseq/{var_callset}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta'
    output:
        delta = 'output/evaluation/known_reference/mummer/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.delta',
    log:
        'log/output/evaluation/known_reference/mummer/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.log',
    benchmark:
        'run/output/evaluation/known_reference/mummer/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap}.{known}/{known}_vs_{hap_reads}.{hap}.rsrc',
    threads: 12
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.haploid_assm} --delta={output.delta} &> {log}'


rule quast_analysis_reference_strandseq_haploid_assembly:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/strandseq/{var_callset}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/icarus.html'
    threads: 12
    log:
        'log/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.haploid_assm}'
        exec += ' &> {log}'
        shell(exec)


rule quast_analysis_reference_strandseq_polished_haploid_assembly:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        haploid_assm = 'output/diploid_assembly/strandseq/{var_callset}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.{polisher}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/icarus.html'
    threads: 12
    log:
        'log/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/run.rsrc',
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+',
        var_callset = '\w+'
    priority: 100
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.known_ref}'
        exec += ' --features gene:{input.genes}'
        exec += ' --conserved-genes-finding'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.haploid_assm}'
        exec += ' &> {log}'
        shell(exec)