
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga.smk'
include: 'racon_polishing.smk'

localrules: master_eval_known_reference

rule master_eval_known_reference:
    input:


rule compute_assembly_delta_reference_squashed:
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
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} {input.known_ref} {input.custom_ref} --delta={output.delta} &> {log}'


rule quast_analysis_reference_squashed:
    input:
        known_ref = 'references/assemblies/{known}.fasta',
        custom_ref = 'references/assemblies/{custom}.fasta',
        genes = 'references/downloads/{genemodel}.gff3.gz'
    output:
        pdf_report = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/report.pdf',
        html_icarus = 'output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/icarus.html',
    log:
        'log/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/{custom}.{known}.{genemodel}/run.rsrc',
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
    wildcard_constraints:
        known = 'GRCh38\w+',
        custom = '[\-\w]+'
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    shell:
        'quast-lg.py --threads {threads} -r {input.known_ref}' \
            ' --features gene:{input.genes} --conserved-genes-finding' \
            ' --output-dir {params.output_dir} {input.custom_ref}' \
            ' &> {log}'


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
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
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
    log:
        'log/output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/canonical.{var_callset}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.rsrc',
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    run:
        'quast-lg.py --threads {threads} -r {input.known_ref}' \
            ' --features gene:{input.genes} --conserved-genes-finding' \
            ' --output-dir {params.output_dir} {input.haploid_assm}' \
            ' &> {log}'


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
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
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
    log:
        'log/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/run.rsrc',
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.pdf_report)
    priority: 100
    run:
        'quast-lg.py --threads {threads} -r {input.known_ref}' \
            ' --features gene:{input.genes} --conserved-genes-finding' \
            ' --output-dir {params.output_dir} {input.haploid_assm}' \
            ' &> {log}'


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
    log:
        'log/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/run.log',
    benchmark:
        'run/output/evaluation/known_reference/quastlg_busco/strandseq.{var_callset}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/run.rsrc',
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 65536
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
        'quast-lg.py --threads {threads} -r {input.known_ref}' \
            ' --features gene:{input.genes} --conserved-genes-finding' \
            ' --output-dir {params.output_dir} {input.haploid_assm}' \
            ' &> {log}'