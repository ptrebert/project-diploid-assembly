
include: 'preprocess_input.smk'

localrules: master_prepare_custom_references

rule master_prepare_custom_references:
    input:


rule ex_nihilo_custom_reference_injections:
    output:
        protected('references/assemblies/HG00733_sra_pbsq1-clr_clustV1.fasta')


rule compute_wtdbg_squashed_assembly_layout:
    input:
        fastq = 'input/fastq/complete/{sample}_1000.fastq.gz'
    output:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
    log: 'log/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.log',
    benchmark: 'run/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.rsrc',
    threads: 48
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.sample]
    run:
        exec = 'wtdbg2 -x {params.param_preset}'
        exec += ' -i {input.fastq}'
        exec += ' -g3g -t {threads}'
        exec += ' -o references/assemblies/squashed_layout/wtdbg2/{wildcards.sample}/{wildcards.sample}'
        exec += ' &> {log}'
        shell(exec)


rule compute_wtdbg_squashed_assembly_consensus:
    input:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz'
    output:
        squashed_assembly = 'references/assemblies/{sample}_sqa.fasta'
    log: 'log/references/assemblies/{sample}_sqa.consensus.log'
    benchmark: 'run/references/assemblies/{sample}_sqa.consensus.rsrc'
    threads: 48
    run:
        exec = 'wtpoa-cns -t {threads}'
        exec += ' -i {input.layout}'
        exec += ' -o {output.squashed_assembly}'
        exec += ' &> {log}'
        shell(exec)


rule filter_squashed_assembly_by_size:
    input:
        'references/assemblies/{sample}_{assembly_type}.fasta'
    output:
        fasta = 'references/assemblies/{sample}_{assembly_type}-100kb.fasta',
        stats = 'output/statistics/assemblies/{sample}_{assembly_type}-100kb.stats.tsv'
    log:
        'log/references/assemblies/{sample}_{assembly_type}-100kb.log'
    wildcard_constraints:
        assembly_type = '(sqa|clustV[0-9])'
    params:
        scriptdir = config['script_dir']
    run:
        exec = '{params.scriptdir}/filter_squashed_assembly.py --debug'
        exec += ' --input-fasta {input}'
        exec += ' --output-fasta {output.fasta}'
        exec += ' --output-metrics {output.stats}'
        exec += ' --min-size 100000 &> {log}'
        shell(exec)