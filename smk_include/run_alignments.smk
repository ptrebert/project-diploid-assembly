
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'canonical_dga.smk'

localrules: master_run_alignments

rule master_run_alignments:
    input:


rule minimap_reads_to_reference_alignment:
    input:
        reads = 'input/fastq/complete/{sample}.fastq.gz',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/alignments/reads_to_reference/{sample}_map-to_{reference}.sam'
    log:
        'log/output/alignments/reads_to_reference/{sample}_map-to_{reference}.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{sample}_map-to_{reference}.rsrc'
    wildcard_constraints:
        reference = '[\w\-]+'
    threads: 64
    run:
        preset = None
        is_pacbio = False
        individual = wildcards.sample.split('_')[0]
        if '_pb' in wildcards.sample:
            is_pacbio = True
            if '-clr' in wildcards.sample:
                preset = 'map-pb'
            elif '-ccs' in wildcards.sample:
                preset = 'asm20'
            else:
                raise ValueError('Unknown minimap2 preset: {}'.format(wildcards.sample))
        if '_ont' in wildcards.sample:
            preset = 'map-ont'

        assert preset is not None, 'No minimap preset selected for sample {}'.format(wildcards.sample)

        exec = 'minimap2 -t {threads} -x ' + preset
        if is_pacbio:
            exec += ' -H '
        exec += ' -a -o {output}'
        exec += ' -R "@RG\\tID:1\\tSM:' + individual + '" '
        exec += ' {input.reference} {input.reads}'
        exec += ' &> {log}'
        shell(exec)


rule minimap_reads_to_contig_aln_polish_round1:
    input:
        reads = 'output/diploid_assembly/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.{hap}.fastq.gz',
        contigs = 'output/diploid_assembly/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/consensus/{readset}.{hap}.fasta'
    output:
        alignments = 'output/diploid_assembly/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/polishing/{polisher}-pass1/{readset}.{hap}.psort.sam',
    log:
        'log/output/diploid_assembly/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/polishing/{polisher}-pass1/{readset}.{hap}.log'
    benchmark:
        'run/output/diploid_assembly/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/polishing/{polisher}-pass1/{readset}.{hap}.run'
    threads: 64
    run:
        preset = None
        is_pacbio = False
        individual = wildcards.readset.split('_')[0]
        if '_pb' in wildcards.readset:
            is_pacbio = True
            if '-clr' in wildcards.readset:
                preset = 'map-pb'
            elif '-ccs' in wildcards.readset:
                preset = 'asm20'
            else:
                raise ValueError('Unknown minimap2 preset: {}'.format(wildcards.readset))
        if '_ont' in wildcards.readset:
            preset = 'map-ont'

        assert preset is not None, 'No minimap preset selected for sample {}'.format(wildcards.readset)

        exec = 'minimap2 -t {threads} -x ' + preset
        if is_pacbio:
            exec += ' -H '
        exec += ' -a'
        exec += ' --eqx -m 5000 --secondary=no'  # parameters recommended by Aaron Wenger
        exec += ' -R "@RG\\tID:1\\tSM:' + individual + '" '
        exec += ' {input.contigs} {input.reads}'
        exec += ' 2> {log}'
        exec += ' | samtools sort'
        exec += ' | samtools view -q 10 -F0x704 /dev/stdin'
        exec += ' > {output}'
        shell(exec)




